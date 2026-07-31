#!/usr/bin/env julia

"""Blocked KPM Hartree--Fock SCF on a power-of-two open rectangle."""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 1:2 || error(
    "usage: run_rectangular_kpm_scf.jl OUTPUT_DIRECTORY [BACKEND]",
)
output = abspath(ARGS[1])
backend = length(ARGS) == 2 ? Symbol(lowercase(ARGS[2])) : :cpu
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")
if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

const NX = parse(Int, get(ENV, "KPM_RECT_NX", "1024"))
const NY = parse(Int, get(ENV, "KPM_RECT_NY", "512"))
ispow2(NX) && ispow2(NY) ||
    error("KPM_RECT_NX and KPM_RECT_NY must be powers of two")
const N = NX * NY
const NE = div(N, 2)
const X_BITS = trailing_zeros(NX)
const Y_BITS = trailing_zeros(NY)
const V2 = 0.5
const U = 1.0
const TAU = sqrt(2.0) - 5.0 / 6.0
const SEED_AMPLITUDE = 0.1
const MOMENTS = 1200
const PROBES = 1024
const PROBE_BLOCK = 128
const PROBE_SEED = 510578
const FINAL_MOMENTS = 1600
const FINAL_PROBES = 2048
const FINAL_SEED = 20260730
const SCF_MIXING = 0.5
const SCF_TOLERANCE = 1e-3
const SCF_MAX_ITERATIONS = 30
const SCF_MIXER = Symbol(lowercase(get(ENV, "KPM_SCF_MIXER", "linear")))
SCF_MIXER in (:linear, :pulay) ||
    error("KPM_SCF_MIXER must be linear or pulay")
const PULAY_HISTORY = parse(Int, get(ENV, "KPM_SCF_PULAY_HISTORY", "6"))
const PULAY_WARMUP = parse(Int, get(ENV, "KPM_SCF_PULAY_WARMUP", "4"))
const PULAY_REGULARIZATION = parse(
    Float64, get(ENV, "KPM_SCF_PULAY_REGULARIZATION", "1e-12"),
)
const PULAY_DAMPING = parse(
    Float64, get(ENV, "KPM_SCF_PULAY_DAMPING", string(SCF_MIXING)),
)
0 < PULAY_DAMPING <= 1 ||
    error("KPM_SCF_PULAY_DAMPING must lie in (0, 1]")
const TRACE_TOLERANCE = 1e-6 * NE
const SPECTRAL_MARGIN = 0.1

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
site_index(x, y) = x + NX * y + 1
tx(x) = -1.0 - V2 * cos(2π * TAU * (Float64(x) + 0.5))
ty(y) = -1.0 - V2 * cos(2π * TAU * (Float64(y) + 0.5))
seed_field(x, y) = iseven(x + y) ? SEED_AMPLITUDE : -SEED_AMPLITUDE
relative_change(candidate, reference) =
    norm(candidate - reference) / max(norm(reference), sqrt(eps(Float64)))

coordinate_code(x, y) =
    coordinate_interleaved_code(x, y, X_BITS, Y_BITS)

function lattice_data()
    codes = Vector{Int}(undef, N)
    seed = Vector{Float64}(undef, N)
    bonds = Tuple{Int,Int,Symbol}[]
    hopping = Float64[]
    for y in 0:(NY - 1), x in 0:(NX - 1)
        site = site_index(x, y)
        codes[site] = coordinate_code(x, y)
        seed[site] = seed_field(x, y)
        if x < NX - 1
            push!(bonds, (site, site_index(x + 1, y), :horizontal))
            push!(hopping, tx(x))
        end
        if y < NY - 1
            push!(bonds, (site, site_index(x, y + 1), :vertical))
            push!(hopping, ty(y))
        end
    end
    sort(codes) == collect(0:(N - 1)) ||
        error("coordinate codes are not a permutation")
    return (; codes, seed, bonds, hopping)
end

make_probe_blocks(codes, probes, block_size, seed) =
    coded_hadamard_block_stream(codes, probes, block_size, seed)

function effective_hamiltonian(data, hartree, fock)
    rows = Vector{Int}(undef, N + 2length(data.bonds))
    columns = similar(rows)
    values = Vector{Float64}(undef, length(rows))
    position = 1
    for site in 1:N
        rows[position] = site
        columns[position] = site
        values[position] = hartree[site]
        position += 1
    end
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        coefficient = data.hopping[index] + fock[index]
        rows[position] = site
        columns[position] = neighbour
        values[position] = coefficient
        position += 1
        rows[position] = neighbour
        columns[position] = site
        values[position] = coefficient
        position += 1
    end
    return sparse(rows, columns, values, N, N)
end

function gershgorin_bounds(data, hartree, fock)
    radius = zeros(Float64, N)
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        magnitude = abs(data.hopping[index] + fock[index])
        radius[site] += magnitude
        radius[neighbour] += magnitude
    end
    return (
        minimum(hartree - radius) - SPECTRAL_MARGIN,
        maximum(hartree + radius) + SPECTRAL_MARGIN,
    )
end

function backend_hamiltonian(H, backend)
    return backend == :cuda ? CUDA.CUSPARSE.CuSparseMatrixCSR(H) : H
end

function blocked_trace_moments(
    H_backend, blocks, degree, backend, synchronize,
)
    moments = zeros(Float64, degree + 1)
    elapsed = 0.0
    total_probes = 0
    for host_probes in blocks
        total_probes += size(host_probes, 2)
        probes_backend =
            backend == :cuda ? CUDA.CuArray(host_probes) : host_probes
        synchronize()
        calculation = @timed kpm_trace_moments(
            H_backend, probes_backend, degree; synchronize=synchronize,
        )
        moments .+= size(host_probes, 2) .* calculation.value
        elapsed += calculation.time
        probes_backend = nothing
        GC.gc(false)
    end
    total_probes > 0 || error("at least one probe block is required")
    moments ./= total_probes
    return moments, elapsed
end

function blocked_local_observables(
    H_backend, blocks, coefficients, data, backend, synchronize,
)
    density_sum = zeros(Float64, N)
    bond_sum = zeros(Float64, length(data.bonds))
    kpm_time = 0.0
    transfer_time = 0.0
    measurement_time = 0.0
    total_probes = 0
    for host_probes in blocks
        block_size = size(host_probes, 2)
        total_probes += block_size
        probes_backend =
            backend == :cuda ? CUDA.CuArray(host_probes) : host_probes
        synchronize()
        calculation = @timed begin
            value = kpm_apply(
                H_backend, probes_backend, coefficients;
                synchronize=synchronize,
            )
            synchronize()
            value
        end
        filtered_backend = calculation.value
        transfer = @timed begin
            value = backend == :cuda ? Array(filtered_backend) :
                filtered_backend
            value
        end
        filtered = transfer.value
        measurement = @timed begin
            density_sum .+= vec(sum(filtered .* host_probes; dims=2))
            Threads.@threads for index in eachindex(data.bonds)
                site, neighbour, _ = data.bonds[index]
                bond_sum[index] +=
                    dot(
                        @view(filtered[site, :]),
                        @view(host_probes[neighbour, :]),
                    ) +
                    dot(
                        @view(filtered[neighbour, :]),
                        @view(host_probes[site, :]),
                    )
            end
        end
        kpm_time += calculation.time
        transfer_time += transfer.time
        measurement_time += measurement.time
        probes_backend = nothing
        filtered_backend = nothing
        filtered = nothing
        GC.gc(false)
    end
    return (
        density=density_sum ./ total_probes,
        bond_order=bond_sum ./ (2total_probes),
        kpm_time=kpm_time,
        transfer_time=transfer_time,
        measurement_time=measurement_time,
    )
end

function mean_fields(density, bond_order, data)
    hartree = zeros(Float64, N)
    for site in 1:N
        x = (site - 1) % NX
        y = div(site - 1, NX)
        x > 0 && (hartree[site] += U * density[site_index(x - 1, y)])
        x < NX - 1 && (hartree[site] += U * density[site_index(x + 1, y)])
        y > 0 && (hartree[site] += U * density[site_index(x, y - 1)])
        y < NY - 1 && (hartree[site] += U * density[site_index(x, y + 1)])
    end
    return hartree, -U .* bond_order
end

function energy(density, bond_order, data)
    kinetic = 2dot(data.hopping, bond_order)
    hartree = 0.0
    fock = 0.0
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        hartree += U * density[site] * density[neighbour]
        fock -= U * abs2(bond_order[index])
    end
    return (; kinetic, hartree, fock, total=kinetic + hartree + fock)
end

function checkerboard_order(density)
    return sum(0:(NY - 1)) do y
        sum(0:(NX - 1)) do x
            (iseven(x + y) ? 1.0 : -1.0) * density[site_index(x, y)]
        end
    end / N
end

"""Persist the production SCF result before the optional independent audit.

The audit deliberately uses more moments and probes, so it can require more
memory than the SCF itself. A failed audit must never discard a converged
production calculation.
"""
function write_production_checkpoint(
    output, density, bond_order, data, chemical_potential, converged,
    termination_reason, completed_iterations,
)
    energy_value = energy(density, bond_order, data)
    open(joinpath(output, "production_observables.toml"), "w") do io
        TOML.print(io, Dict(
            "scf_converged" => converged,
            "scf_termination_reason" => string(termination_reason),
            "scf_iterations" => completed_iterations,
            "particle_number" => sum(density),
            "chemical_potential" => chemical_potential,
            "checkerboard_order" => checkerboard_order(density),
            "energy_kinetic" => energy_value.kinetic,
            "energy_hartree" => energy_value.hartree,
            "energy_fock" => energy_value.fock,
            "energy_total" => energy_value.total,
            "final_audit_status" => "pending",
        ))
    end
    open(joinpath(output, "production_density.csv"), "w") do io
        write_csv_row(io, ("site", "x", "y", "density"))
        for y in 0:(NY - 1), x in 0:(NX - 1)
            site = site_index(x, y)
            write_csv_row(io, (site, x, y, density[site]))
        end
    end
    return energy_value
end

function main()
    synchronize = backend == :cuda ? CUDA.synchronize : () -> nothing
    device_name = backend == :cuda ? CUDA.name(CUDA.device()) : "CPU"
    device_total_memory =
        backend == :cuda ? CUDA.total_memory() : 0
    device_free_memory_before =
        backend == :cuda ? CUDA.free_memory() : 0
    data = lattice_data()
    println("Generating fixed production probe blocks")
    flush(stdout)
    production_blocks =
        make_probe_blocks(data.codes, PROBES, PROBE_BLOCK, PROBE_SEED)

    mkpath(output)
    history_path = joinpath(output, "scf_history.csv")
    open(history_path, "w") do io
        write_csv_row(io, (
            "iteration", "spectral_lower", "spectral_upper",
            "chemical_potential", "trace", "trace_error", "vh_residual",
            "vf_residual", "density_residual", "bond_residual",
            "two_cycle_residual", "mixing_method", "energy_total", "checkerboard_order",
            "trace_moment_time_s", "local_kpm_time_s",
            "transfer_to_host_time_s", "measurement_time_s",
        ))
    end

    hartree = copy(data.seed)
    fock = zeros(Float64, length(data.bonds))
    previous_density = nothing
    previous_bond_order = nothing
    density_two_steps_ago = nothing
    last_density = zeros(Float64, N)
    last_bond_order = zeros(Float64, length(data.bonds))
    stable_count = 0
    converged = false
    termination_reason = :max_iterations
    completed_iterations = 0
    last_mu = NaN
    mixer = SCF_MIXER == :pulay ? PulayMixer(
        history=PULAY_HISTORY, warmup=PULAY_WARMUP,
        regularization=PULAY_REGULARIZATION,
    ) : nothing
    input_mixing_method = :seed

    println(
        "Blocked KPM SCF: $(NX)x$(NY), N=$N, M=$MOMENTS, R=$PROBES, " *
        "block=$PROBE_BLOCK mixer=$SCF_MIXER",
    )
    flush(stdout)
    for iteration in 1:SCF_MAX_ITERATIONS
        completed_iterations = iteration
        lower, upper = gershgorin_bounds(data, hartree, fock)
        center = (lower + upper) / 2
        halfwidth = (upper - lower) / 2
        H = effective_hamiltonian(data, hartree, fock)
        scaled_H = center == 0 ? H / halfwidth :
            (H - center * sparse(I, N, N)) / halfwidth
        H_backend = backend_hamiltonian(scaled_H, backend)
        trace_moments, trace_time = blocked_trace_moments(
            H_backend, production_blocks, MOMENTS, backend, synchronize,
        )
        mu_result = find_scaled_chemical_potential(
            trace_moments, NE; tolerance=TRACE_TOLERANCE,
        )
        last_mu = center + halfwidth * mu_result.scaled_mu
        coefficients =
            projector_coefficients(MOMENTS, mu_result.scaled_mu)
        local_observables = blocked_local_observables(
            H_backend, production_blocks, coefficients, data, backend,
            synchronize,
        )
        density = local_observables.density
        bond_order = local_observables.bond_order
        new_hartree, new_fock = mean_fields(density, bond_order, data)
        vh_residual = iteration == 1 ? Inf :
            relative_change(new_hartree, hartree)
        vf_residual = iteration == 1 ? Inf :
            relative_change(new_fock, fock)
        density_residual = iteration == 1 ? Inf :
            relative_change(density, previous_density)
        bond_residual = iteration == 1 ? Inf :
            relative_change(bond_order, previous_bond_order)
        two_cycle_residual = iteration < 3 ? Inf :
            relative_change(density, density_two_steps_ago)
        energy_value = energy(density, bond_order, data)
        order = checkerboard_order(density)
        trace_value = sum(density)
        open(history_path, "a") do io
            write_csv_row(io, (
                iteration, lower, upper, last_mu, trace_value,
                abs(trace_value - NE), vh_residual, vf_residual,
                density_residual, bond_residual, two_cycle_residual,
                input_mixing_method,
                energy_value.total, order, trace_time,
                local_observables.kpm_time,
                local_observables.transfer_time,
                local_observables.measurement_time,
            ))
        end
        @printf(
            "SCF %d/%d | Tr=%.6f | μ=%.6f | VH=%.3e | VF=%.3e | n=%.3e | b=%.3e | mix=%s | E=%.6f\n",
            iteration, SCF_MAX_ITERATIONS, trace_value, last_mu,
            vh_residual, vf_residual, density_residual, bond_residual,
            string(input_mixing_method), energy_value.total,
        )
        flush(stdout)
        residuals =
            (vh_residual, vf_residual, density_residual, bond_residual)
        stable_count = all(isfinite, residuals) &&
                       maximum(residuals) <= SCF_TOLERANCE ?
            stable_count + 1 : 0
        last_density = density
        last_bond_order = bond_order
        if stable_count >= 2
            converged = true
            termination_reason = :converged
            break
        elseif iteration >= 3 &&
               two_cycle_residual <= SCF_TOLERANCE &&
               density_residual > SCF_TOLERANCE
            termination_reason = :two_cycle_detected
            break
        end
        density_two_steps_ago = previous_density
        previous_density = density
        previous_bond_order = bond_order
        if iteration == 1
            hartree = new_hartree
            fock = new_fock
            input_mixing_method = :direct
        else
            if isnothing(mixer)
                hartree = SCF_MIXING .* new_hartree .+
                          (1 - SCF_MIXING) .* hartree
                fock = SCF_MIXING .* new_fock .+
                       (1 - SCF_MIXING) .* fock
                input_mixing_method = :linear
            else
                input_fields = vcat(hartree, fock)
                output_fields = vcat(new_hartree, new_fock)
                mixed_fields, input_mixing_method = pulay_update!(
                    mixer, input_fields, output_fields;
                    damping=PULAY_DAMPING, linear_damping=SCF_MIXING,
                )
                hartree = mixed_fields[1:N]
                fock = mixed_fields[(N + 1):end]
            end
        end
        H_backend = nothing
        GC.gc()
        backend == :cuda && CUDA.reclaim()
    end

    println("Writing production checkpoint before independent final audit")
    flush(stdout)
    production_energy = write_production_checkpoint(
        output, last_density, last_bond_order, data, last_mu, converged,
        termination_reason, completed_iterations,
    )
    production_blocks = nothing
    GC.gc()
    backend == :cuda && CUDA.reclaim()
    println("Generating independent final-audit probe blocks")
    flush(stdout)
    audit_blocks =
        make_probe_blocks(data.codes, FINAL_PROBES, PROBE_BLOCK, FINAL_SEED)
    lower, upper = gershgorin_bounds(data, hartree, fock)
    center = (lower + upper) / 2
    halfwidth = (upper - lower) / 2
    H = effective_hamiltonian(data, hartree, fock)
    scaled_H = center == 0 ? H / halfwidth :
        (H - center * sparse(I, N, N)) / halfwidth
    H_backend = backend_hamiltonian(scaled_H, backend)
    audit_moments, _ = blocked_trace_moments(
        H_backend, audit_blocks, FINAL_MOMENTS, backend, synchronize,
    )
    audit_mu_result = find_scaled_chemical_potential(
        audit_moments, NE; tolerance=TRACE_TOLERANCE,
    )
    audit_mu = center + halfwidth * audit_mu_result.scaled_mu
    audit = blocked_local_observables(
        H_backend, audit_blocks,
        projector_coefficients(FINAL_MOMENTS, audit_mu_result.scaled_mu),
        data, backend, synchronize,
    )

    audit_summary = Dict(
        "moments" => FINAL_MOMENTS,
        "probes" => FINAL_PROBES,
        "seed" => FINAL_SEED,
        "chemical_potential" => audit_mu,
        "trace" => sum(audit.density),
        "trace_error" => abs(sum(audit.density) - NE),
        "density_max_abs_difference" =>
            maximum(abs, audit.density - last_density),
        "density_rms_difference" =>
            norm(audit.density - last_density) / sqrt(N),
        "bond_max_abs_difference" =>
            maximum(abs, audit.bond_order - last_bond_order),
        "bond_rms_difference" =>
            norm(audit.bond_order - last_bond_order) /
            sqrt(length(data.bonds)),
    )
    open(joinpath(output, "final_audit.toml"), "w") do io
        TOML.print(io, audit_summary)
    end
    audited_energy = energy(audit.density, audit.bond_order, data)
    open(joinpath(output, "observables.toml"), "w") do io
        TOML.print(io, Dict(
            "scf_converged" => converged,
            "scf_termination_reason" => string(termination_reason),
            "particle_number" => sum(last_density),
            "chemical_potential" => last_mu,
            "checkerboard_order" => checkerboard_order(last_density),
            "energy_kinetic" => production_energy.kinetic,
            "energy_hartree" => production_energy.hartree,
            "energy_fock" => production_energy.fock,
            "energy_total" => production_energy.total,
            "audited_energy_total" => audited_energy.total,
        ))
    end
    open(joinpath(output, "density.csv"), "w") do io
        write_csv_row(io, ("site", "x", "y", "production", "audit"))
        for y in 0:(NY - 1), x in 0:(NX - 1)
            site = site_index(x, y)
            write_csv_row(io, (
                site, x, y, last_density[site], audit.density[site],
            ))
        end
    end
    open(joinpath(output, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "solver" => "blocked_sparse_kpm_rectangular_hf",
            "backend" => string(backend),
            "device_name" => device_name,
            "device_total_memory_bytes" => device_total_memory,
            "device_free_memory_before_bytes" => device_free_memory_before,
            "nx" => NX,
            "ny" => NY,
            "matrix_dimension" => N,
            "target_particles" => NE,
            "moments" => MOMENTS,
            "probes" => PROBES,
            "probe_block" => PROBE_BLOCK,
            "scf_mixer" => string(SCF_MIXER),
            "pulay_history" => PULAY_HISTORY,
            "pulay_warmup" => PULAY_WARMUP,
            "pulay_regularization" => PULAY_REGULARIZATION,
            "pulay_damping" => PULAY_DAMPING,
            "scf_converged" => converged,
            "scf_termination_reason" => string(termination_reason),
            "finished_at" => string(now(UTC)),
        ))
    end
    println("Blocked KPM SCF complete: $output")
    println("  converged=$converged termination=$termination_reason")
    return converged
end

exit(main() ? 0 : 2)
