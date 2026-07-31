#!/usr/bin/env julia

"""Sparse-vector KPM Hartree--Fock SCF for an open square lattice.

Only the density diagonal and real nearest-neighbour density-matrix entries
are estimated. They are sufficient for the project's spinless
nearest-neighbour Hartree and real-exchange Fock fields. Fixed Hadamard probes
make the SCF map deterministic. A trace-moment pass determines the canonical
chemical potential before the local-observable pass.
"""

using Dates
using LinearAlgebra
using MPO_MeanField
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 2:3 || error(
    "usage: run_square_kpm_scf.jl CAMPAIGN_FILE OUTPUT_DIRECTORY [BACKEND]",
)
campaign_file = abspath(ARGS[1])
output = abspath(ARGS[2])
backend = length(ARGS) == 3 ? Symbol(lowercase(ARGS[3])) : :cpu
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
length(campaign.runs) == 1 || error("KPM SCF pilot expects one campaign run")
spec = only(campaign.runs)
params = spec.params
params isa ParametersSquare || error("KPM SCF requires ParametersSquare")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

const MOMENTS = parse(Int, get(ENV, "KPM_SCF_MOMENTS", "1200"))
const PROBES = parse(Int, get(ENV, "KPM_SCF_PROBES", "1024"))
const PROBE_SEED = parse(Int, get(ENV, "KPM_SCF_PROBE_SEED", "510578"))
const FINAL_MOMENTS =
    parse(Int, get(ENV, "KPM_SCF_FINAL_MOMENTS", "1600"))
const FINAL_PROBES =
    parse(Int, get(ENV, "KPM_SCF_FINAL_PROBES", "4096"))
const FINAL_SEED =
    parse(Int, get(ENV, "KPM_SCF_FINAL_SEED", "20260730"))
const SPECTRAL_MARGIN = 0.1
const SCF_MIXER = Symbol(lowercase(get(ENV, "KPM_SCF_MIXER", "linear")))
SCF_MIXER in (:linear, :pulay) ||
    error("KPM_SCF_MIXER must be linear or pulay")
const PULAY_HISTORY = parse(Int, get(ENV, "KPM_SCF_PULAY_HISTORY", "6"))
const PULAY_WARMUP = parse(Int, get(ENV, "KPM_SCF_PULAY_WARMUP", "4"))
const PULAY_REGULARIZATION = parse(
    Float64, get(ENV, "KPM_SCF_PULAY_REGULARIZATION", "1e-12"),
)
const PULAY_COEFFICIENT_LIMIT = parse(
    Float64, get(ENV, "KPM_SCF_PULAY_COEFFICIENT_LIMIT", "8.0"),
)
PULAY_COEFFICIENT_LIMIT > 0 ||
    error("KPM_SCF_PULAY_COEFFICIENT_LIMIT must be positive")
const PULAY_DIAGNOSTICS = lowercase(get(
    ENV, "KPM_SCF_PULAY_DIAGNOSTICS", "false",
)) in ("1", "true", "yes")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
relative_change(candidate, reference) =
    norm(candidate - reference) / max(norm(reference), sqrt(eps(Float64)))

function lattice_data(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    onsite = zeros(Float64, N)
    seed = zeros(Float64, N)
    bonds = Tuple{Int,Int,Symbol}[]
    hopping = Float64[]
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) :
        Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) :
        Float64(params.t[2](x, y))
    for y in 0:(side - 1), x in 0:(side - 1)
        site = square_lattice_index(x, y, params.L)
        onsite[site] = isnothing(params.W) ? 0.0 : Float64(params.W(x, y))
        seed[site] = isnothing(params.S) ? 0.0 : Float64(params.S(x, y))
        if x < side - 1
            push!(bonds, (
                site, square_lattice_index(x + 1, y, params.L), :horizontal,
            ))
            push!(hopping, tx(x, y))
        end
        if y < side - 1
            push!(bonds, (
                site, square_lattice_index(x, y + 1, params.L), :vertical,
            ))
            push!(hopping, ty(x, y))
        end
    end
    return (; side, N, onsite, seed, bonds, hopping)
end

function effective_hamiltonian(data, hartree, fock)
    rows = Vector{Int}(undef, data.N + 2length(data.bonds))
    columns = similar(rows)
    values = Vector{Float64}(undef, length(rows))
    position = 1
    for site in 1:data.N
        rows[position] = site
        columns[position] = site
        values[position] = data.onsite[site] + hartree[site]
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
    return sparse(rows, columns, values, data.N, data.N)
end

function gershgorin_bounds(data, hartree, fock)
    radius = zeros(Float64, data.N)
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        magnitude = abs(data.hopping[index] + fock[index])
        radius[site] += magnitude
        radius[neighbour] += magnitude
    end
    diagonal = data.onsite + hartree
    return (
        minimum(diagonal - radius) - SPECTRAL_MARGIN,
        maximum(diagonal + radius) + SPECTRAL_MARGIN,
    )
end

function local_estimates(filtered, probes, bonds)
    R = size(probes, 2)
    density = vec(mean(filtered .* probes; dims=2))
    bond_order = Vector{Float64}(undef, length(bonds))
    Threads.@threads for index in eachindex(bonds)
        site, neighbour, _ = bonds[index]
        bond_order[index] = (
            dot(@view(filtered[site, :]), @view(probes[neighbour, :])) +
            dot(@view(filtered[neighbour, :]), @view(probes[site, :]))
        ) / (2R)
    end
    return (; density, bond_order)
end

function mean_fields(density, bond_order, data, U)
    hartree = zeros(Float64, data.N)
    for (site, neighbour, _) in data.bonds
        hartree[site] += U * density[neighbour]
        hartree[neighbour] += U * density[site]
    end
    return hartree, -U .* bond_order
end

function hf_energy(data, density, bond_order, U)
    kinetic = dot(data.onsite, density) +
        2dot(data.hopping, bond_order)
    hartree = 0.0
    fock = 0.0
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        hartree += U * density[site] * density[neighbour]
        fock -= U * abs2(bond_order[index])
    end
    return (;
        kinetic, hartree, fock, interaction=hartree + fock,
        total=kinetic + hartree + fock,
    )
end

function checkerboard_order(density, data, L)
    return sum(0:(data.side - 1)) do y
        sum(0:(data.side - 1)) do x
            sign = iseven(x + y) ? 1.0 : -1.0
            sign * density[square_lattice_index(x, y, L)]
        end
    end / data.N
end

function main()
device_name = "CPU"
device_total_memory = 0
device_free_memory_before = 0
synchronize = () -> nothing
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    device_total_memory = CUDA.total_memory()
    device_free_memory_before = CUDA.free_memory()
    synchronize = CUDA.synchronize
end

data = lattice_data(params)
PROBES <= data.N || error("KPM_SCF_PROBES=$PROBES exceeds N=$(data.N)")
FINAL_PROBES <= data.N ||
    error("KPM_SCF_FINAL_PROBES=$FINAL_PROBES exceeds N=$(data.N)")
MOMENTS >= 2 || error("KPM_SCF_MOMENTS must be at least 2")
FINAL_MOMENTS >= 2 || error("KPM_SCF_FINAL_MOMENTS must be at least 2")
Ne = round(Int, params.density * data.N)
trace_tolerance = max(1e-6 * Ne, 1e-6)
scf_tolerance = params.scf_tol / 100
max_scf_iterations = min(
    params.scf_max_iterations,
    parse(Int, get(
        ENV, "KPM_SCF_MAX_ITERATIONS", string(params.scf_max_iterations),
    )),
)
max_scf_iterations > 0 || error("KPM_SCF_MAX_ITERATIONS must be positive")
pulay_damping = parse(
    Float64, get(ENV, "KPM_SCF_PULAY_DAMPING", string(params.scf_mixing)),
)
0 < pulay_damping <= 1 ||
    error("KPM_SCF_PULAY_DAMPING must lie in (0, 1]")
host_probes = probing_matrix(data.N, PROBES, :hadamard, PROBE_SEED)
probes_backend = backend == :cuda ? CUDA.CuArray(host_probes) : host_probes
synchronize()

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
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
last_density = zeros(Float64, data.N)
last_bond_order = zeros(Float64, length(data.bonds))
converged = false
termination_reason = :max_iterations
stable_count = 0
last_mu = NaN
mixer = SCF_MIXER == :pulay ? PulayMixer(
    history=PULAY_HISTORY, warmup=PULAY_WARMUP,
    regularization=PULAY_REGULARIZATION,
    coefficient_limit=PULAY_COEFFICIENT_LIMIT,
) : nothing
pulay_diagnostics_path = joinpath(output, "pulay_diagnostics.csv")
if PULAY_DIAGNOSTICS && !isnothing(mixer)
    open(pulay_diagnostics_path, "w") do io
        write_csv_row(io, (
            "iteration", "applied_next_iteration", "status", "history_size",
            "gram_condition", "max_abs_coefficient",
            "candidate_to_linear_step_ratio",
        ))
    end
end
input_mixing_method = :seed

println(
    "KPM SCF: side=$(data.side) N=$(data.N) Ne=$Ne " *
    "M=$MOMENTS R=$PROBES mixing=$(params.scf_mixing) mixer=$SCF_MIXER",
)
flush(stdout)

for iteration in 1:max_scf_iterations
    lower, upper = gershgorin_bounds(data, hartree, fock)
    center = (lower + upper) / 2
    halfwidth = (upper - lower) / 2
    H = effective_hamiltonian(data, hartree, fock)
    scaled_H = center == 0 ? H / halfwidth :
        (H - center * sparse(I, data.N, data.N)) / halfwidth
    H_backend = backend == :cuda ?
        CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H) : scaled_H

    synchronize()
    trace_calculation = @timed kpm_trace_moments(
        H_backend, probes_backend, MOMENTS; synchronize=synchronize,
    )
    trace_moments = trace_calculation.value
    chemical_potential = find_scaled_chemical_potential(
        trace_moments, Ne; tolerance=trace_tolerance,
    )
    scaled_mu = chemical_potential.scaled_mu
    last_mu = center + halfwidth * scaled_mu
    coefficients = projector_coefficients(MOMENTS, scaled_mu)
    GC.gc()
    backend == :cuda && CUDA.reclaim()

    synchronize()
    kpm_calculation = @timed begin
        value = kpm_apply(
            H_backend, probes_backend, coefficients;
            synchronize=synchronize,
        )
        synchronize()
        value
    end
    filtered_backend = kpm_calculation.value
    transfer = @timed begin
        value = backend == :cuda ? Array(filtered_backend) : filtered_backend
        value
    end
    filtered = transfer.value
    measurement = @timed local_estimates(filtered, host_probes, data.bonds)
    density = measurement.value.density
    bond_order = measurement.value.bond_order
    new_hartree, new_fock =
        mean_fields(density, bond_order, data, Float64(params.U))

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
    energy = hf_energy(data, density, bond_order, Float64(params.U))
    order = checkerboard_order(density, data, params.L)
    trace_value = sum(density)

    open(history_path, "a") do io
        write_csv_row(io, (
            iteration, lower, upper, last_mu, trace_value,
            abs(trace_value - Ne), vh_residual, vf_residual,
            density_residual, bond_residual, two_cycle_residual,
            input_mixing_method, energy.total,
            order, trace_calculation.time, kpm_calculation.time,
            transfer.time, measurement.time,
        ))
    end
    @printf(
        "SCF %d/%d | Tr=%.9f | μ=%.6f | VH=%.3e | VF=%.3e | n=%.3e | b=%.3e | mix=%s | E=%.9f | bounds=[%.3f,%.3f]\n",
        iteration, max_scf_iterations, trace_value, last_mu,
        vh_residual, vf_residual, density_residual, bond_residual,
        string(input_mixing_method), energy.total, lower, upper,
    )
    flush(stdout)

    residuals = (vh_residual, vf_residual, density_residual, bond_residual)
    if all(isfinite, residuals) && maximum(residuals) <= scf_tolerance
        stable_count += 1
    else
        stable_count = 0
    end
    last_density = density
    last_bond_order = bond_order
    if stable_count >= 2
        converged = true
        termination_reason = :converged
        break
    elseif iteration >= 3 && two_cycle_residual <= scf_tolerance &&
           density_residual > scf_tolerance
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
        mixing = params.scf_mixing
        if isnothing(mixer)
            hartree = mixing .* new_hartree .+ (1 - mixing) .* hartree
            fock = mixing .* new_fock .+ (1 - mixing) .* fock
            input_mixing_method = :linear
        else
            input_fields = vcat(hartree, fock)
            output_fields = vcat(new_hartree, new_fock)
            mixed_fields, input_mixing_method = pulay_update!(
                mixer, input_fields, output_fields;
                damping=pulay_damping, linear_damping=mixing,
                diagnostics=PULAY_DIAGNOSTICS,
            )
            if PULAY_DIAGNOSTICS
                diagnostic = pulay_last_diagnostics(mixer)
                open(pulay_diagnostics_path, "a") do io
                    write_csv_row(io, (
                        iteration, iteration + 1, diagnostic.status,
                        diagnostic.history_size, diagnostic.gram_condition,
                        diagnostic.max_abs_coefficient,
                        diagnostic.candidate_to_linear_step_ratio,
                    ))
                end
            end
            hartree = mixed_fields[1:data.N]
            fock = mixed_fields[(data.N + 1):end]
        end
    end

    filtered_backend = nothing
    filtered = nothing
    H_backend = nothing
    GC.gc()
    backend == :cuda && CUDA.reclaim()
end

println("Final independent fixed-field audit: M=$FINAL_MOMENTS R=$FINAL_PROBES")
flush(stdout)
final_host_probes =
    probing_matrix(data.N, FINAL_PROBES, :hadamard, FINAL_SEED)
final_probes_backend = backend == :cuda ?
    CUDA.CuArray(final_host_probes) : final_host_probes
lower, upper = gershgorin_bounds(data, hartree, fock)
center = (lower + upper) / 2
halfwidth = (upper - lower) / 2
H = effective_hamiltonian(data, hartree, fock)
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, data.N, data.N)) / halfwidth
H_backend = backend == :cuda ?
    CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H) : scaled_H
final_trace_moments = kpm_trace_moments(
    H_backend, final_probes_backend, FINAL_MOMENTS; synchronize=synchronize,
)
final_mu_result = find_scaled_chemical_potential(
    final_trace_moments, Ne; tolerance=trace_tolerance,
)
final_mu = center + halfwidth * final_mu_result.scaled_mu
final_coefficients =
    projector_coefficients(FINAL_MOMENTS, final_mu_result.scaled_mu)
GC.gc()
backend == :cuda && CUDA.reclaim()
final_filtered_backend = kpm_apply(
    H_backend, final_probes_backend, final_coefficients;
    synchronize=synchronize,
)
final_filtered = backend == :cuda ?
    Array(final_filtered_backend) : final_filtered_backend
final_measurement =
    local_estimates(final_filtered, final_host_probes, data.bonds)
final_density = final_measurement.density
final_bond_order = final_measurement.bond_order

audit = Dict(
    "moments" => FINAL_MOMENTS,
    "probes" => FINAL_PROBES,
    "seed" => FINAL_SEED,
    "chemical_potential" => final_mu,
    "trace" => sum(final_density),
    "trace_error" => abs(sum(final_density) - Ne),
    "density_max_abs_difference" =>
        maximum(abs, final_density - last_density),
    "density_rms_difference" =>
        norm(final_density - last_density) / sqrt(data.N),
    "bond_max_abs_difference" =>
        maximum(abs, final_bond_order - last_bond_order),
    "bond_rms_difference" =>
        norm(final_bond_order - last_bond_order) /
        sqrt(length(data.bonds)),
)
open(joinpath(output, "final_audit.toml"), "w") do io
    TOML.print(io, audit)
end

final_energy =
    hf_energy(data, final_density, final_bond_order, Float64(params.U))
open(joinpath(output, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "scf_converged" => converged,
        "scf_termination_reason" => string(termination_reason),
        "particle_number" => sum(last_density),
        "chemical_potential" => last_mu,
        "checkerboard_order" =>
            checkerboard_order(last_density, data, params.L),
        "energy_kinetic" => hf_energy(
            data, last_density, last_bond_order, Float64(params.U),
        ).kinetic,
        "energy_hartree" => hf_energy(
            data, last_density, last_bond_order, Float64(params.U),
        ).hartree,
        "energy_fock" => hf_energy(
            data, last_density, last_bond_order, Float64(params.U),
        ).fock,
        "energy_total" => hf_energy(
            data, last_density, last_bond_order, Float64(params.U),
        ).total,
        "audited_energy_total" => final_energy.total,
    ))
end
open(joinpath(output, "density.csv"), "w") do io
    write_csv_row(io, ("site", "x", "y", "production", "audit"))
    for y in 0:(data.side - 1), x in 0:(data.side - 1)
        site = square_lattice_index(x, y, params.L)
        write_csv_row(io, (
            site, x, y, last_density[site], final_density[site],
        ))
    end
end
open(joinpath(output, "bond_order.csv"), "w") do io
    write_csv_row(io, (
        "site", "neighbour", "orientation", "production", "audit",
    ))
    for index in eachindex(data.bonds)
        site, neighbour, orientation = data.bonds[index]
        write_csv_row(io, (
            site, neighbour, orientation, last_bond_order[index],
            final_bond_order[index],
        ))
    end
end

open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => campaign.name,
        "label" => spec.label,
        "solver" => "sparse_kpm_square_hf",
        "backend" => string(backend),
        "device_name" => device_name,
        "device_total_memory_bytes" => device_total_memory,
        "device_free_memory_before_bytes" => device_free_memory_before,
        "matrix_dimension" => data.N,
        "target_particles" => Ne,
        "moments" => MOMENTS,
        "probes" => PROBES,
        "probe_seed" => PROBE_SEED,
        "trace_relative_tolerance" => 1e-6,
        "scf_relative_tolerance" => scf_tolerance,
        "scf_max_iterations" => max_scf_iterations,
        "scf_mixer" => string(SCF_MIXER),
        "pulay_history" => PULAY_HISTORY,
        "pulay_warmup" => PULAY_WARMUP,
            "pulay_regularization" => PULAY_REGULARIZATION,
            "pulay_coefficient_limit" => PULAY_COEFFICIENT_LIMIT,
            "pulay_damping" => pulay_damping,
        "spectral_policy" => "per_iteration_gershgorin_with_margin",
        "spectral_margin" => SPECTRAL_MARGIN,
        "scf_converged" => converged,
        "scf_termination_reason" => string(termination_reason),
        "finished_at" => string(now(UTC)),
    ))
end

println("KPM SCF complete: $output")
println("  converged=$converged termination=$termination_reason")
return converged
end

exit(main() ? 0 : 2)
