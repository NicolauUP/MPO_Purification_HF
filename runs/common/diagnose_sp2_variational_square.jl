#!/usr/bin/env julia

"""Compare zip-up and globally fitted MPO squares along fixed-H SP2 trajectories.

This is an isolated small-system diagnostic, not a production purification
method. At every iteration it first forms a higher-cap reference square. The
normal `apply(rho, rho)` square and a variational compression of that reference
are then compared at the same target bond dimension. Separate SP2 trajectories
are advanced with the two compressed squares.
"""

using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using TOML

length(ARGS) in 8:10 || error(
    "usage: diagnose_sp2_variational_square.jl CAMPAIGN_FILE TASK_INDEX " *
    "TARGET_MAXDIM REFERENCE_MAXDIM CUTOFF NSWEEPS STEPS OUTPUT_DIRECTORY " *
    "[SPECTRAL_BOUND [BACKEND]]",
)

campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
target_maxdim = parse(Int, ARGS[3])
reference_maxdim = parse(Int, ARGS[4])
cutoff = parse(Float64, ARGS[5])
nsweeps = parse(Int, ARGS[6])
steps = parse(Int, ARGS[7])
output = abspath(ARGS[8])
bound_override = length(ARGS) == 9 ? parse(Float64, ARGS[9]) : nothing
if length(ARGS) == 10
    bound_override = parse(Float64, ARGS[9])
end
backend = length(ARGS) == 10 ? Symbol(lowercase(ARGS[10])) : :cpu

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
target_maxdim > 0 || error("TARGET_MAXDIM must be positive")
reference_maxdim >= target_maxdim ||
    error("REFERENCE_MAXDIM must be at least TARGET_MAXDIM")
isfinite(cutoff) && cutoff >= 0 || error("CUTOFF must be finite and nonnegative")
nsweeps > 0 || error("NSWEEPS must be positive")
steps > 0 || error("STEPS must be positive")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("diagnostic requires ParametersSquare")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

to_device = backend == :cuda ?
    (value -> ITensors.adapt(CUDA.CuArray, value)) : identity
synchronize_backend() = backend == :cuda ? CUDA.synchronize() : nothing
device_name = backend == :cuda ? CUDA.name(CUDA.device()) : "CPU"
device_total_memory = backend == :cuda ? CUDA.total_memory() : 0
device_free_memory_before = backend == :cuda ? CUDA.free_memory() : 0

function backend_timed(f)
    synchronize_backend()
    measurement = @timed begin
        value = f()
        synchronize_backend()
        value
    end
    return measurement
end

function with_numerics(params::ParametersSquare, maxdim::Int, cutoff::Float64)
    return ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=cutoff, itensors_maxdim=maxdim,
        density=params.density, purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end

mean_chi(mpo::MPO) = length(mpo) <= 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) /
    (length(mpo) - 1)

function direct_initial_hamiltonian(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    H = zeros(Float64, N, N)
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) :
        Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) :
        Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, params.L)
        H[site, site] = (isnothing(params.W) ? 0.0 : Float64(params.W(x, y))) +
                        (isnothing(params.S) ? 0.0 : Float64(params.S(x, y)))
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, params.L)
            H[site, neighbour] = H[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, params.L)
            H[site, neighbour] = H[neighbour, site] = ty(x, y)
        end
    end
    return H
end

function matrix_to_interleaved_mpo(matrix, sites; cutoff::Float64, maxdim::Int)
    L = length(sites)
    N = 2^L
    values = zeros(eltype(matrix), ntuple(_ -> 2, 2L))
    coordinates = Vector{Int}(undef, 2L)
    for row in 1:N, column in 1:N
        for site in 1:L
            shift = L - site
            coordinates[2site - 1] = ((row - 1) >> shift & 1) + 1
            coordinates[2site] = ((column - 1) >> shift & 1) + 1
        end
        values[coordinates...] = matrix[row, column]
    end
    physical = Index[]
    for site in sites
        push!(physical, prime(site), dag(site))
    end
    return MPO(
        ITensor(values, physical...), sites;
        cutoff=cutoff, maxdim=maxdim,
    )
end

function vectorize_mpo(mpo::MPO, vector_sites)
    tensors = ITensor[]
    combiners = ITensor[]
    for site in eachindex(mpo)
        physical = siteinds(mpo, site)
        comb = combiner(physical...)
        tensor = mpo[site] * comb
        combined = only(uniqueinds(comb, physical))
        push!(tensors, replaceind(tensor, combined => vector_sites[site]))
        push!(combiners, comb)
    end
    return MPS(tensors), combiners
end

function devectorize_mps(state::MPS, combiners)
    tensors = ITensor[]
    for site in eachindex(state)
        combined = first(inds(combiners[site]))
        physical_site = only(siteinds(state, site))
        tensor = replaceind(state[site], physical_site => combined) * dag(combiners[site])
        push!(tensors, tensor)
    end
    return MPO(tensors)
end

function left_multiplication_superoperator(
    left::MPO,
    input_sites,
    output_sites,
)
    tensors = ITensor[]
    output_combiners = ITensor[]
    spectator_outputs = Index[]
    for site in eachindex(left)
        physical = siteinds(left, site)
        output = only(filter(index -> plev(index) > 0, physical))
        input = only(filter(index -> plev(index) == 0, physical))
        spectator_input = sim(input; tags="VectorInputSpectator,n=$site")
        spectator_output = sim(input; tags="VectorOutputSpectator,n=$site")
        # CUDA cannot contract a device DenseTensor with ITensor's
        # scalar-backed DiagTensor storage. Densify these small structural
        # tensors on the host before transferring them to the selected backend.
        input_combiner =
            to_device(dense(combiner(input, spectator_input)))
        output_combiner =
            to_device(dense(combiner(output, spectator_output)))
        spectator_identity =
            to_device(dense(delta(spectator_input, spectator_output)))
        tensor = left[site] * spectator_identity *
                 input_combiner * output_combiner
        combined_input = only(uniqueinds(
            input_combiner, (input, spectator_input),
        ))
        combined_output = only(uniqueinds(
            output_combiner, (output, spectator_output),
        ))
        tensor = replaceinds(
            tensor,
            (combined_input, combined_output) =>
                (input_sites[site], output_sites[site]),
        )
        push!(tensors, tensor)
        push!(output_combiners, output_combiner)
        push!(spectator_outputs, spectator_output)
    end
    # `left` and the vectorized right factor originate from the same MPO and
    # therefore initially share link-index identities. They are independent
    # tensor-network layers and must not contract those links with each other.
    return sim(linkinds, MPO(tensors)), output_combiners, spectator_outputs
end

function devectorize_left_product(
    state::MPS,
    template::MPO,
    output_combiners,
    spectator_outputs,
)
    tensors = ITensor[]
    for site in eachindex(state)
        physical = siteinds(template, site)
        output = only(filter(index -> plev(index) > 0, physical))
        input = only(filter(index -> plev(index) == 0, physical))
        combined_output = only(uniqueinds(
            output_combiners[site], (output, spectator_outputs[site]),
        ))
        state_site = only(siteinds(state, site))
        tensor = replaceind(state[site], state_site => combined_output) *
                 dag(output_combiners[site])
        tensor = replaceind(tensor, spectator_outputs[site] => input)
        push!(tensors, tensor)
    end
    return MPO(tensors)
end

function variational_compress(
    target::MPO,
    initial::MPO,
    vector_sites,
    identity_operator;
    cutoff::Float64,
    maxdim::Int,
    nsweeps::Int,
)
    target_vector, combiners = vectorize_mpo(target, vector_sites)
    initial_vector, _ = vectorize_mpo(initial, vector_sites)
    fitted = apply(
        identity_operator, target_vector;
        alg="fit", init=initial_vector, nsweeps=nsweeps,
        cutoff=cutoff, maxdim=maxdim,
    )
    return devectorize_mps(fitted, combiners)
end

function implicit_variational_square(
    rho::MPO,
    initial::MPO;
    cutoff::Float64,
    maxdim::Int,
    nsweeps::Int,
)
    input_sites = [
        Index(4, "VectorInput,Site,n=$site") for site in eachindex(rho)
    ]
    output_sites = [
        Index(4, "VectorOutput,Site,n=$site") for site in eachindex(rho)
    ]
    input_state, _ = vectorize_mpo(rho, input_sites)
    initial_state, _ = vectorize_mpo(initial, output_sites)
    superoperator, output_combiners, spectator_outputs =
        left_multiplication_superoperator(rho, input_sites, output_sites)
    fitted = contract(
        superoperator, input_state;
        alg="fit", init=initial_state, nsweeps=nsweeps,
        cutoff=cutoff, maxdim=maxdim,
    )
    return devectorize_left_product(
        fitted, rho, output_combiners, spectator_outputs,
    )
end

function relative_difference(left::MPO, right::MPO)
    denominator = max(real(inner(right, right)), eps(Float64))
    numerator = max(real(inner(left, left)) + real(inner(right, right)) -
                    2real(inner(left, right)), 0.0)
    return sqrt(numerator / denominator)
end

function state_metrics(rho::MPO, exact_projector::MPO, params, Ne)
    rho_squared = apply(rho, rho; cutoff=cutoff, maxdim=target_maxdim)
    return (
        trace=real(tr(rho)),
        trace_error=abs(real(tr(rho)) - Ne),
        idempotency=MPO_MeanField.idempotency_residual(rho, rho_squared),
        exact_error=relative_difference(rho, exact_projector),
        max_chi=maxlinkdim(rho),
        mean_chi=mean_chi(rho),
    )
end

params = with_numerics(spec.params, target_maxdim, cutoff)
N = 2^params.L
N <= 1024 || error("dense exact-projector diagnostic is restricted to N <= 1024")
Ne = round(Int, N * params.density)
campaign_bounds = validate_spectral_bounds(spec.spectral_bounds...)
bounds = isnothing(bound_override) ? campaign_bounds :
    validate_spectral_bounds(-bound_override, bound_override)

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
started_at = now(UTC)

H = direct_initial_hamiltonian(params)
eigenpairs = eigen(Symmetric(H))
occupied = @view eigenpairs.vectors[:, 1:Ne]
exact_matrix = Matrix(occupied * occupied')

sys = System(params)
exact_projector = matrix_to_interleaved_mpo(
    exact_matrix, sys.sites; cutoff=min(cutoff, 1e-14), maxdim=reference_maxdim,
)
rho0_cpu = construct_rho_0(
    sys, params, bounds...; method=:sp2, verify_spectral_bounds=false,
)
exact_projector = to_device(exact_projector)
rho0 = to_device(rho0_cpu)
device_scalar_type = ITensors.scalartype(rho0)
backend == :cuda && device_scalar_type != Float64 && error(
    "CUDA diagnostic requires Float64 MPO tensors, got $device_scalar_type",
)
rho_zip = deepcopy(rho0)
rho_fit = deepcopy(rho0)
rho_implicit = deepcopy(rho0)
vector_sites = siteinds("Qudit", params.L; dim=4)
identity_operator = to_device(MPO(vector_sites, "Id"))
trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)

open(joinpath(output, "iterations.csv"), "w") do io
    write_csv_row(io, (
        "method", "iteration", "branch", "trace", "trace_error",
        "idempotency_residual", "relative_exact_projector_error",
        "rho_max_chi", "rho_mean_chi", "reference_square_max_chi",
        "compressed_square_max_chi", "square_relative_error",
        "reference_square_time_s", "compression_time_s",
        "compression_allocations_bytes",
    ))
    for iteration in 1:steps
        for method in (:zipup, :variational, :implicit)
            global rho_zip, rho_fit, rho_implicit
            rho = method == :zipup ? rho_zip :
                (method == :variational ? rho_fit : rho_implicit)
            reference_timing = backend_timed() do
                apply(
                    rho, rho; cutoff=min(cutoff, 1e-14),
                    maxdim=reference_maxdim,
                )
            end
            reference_square = reference_timing.value
            compression = backend_timed() do
                if method == :zipup
                    apply(rho, rho; cutoff=cutoff, maxdim=target_maxdim)
                elseif method == :variational
                    initial = apply(rho, rho; cutoff=cutoff, maxdim=target_maxdim)
                    variational_compress(
                        reference_square, initial, vector_sites, identity_operator;
                        cutoff=cutoff, maxdim=target_maxdim, nsweeps=nsweeps,
                    )
                else
                    initial = apply(rho, rho; cutoff=cutoff, maxdim=target_maxdim)
                    implicit_variational_square(
                        rho, initial;
                        cutoff=cutoff, maxdim=target_maxdim, nsweeps=nsweeps,
                    )
                end
            end
            compressed_square = compression.value
            metrics = state_metrics(rho, exact_projector, params, Ne)
            trace_squared = real(tr(compressed_square))
            branch = MPO_MeanField._sp2_branch(
                metrics.trace, trace_squared, Ne, trace_tolerance,
            )
            write_csv_row(io, (
                method, iteration, branch, metrics.trace, metrics.trace_error,
                metrics.idempotency, metrics.exact_error,
                metrics.max_chi, metrics.mean_chi, maxlinkdim(reference_square),
                maxlinkdim(compressed_square),
                relative_difference(compressed_square, reference_square),
                reference_timing.time, compression.time, compression.bytes,
            ))
            flush(io)
            next_rho = branch == :square ? compressed_square : +(
                2.0 * rho, -compressed_square;
                cutoff=cutoff, maxdim=target_maxdim,
            )
            if method == :zipup
                rho_zip = next_rho
            elseif method == :variational
                rho_fit = next_rho
            else
                rho_implicit = next_rho
            end
        end
    end
end

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "diagnostic" => "sp2_zipup_vs_variational_square",
        "backend" => string(backend), "device_name" => device_name,
        "device_scalar_type" => string(device_scalar_type),
        "device_total_memory_bytes" => device_total_memory,
        "device_free_memory_before_bytes" => device_free_memory_before,
        "device_free_memory_after_bytes" =>
            (backend == :cuda ? CUDA.free_memory() : 0),
        "campaign" => string(campaign.name), "label" => string(spec.label),
        "task_index" => task_index, "matrix_dimension" => N,
        "target_particles" => Ne, "target_maxdim" => target_maxdim,
        "reference_maxdim" => reference_maxdim, "cutoff" => cutoff,
        "variational_sweeps" => nsweeps, "steps" => steps,
        "spectral_lower" => bounds[1], "spectral_upper" => bounds[2],
        "exact_projector_max_chi" => maxlinkdim(exact_projector),
        "exact_fermi_gap" => eigenpairs.values[Ne + 1] - eigenpairs.values[Ne],
        "started_at" => string(started_at), "finished_at" => string(now(UTC)),
    ))
end

println("SP2 square-compression diagnostic written to $output")
