#!/usr/bin/env julia

"""Run one fixed-H SP2 trajectory with an implicit variational MPO square.

The variational fit targets `rho * rho` through a left-multiplication
superoperator. It never constructs a higher-cap reference square or an exact
projector. The fitted square is reused for the SP2 branch decision and the
idempotency residual.
"""

using Dates
using ITensors
using ITensorMPS
using MPO_MeanField
using TOML

length(ARGS) in 7:11 || error(
    "usage: run_sp2_implicit_square.jl CAMPAIGN_FILE TASK_INDEX MAXDIM " *
    "CUTOFF NSWEEPS STEPS OUTPUT_DIRECTORY " *
    "[SPECTRAL_BOUND [BACKEND [IDEMPOTENCY_TOLERANCE " *
    "[RELATIVE_TRACE_TOLERANCE]]]]",
)

campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
maxdim = parse(Int, ARGS[3])
cutoff = parse(Float64, ARGS[4])
nsweeps = parse(Int, ARGS[5])
steps = parse(Int, ARGS[6])
output = abspath(ARGS[7])
bound_override = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : nothing
backend = length(ARGS) >= 9 ? Symbol(lowercase(ARGS[9])) : :cpu
idempotency_tolerance = length(ARGS) >= 10 ? parse(Float64, ARGS[10]) : 2e-4
relative_trace_tolerance = length(ARGS) >= 11 ? parse(Float64, ARGS[11]) : 1e-6

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
maxdim > 0 || error("MAXDIM must be positive")
isfinite(cutoff) && cutoff >= 0 || error("CUTOFF must be finite and nonnegative")
nsweeps > 0 || error("NSWEEPS must be positive")
steps > 0 || error("STEPS must be positive")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
isfinite(idempotency_tolerance) && idempotency_tolerance > 0 ||
    error("IDEMPOTENCY_TOLERANCE must be finite and positive")
isfinite(relative_trace_tolerance) && relative_trace_tolerance > 0 ||
    error("RELATIVE_TRACE_TOLERANCE must be finite and positive")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("diagnostic requires ParametersSquare")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

to_device = backend == :cuda ?
    (value -> ITensors.adapt(CUDA.CuArray, value)) : identity
synchronize_backend() = backend == :cuda ? CUDA.synchronize() : nothing
device_free_memory() = backend == :cuda ? CUDA.free_memory() : 0

function backend_timed(f)
    synchronize_backend()
    measurement = @timed begin
        value = f()
        synchronize_backend()
        value
    end
    return measurement
end

function with_numerics(params::ParametersSquare)
    return ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=cutoff, itensors_maxdim=maxdim,
        density=params.density, purification_steps=steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end

mean_chi(mpo::MPO) = length(mpo) <= 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) /
    (length(mpo) - 1)

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

function vectorize_mpo(mpo::MPO, vector_sites)
    tensors = ITensor[]
    for site in eachindex(mpo)
        physical = siteinds(mpo, site)
        comb = combiner(physical...)
        tensor = mpo[site] * comb
        combined = only(uniqueinds(comb, physical))
        push!(tensors, replaceind(tensor, combined => vector_sites[site]))
    end
    return MPS(tensors)
end

function left_multiplication_superoperator(left::MPO, input_sites, output_sites)
    tensors = ITensor[]
    output_combiners = ITensor[]
    spectator_outputs = Index[]
    for site in eachindex(left)
        physical = siteinds(left, site)
        output = only(filter(index -> plev(index) > 0, physical))
        input = only(filter(index -> plev(index) == 0, physical))
        spectator_input = sim(input; tags="VectorInputSpectator,n=$site")
        spectator_output = sim(input; tags="VectorOutputSpectator,n=$site")
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
    # The two factors originate from the same MPO but represent independent
    # tensor-network layers, so their link-index identities must be distinct.
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

function implicit_variational_square(
    rho::MPO;
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
    input_state = vectorize_mpo(rho, input_sites)
    # Near a projector, rho itself is a physically motivated initial guess for
    # rho² and avoids constructing even a target-cap zip-up product.
    initial_state = vectorize_mpo(deepcopy(rho), output_sites)
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

params = with_numerics(spec.params)
N = 2^params.L
Ne = round(Int, N * params.density)
0 < Ne < N || error("SP2 requires 0 < Ne < N")
campaign_bounds = validate_spectral_bounds(spec.spectral_bounds...)
bounds = isnothing(bound_override) ? campaign_bounds :
    validate_spectral_bounds(-bound_override, bound_override)
branch_trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)
convergence_trace_tolerance = relative_trace_tolerance * max(1, Ne)
hermiticity_tolerance = MPO_MeanField._sp2_hermiticity_tolerance(params)

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
started_at = now(UTC)
sys = System(params)
initialization = backend_timed() do
    to_device(construct_rho_0(
        sys, params, bounds...;
        method=:sp2, verify_spectral_bounds=false,
    ))
end
rho = initialization.value
device_scalar_type = ITensors.scalartype(rho)
backend == :cuda && device_scalar_type != Float64 && error(
    "CUDA diagnostic requires Float64 MPO tensors, got $device_scalar_type",
)

converged = false
termination_reason = :max_iterations
completed_iterations = 0

open(joinpath(output, "iterations.csv"), "w") do io
    write_csv_row(io, (
        "iteration", "branch", "trace", "trace_squared", "trace_hole",
        "trace_error", "relative_trace_error", "idempotency_residual",
        "hermiticity_residual", "rho_max_chi", "rho_mean_chi",
        "square_max_chi", "square_mean_chi", "fit_time_s",
        "fit_allocations_bytes", "fit_gc_time_s", "update_time_s",
        "update_allocations_bytes", "update_gc_time_s",
        "device_free_memory_bytes",
    ))
    for iteration in 1:steps
        global rho, converged, termination_reason, completed_iterations
        fit = backend_timed() do
            implicit_variational_square(
                rho; cutoff=cutoff, maxdim=maxdim, nsweeps=nsweeps,
            )
        end
        rho_squared = fit.value
        trace_value = real(tr(rho))
        trace_squared = real(tr(rho_squared))
        trace_hole = 2trace_value - trace_squared
        trace_error = abs(trace_value - Ne)
        relative_trace_error = trace_error / max(1, Ne)
        idempotency =
            MPO_MeanField.idempotency_residual(rho, rho_squared)
        hermiticity =
            MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
        branch = MPO_MeanField._sp2_branch(
            trace_value, trace_squared, Ne, branch_trace_tolerance,
        )

        if trace_error <= convergence_trace_tolerance &&
           idempotency <= idempotency_tolerance &&
           hermiticity <= hermiticity_tolerance
            update = (value=rho, time=0.0, bytes=0, gctime=0.0)
            converged = true
            termination_reason = :idempotency_threshold
        else
            update = backend_timed() do
                branch == :square ? rho_squared : +(
                    2.0 * rho, -rho_squared;
                    cutoff=cutoff, maxdim=maxdim,
                )
            end
        end

        write_csv_row(io, (
            iteration, branch, trace_value, trace_squared, trace_hole,
            trace_error, relative_trace_error, idempotency, hermiticity,
            maxlinkdim(rho), mean_chi(rho), maxlinkdim(rho_squared),
            mean_chi(rho_squared), fit.time, fit.bytes, fit.gctime,
            update.time, update.bytes, update.gctime, device_free_memory(),
        ))
        flush(io)
        completed_iterations = iteration
        if converged || iteration == steps
            break
        end
        rho = update.value
    end
end

final_trace = real(tr(rho))
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "diagnostic" => "fixed_h_sp2_implicit_square",
        "backend" => string(backend),
        "device_name" => backend == :cuda ? CUDA.name(CUDA.device()) : "CPU",
        "device_scalar_type" => string(device_scalar_type),
        "device_total_memory_bytes" =>
            (backend == :cuda ? CUDA.total_memory() : 0),
        "device_free_memory_after_bytes" => device_free_memory(),
        "campaign" => string(campaign.name),
        "label" => string(spec.label),
        "task_index" => task_index,
        "matrix_dimension" => N,
        "target_particles" => Ne,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "variational_sweeps" => nsweeps,
        "maximum_steps" => steps,
        "completed_iterations" => completed_iterations,
        "converged" => converged,
        "termination_reason" => string(termination_reason),
        "idempotency_tolerance" => idempotency_tolerance,
        "relative_trace_tolerance" => relative_trace_tolerance,
        "spectral_lower" => bounds[1],
        "spectral_upper" => bounds[2],
        "initialization_time_s" => initialization.time,
        "initialization_allocations_bytes" => initialization.bytes,
        "final_trace" => final_trace,
        "final_trace_error" => abs(final_trace - Ne),
        "final_rho_max_chi" => maxlinkdim(rho),
        "final_rho_mean_chi" => mean_chi(rho),
        "started_at" => string(started_at),
        "finished_at" => string(now(UTC)),
    ))
end

println("Implicit-only SP2 diagnostic written to $output")
