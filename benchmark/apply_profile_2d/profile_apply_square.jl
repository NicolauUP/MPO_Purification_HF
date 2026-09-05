#!/usr/bin/env julia

"""Profile MPO squaring with the installed ITensorMPS zip-up algorithm.

This diagnostic does not modify `ITensorMPS.apply` or the purification code.
It first times the production operation

    apply(rho, rho; cutoff, maxdim)

on either a controlled synthetic MPO or a representative SP2 trajectory state.
It then mirrors the installed
`ITensors.contract(::Algorithm\"zipup\", ::MPO, ::AbstractMPS)` implementation,
placing synchronized timing boundaries around orthogonalization, every local
contraction, every `factorize`, and the final `truncate!` sweep. The resulting
MPO is checked against the production `apply` result before timings are
reported.

The phase timings deliberately synchronize after each phase. They are useful
for attribution, but the unsplit production `apply` timing is the performance
number to compare between versions or settings.

Usage:
  julia --project=. benchmark/apply_profile_2d/profile_apply_square.jl \
    CAMPAIGN_FILE TASK_INDEX OUTPUT_DIRECTORY [options]

Options:
  --backend cpu|cuda       (default: cpu)
  --input synthetic|sp2    controlled synthetic MPO or SP2 state (default: synthetic)
  --input-chi N            synthetic input cap; defaults to the case maxdim
  --seed N                 synthetic MPO seed (default: 510578)
  --sp2-steps N            SP2 steps used with --input sp2 (default: 8)
  --warmups N              unrecorded production applies (default: 1)
  --repetitions N          measured repeats on the same input MPO (default: 3)
  --check-tolerance X      relative Hilbert-norm agreement tolerance (default: 1e-11)
"""

using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using Random
using TOML

# Campaigns deliberately remain ordinary Julia files. Load the selected one
# before this file defines compiled functions so the `campaign` binding is in
# the current world age (the same rule used by the existing run diagnostics).
if !isempty(ARGS) && ARGS[1] != "--help"
    Base.include(Main, abspath(ARGS[1]))
end

# Import CUDA at global scope when requested. Loading it inside a compiled
# runtime helper creates a new-world-age binding, which is precisely the
# failure mode avoided by the public campaign launcher.
backend_argument = findfirst(==( "--backend"), ARGS)
if !isnothing(backend_argument) && backend_argument < length(ARGS) &&
   lowercase(ARGS[backend_argument + 1]) == "cuda"
    @eval using CUDA
end

function usage(io::IO=stdout)
    println(io, "Usage: profile_apply_square.jl CAMPAIGN_FILE TASK_INDEX OUTPUT_DIRECTORY [options]")
    println(io, "  --backend cpu|cuda --input synthetic --input-chi 512 --seed 510578")
    println(io, "  --sp2-steps 8 --warmups 1 --repetitions 3 --check-tolerance 1e-11")
end

function parse_arguments(arguments)
    length(arguments) >= 3 || return nothing
    configuration = (
        campaign=abspath(arguments[1]), task=parse(Int, arguments[2]), output=abspath(arguments[3]),
        backend=:cpu, input=:synthetic, input_chi=nothing, seed=510578,
        sp2_steps=8, warmups=1, repetitions=3, check_tolerance=1e-11,
    )
    index = 4
    while index <= length(arguments)
        option = arguments[index]
        option == "--backend" && (configuration = merge(configuration, (backend=Symbol(lowercase(arguments[index + 1])),)); index += 2; continue)
        option == "--input" && (configuration = merge(configuration, (input=Symbol(lowercase(arguments[index + 1])),)); index += 2; continue)
        option == "--input-chi" && (configuration = merge(configuration, (input_chi=parse(Int, arguments[index + 1]),)); index += 2; continue)
        option == "--seed" && (configuration = merge(configuration, (seed=parse(Int, arguments[index + 1]),)); index += 2; continue)
        option == "--sp2-steps" && (configuration = merge(configuration, (sp2_steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        option == "--warmups" && (configuration = merge(configuration, (warmups=parse(Int, arguments[index + 1]),)); index += 2; continue)
        option == "--repetitions" && (configuration = merge(configuration, (repetitions=parse(Int, arguments[index + 1]),)); index += 2; continue)
        option == "--check-tolerance" && (configuration = merge(configuration, (check_tolerance=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        option == "--help" && return nothing
        throw(ArgumentError("unknown option $option"))
    end
    isfile(configuration.campaign) || throw(ArgumentError("campaign file does not exist: $(configuration.campaign)"))
    configuration.backend in (:cpu, :cuda) || throw(ArgumentError("--backend must be cpu or cuda"))
    configuration.input in (:synthetic, :sp2) || throw(ArgumentError("--input must be synthetic or sp2"))
    (isnothing(configuration.input_chi) || configuration.input_chi > 0) ||
        throw(ArgumentError("--input-chi must be positive"))
    configuration.sp2_steps >= 0 || throw(ArgumentError("--sp2-steps must be nonnegative"))
    configuration.warmups >= 0 || throw(ArgumentError("--warmups must be nonnegative"))
    configuration.repetitions > 0 || throw(ArgumentError("--repetitions must be positive"))
    isfinite(configuration.check_tolerance) && configuration.check_tolerance >= 0 ||
        throw(ArgumentError("--check-tolerance must be finite and nonnegative"))
    ispath(configuration.output) && throw(ArgumentError("refusing to overwrite existing output directory: $(configuration.output)"))
    return configuration
end

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

mean_chi(mpo::MPO) = length(mpo) <= 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) / (length(mpo) - 1)

function runtime_callbacks(backend::Symbol)
    backend == :cpu && return (identity, () -> nothing, nothing, "CPU")
    isdefined(Main, :CUDA) || error("CUDA was not loaded; pass --backend cuda at process startup")
    CUDA.functional() || error("CUDA is not functional on this node")
    return (
        value -> ITensors.adapt(CUDA.CuArray, value), CUDA.synchronize, CUDA,
        CUDA.name(CUDA.device()),
    )
end

"""Synchronized phase timing. `value` may contain GPU-resident ITensors."""
function phase_timed(f, synchronize)
    synchronize()
    timed = @timed begin
        value = f()
        synchronize()
        value
    end
    return timed
end

function relative_mpo_error(reference::MPO, candidate::MPO)
    reference_norm_squared = real(inner(reference, reference))
    candidate_norm_squared = real(inner(candidate, candidate))
    overlap = real(inner(reference, candidate))
    numerator_squared = max(0.0, reference_norm_squared + candidate_norm_squared - 2overlap)
    return sqrt(numerator_squared / max(reference_norm_squared, eps(Float64)))
end

"""Exact phase-labelled mirror of the installed MPO/MPO zip-up implementation.

Keep this in lockstep with ITensorMPS `src/mpo.jl`. It is intentionally local
to the benchmark so production calls continue to use `apply` unchanged.
"""
function profiled_zipup_contract(
        A::MPO,
        B::MPO;
        cutoff::Float64,
        maxdim::Int,
        synchronize,
        record_phase,
    )
    N = length(A)
    N == length(B) || throw(DimensionMismatch("MPO lengths do not match"))
    N == 1 && return MPO([A[1] * B[1]])

    A_timing = phase_timed(synchronize) do
        orthogonalize(A, 1)
    end
    A = A_timing.value
    record_phase(:zipup_orthogonalize_left, 0, A_timing)
    B_timing = phase_timed(synchronize) do
        orthogonalize(B, 1)
    end
    B = B_timing.value
    record_phase(:zipup_orthogonalize_right, 0, B_timing)

    A = sim(linkinds, A)
    sA = siteinds(uniqueinds, A, B)
    sB = siteinds(uniqueinds, B, A)
    C = typeof(B)(N)
    lC_i = Index[]
    R = ITensor(true)

    for i in 1:(N - 2)
        contraction = phase_timed(synchronize) do
            R * A[i] * B[i]
        end
        RAB_i = contraction.value
        record_phase(:zipup_local_contract, i, contraction)
        left_inds = [sA[i]..., sB[i]..., lC_i...]
        factorization = phase_timed(synchronize) do
            factorize(
                RAB_i,
                left_inds;
                ortho="left",
                tags=ITensors.commontags(linkinds(A, i)),
                cutoff=cutoff,
                maxdim=maxdim,
                mindim=1,
            )
        end
        C[i], R = factorization.value
        record_phase(:zipup_factorize, i, factorization)
        lC_i = dag(commoninds(C[i], R))
    end

    i = N - 1
    final_contraction = phase_timed(synchronize) do
        R * A[i] * B[i] * A[i + 1] * B[i + 1]
    end
    RAB_i = final_contraction.value
    record_phase(:zipup_final_local_contract, i, final_contraction)
    left_inds = [sA[i]..., sB[i]..., lC_i...]
    final_factorization = phase_timed(synchronize) do
        factorize(
            RAB_i,
            left_inds;
            ortho="right",
            tags=ITensors.commontags(linkinds(A, i)),
            cutoff=cutoff,
            maxdim=maxdim,
            mindim=1,
        )
    end
    C[N - 1], C[N] = final_factorization.value
    record_phase(:zipup_final_factorize, i, final_factorization)

    final_truncation = phase_timed(synchronize) do
        truncate!(C; cutoff=cutoff, maxdim=maxdim, mindim=1)
    end
    record_phase(:final_truncate, 0, final_truncation)
    return C
end

function prepare_sp2_state(case, configuration, to_device, synchronize)
    params = legacy_parameters(case.model, case.representation, case.solver)
    params isa ParametersSquare || throw(ArgumentError("this benchmark currently requires SquareModel"))
    bounds = case.spectral_bounds
    isnothing(bounds) && throw(ArgumentError("the case must provide explicit spectral_bounds"))
    sys = System(params)
    # Match `run_scf!`: construct_rho_0 combines H0, VH, VF, and the identity
    # MPO, so all four must live on the selected backend before initialization.
    # Leaving the three fields on CPU causes a host Matrix to enter a GPU
    # broadcast during the initial affine spectral map.
    sys.H0 = to_device(sys.H0)
    sys.VH = to_device(sys.VH)
    sys.VF = to_device(sys.VF)
    rho = construct_rho_0(sys, params, bounds...; method=:sp2, to_gpu=to_device)
    synchronize()
    Ne = round(Int, (2 ^ params.L) * params.density)
    trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)
    for step in 1:configuration.sp2_steps
        rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        synchronize()
        branch = MPO_MeanField._sp2_branch(
            real(tr(rho)), real(tr(rho_squared)), Ne, trace_tolerance,
        )
        rho = branch == :square ? rho_squared : +(
            2.0 * rho, -rho_squared; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        synchronize()
    end
    return rho, params
end

"""Build a controlled, Hermitian MPO without Hamiltonian or SCF work.

The outer product has link dimension approximately `mps_linkdim^2`. Choosing
`ceil(sqrt(input_chi))` therefore targets the desired MPO bond scale while
leaving the input normalized and reproducible.
"""
function prepare_synthetic_state(params, input_chi::Int, seed::Int, to_device, synchronize)
    sites = siteinds("Qubit", params.L)
    mps_linkdim = ceil(Int, sqrt(input_chi))
    psi = random_mps(MersenneTwister(seed), sites; linkdims=mps_linkdim)
    rho = to_device(outer(psi, dag(prime(psi))))
    synchronize()
    return rho
end

function write_metadata(output, configuration, case, params, device_name, cuda)
    open(joinpath(output, "metadata.toml"), "w") do io
        println(io, "created_at = \"", Dates.now(), "\"")
        println(io, "campaign = \"", replace(configuration.campaign, '\\' => "/"), "\"")
        println(io, "case_label = \"", case.label, "\"")
        println(io, "backend = \"", configuration.backend, "\"")
        println(io, "device_name = \"", device_name, "\"")
        println(io, "cuda_runtime = \"", isnothing(cuda) ? "not_applicable" : cuda.runtime_version(), "\"")
        println(io, "input = \"", configuration.input, "\"")
        println(io, "input_chi_requested = ", something(configuration.input_chi, params.itensors_maxdim))
        println(io, "synthetic_seed = ", configuration.seed)
        println(io, "sp2_preparation_steps = ", configuration.sp2_steps)
        println(io, "warmups = ", configuration.warmups)
        println(io, "repetitions = ", configuration.repetitions)
        println(io, "cutoff = ", params.itensors_tol)
        println(io, "maxdim = ", params.itensors_maxdim)
        println(io, "qtt_length = ", params.L)
        println(io, "matrix_dimension = ", 2 ^ params.L)
        println(io, "check_tolerance = ", configuration.check_tolerance)
        zipup_decomposition = params.itensors_tol <= 1e-12 ? "svd" : "eigen"
        println(io, "zipup_automatic_factorization = \"", zipup_decomposition, "\"")
        println(io, "final_truncate_decomposition = \"svd\"")
        println(io, "implementation = \"production apply versus phase-labelled mirror of installed zipup\"")
        println(io, "phase_timing_note = \"phase rows synchronize after each phase; use production_apply rows for end-to-end comparisons\"")
    end
end

function run_profile(configuration, campaign)
    1 <= configuration.task <= length(campaign.cases) || error("task index is outside campaign")
    case = campaign.cases[configuration.task]
    to_device, synchronize, cuda, device_name = runtime_callbacks(configuration.backend)
    if configuration.input == :sp2
        rho, params = prepare_sp2_state(case, configuration, to_device, synchronize)
    else
        params = legacy_parameters(case.model, case.representation, case.solver)
        params isa ParametersSquare || throw(ArgumentError("synthetic benchmark requires SquareModel"))
        rho = prepare_synthetic_state(
            params, something(configuration.input_chi, params.itensors_maxdim),
            configuration.seed, to_device, synchronize,
        )
    end
    mkpath(configuration.output)
    write_metadata(configuration.output, configuration, case, params, device_name, cuda)

    open(joinpath(configuration.output, "phase_times.csv"), "w") do phases
        write_csv_row(phases, (
            "repetition", "kind", "site", "wall_time_s", "host_allocations_bytes", "host_gc_time_s",
            "rho_max_chi", "rho_mean_chi", "gpu_free_before_bytes", "gpu_free_after_bytes",
            "gpu_used_before_bytes", "gpu_used_after_bytes", "gpu_cached_before_bytes", "gpu_cached_after_bytes",
        ))
        open(joinpath(configuration.output, "summary.csv"), "w") do summary
            write_csv_row(summary, (
                "repetition", "production_apply_time_s", "instrumented_zipup_time_s",
                "production_max_chi", "instrumented_max_chi", "production_mean_chi", "instrumented_mean_chi",
                "relative_mpo_error", "trace_difference", "idempotency_difference",
                "production_host_allocations_bytes", "instrumented_host_allocations_bytes",
            ))
            for _ in 1:configuration.warmups
                warmup = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
                synchronize()
                warmup = nothing
                GC.gc()
            end
            for repetition in 1:configuration.repetitions
                gpu_free_before = isnothing(cuda) ? 0 : cuda.free_memory()
                gpu_used_before = isnothing(cuda) ? 0 : cuda.used_memory()
                gpu_cached_before = isnothing(cuda) ? 0 : cuda.cached_memory()
                production = phase_timed(synchronize) do
                    apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
                end
                production_square = production.value
                gpu_free_after = isnothing(cuda) ? 0 : cuda.free_memory()
                gpu_used_after = isnothing(cuda) ? 0 : cuda.used_memory()
                gpu_cached_after = isnothing(cuda) ? 0 : cuda.cached_memory()
                write_csv_row(phases, (
                    repetition, "production_apply", 0, production.time, production.bytes, production.gctime,
                    maxlinkdim(production_square), mean_chi(production_square), gpu_free_before, gpu_free_after,
                    gpu_used_before, gpu_used_after, gpu_cached_before, gpu_cached_after,
                ))

                phase_total = Ref(0.0)
                phase_allocations = Ref(0)
                record_phase = function (kind, site, timed)
                    phase_total[] += timed.time
                    phase_allocations[] += timed.bytes
                    write_csv_row(phases, (
                        repetition, kind, site, timed.time, timed.bytes, timed.gctime,
                        maxlinkdim(rho), mean_chi(rho), 0, 0, 0, 0, 0, 0,
                    ))
                    flush(phases)
                end
                instrumented = profiled_zipup_contract(
                    rho', rho;
                    cutoff=params.itensors_tol,
                    maxdim=params.itensors_maxdim,
                    synchronize=synchronize,
                    record_phase=record_phase,
                )
                # This is the final `apply` index convention, outside the zip-up contract.
                instrumented = replaceprime(instrumented, 2 => 1)
                synchronize()
                relative_error = relative_mpo_error(production_square, instrumented)
                relative_error <= configuration.check_tolerance || error(
                    "instrumented zip-up differs from production apply: relative error=$relative_error " *
                    "> $(configuration.check_tolerance)",
                )
                production_idempotency = MPO_MeanField.idempotency_residual(rho, production_square)
                instrumented_idempotency = MPO_MeanField.idempotency_residual(rho, instrumented)
                write_csv_row(summary, (
                    repetition, production.time, phase_total[], maxlinkdim(production_square), maxlinkdim(instrumented),
                    mean_chi(production_square), mean_chi(instrumented), relative_error,
                    real(tr(production_square) - tr(instrumented)), production_idempotency - instrumented_idempotency,
                    production.bytes, phase_allocations[],
                ))
                flush(summary)
                production_square = nothing
                instrumented = nothing
                GC.gc()
            end
        end
    end
    println("Apply-profile artefacts written to $(configuration.output)")
end

configuration = parse_arguments(ARGS)
if isnothing(configuration)
    usage()
else
    isdefined(Main, :campaign) || error("campaign file must define `campaign`")
    run_profile(configuration, getfield(Main, :campaign))
end
