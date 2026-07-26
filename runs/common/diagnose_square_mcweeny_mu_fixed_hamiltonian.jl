#!/usr/bin/env julia

"""Fixed-Hamiltonian direct-McWeeny control for a square campaign.

This uses the production `:mcweeny_mu` initialization and recurrence on the
first SCF Hamiltonian H0 + S. Unlike canonical SP2 it has no trace-selected
branch: P -> 3P^2 - 2P^3. Trace is recorded as an invariant only.
"""

using Dates
using ITensors
using ITensorMPS
using MPO_MeanField
using SHA
using TOML

length(ARGS) in (6, 7, 8) || error("usage: diagnose_square_mcweeny_mu_fixed_hamiltonian.jl CAMPAIGN_FILE TASK_INDEX MAXDIM OUTPUT_DIRECTORY SPECTRAL_LOWER SPECTRAL_UPPER [MU [ITENSORS_TOL]]")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
maxdim = tryparse(Int, ARGS[3])
output = abspath(ARGS[4])
lower = tryparse(Float64, ARGS[5])
upper = tryparse(Float64, ARGS[6])
chemical_potential = length(ARGS) >= 7 ? tryparse(Float64, ARGS[7]) : 0.0
itensors_tol_override = length(ARGS) == 8 ? tryparse(Float64, ARGS[8]) : nothing
isnothing(task_index) && error("TASK_INDEX must be an integer")
isnothing(maxdim) && error("MAXDIM must be an integer")
isnothing(lower) && error("SPECTRAL_LOWER must be a float")
isnothing(upper) && error("SPECTRAL_UPPER must be a float")
isnothing(chemical_potential) && error("MU must be a float")
length(ARGS) == 8 && isnothing(itensors_tol_override) && error("ITENSORS_TOL must be a float")
maxdim > 0 || error("MAXDIM must be positive")
isfinite(lower) && isfinite(upper) && lower < upper || error("spectral bounds must be finite and strictly ordered")
lower < chemical_potential < upper || error("MU must lie strictly within the spectral interval")
isnothing(itensors_tol_override) || (isfinite(itensors_tol_override) && itensors_tol_override > 0) || error("ITENSORS_TOL must be finite and positive")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("fixed square McWeeny diagnostic requires ParametersSquare")
ispath(output) && error("refusing to overwrite existing output directory: $output")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end

function with_truncation(params::ParametersSquare, maxdim::Int, itensors_tol::Float64)
    ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=itensors_tol,
        itensors_maxdim=maxdim, density=params.density,
        purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end

function mean_bond_dimension(mpo::MPO)
    length(mpo) <= 1 && return 1.0
    dimensions = [dim(commonind(mpo[index], mpo[index + 1])) for index in 1:(length(mpo) - 1)]
    isempty(dimensions) ? 1.0 : sum(dimensions) / length(dimensions)
end

params = with_truncation(spec.params, maxdim, something(itensors_tol_override, spec.params.itensors_tol))
bounds = (lower, upper)
campaign_bounds = (Float64(spec.spectral_bounds[1]), Float64(spec.spectral_bounds[2]))
N = 2 ^ params.L
Ne = round(Int, N * params.density)
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
started_at = now(UTC)
mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))

sys = System(params)
initialization = @timed construct_rho_0(
    sys, params, bounds[1], bounds[2]; method=:mcweeny_mu,
    chemical_potential=chemical_potential,
)
rho = initialization.value
idempotency_tolerance = 1e-6
converged = false
last_iteration = 0

open(joinpath(output, "iterations.csv"), "w") do io
    write_csv_row(io, (
        "iteration", "trace", "trace_error", "idempotency_residual",
        "hermiticity_residual", "rho_max_chi", "rho_squared_max_chi",
        "rho_cubed_max_chi", "rho_mean_chi", "rho_squared_mean_chi",
        "rho_cubed_mean_chi", "cap_reached", "step_time_s",
        "step_allocations_bytes", "step_gc_time_s",
    ))
    for iteration in 1:params.purification_steps
        global rho, converged, last_iteration
        step = @timed begin
            rho_squared = apply(rho, rho;
                cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
            )
            trace_value = real(tr(rho))
            idempotency = MPO_MeanField.idempotency_residual(rho, rho_squared)
            hermiticity = MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
            if idempotency < idempotency_tolerance
                (rho_squared, nothing, trace_value, idempotency, hermiticity)
            else
                rho_cubed = apply(rho, rho_squared;
                    cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
                )
                (rho_squared, rho_cubed, trace_value, idempotency, hermiticity)
            end
        end
        rho_squared, rho_cubed, trace_value, idempotency, hermiticity = step.value
        rho_chi = maxlinkdim(rho)
        rho_squared_chi = maxlinkdim(rho_squared)
        rho_cubed_chi = isnothing(rho_cubed) ? 0 : maxlinkdim(rho_cubed)
        write_csv_row(io, (
            iteration, trace_value, abs(trace_value - Ne), idempotency, hermiticity,
            rho_chi, rho_squared_chi, rho_cubed_chi,
            mean_bond_dimension(rho), mean_bond_dimension(rho_squared),
            isnothing(rho_cubed) ? 0.0 : mean_bond_dimension(rho_cubed),
            rho_chi >= maxdim || rho_squared_chi >= maxdim || rho_cubed_chi >= maxdim,
            step.time, step.bytes, step.gctime,
        ))
        flush(io)
        last_iteration = iteration
        if idempotency < idempotency_tolerance
            converged = true
            break
        end
        rho = +(
            3.0 * rho_squared, -2.0 * rho_cubed;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
    end
end

final_measurement = @timed begin
    rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    (
        trace=real(tr(rho)),
        idempotency=MPO_MeanField.idempotency_residual(rho, rho_squared),
        hermiticity=MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params),
        rho_squared_chi=maxlinkdim(rho_squared),
    )
end
final = final_measurement.value
termination_reason = converged ? :idempotency_threshold : :max_iterations

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => string(campaign.name), "label" => string(spec.label),
        "task_index" => task_index, "diagnostic" => "fixed_initial_hamiltonian_mcweeny_mu",
        "matrix_dimension" => N, "target_particles" => Ne,
        "chemical_potential" => chemical_potential,
        "spectral_lower" => bounds[1], "spectral_upper" => bounds[2],
        "campaign_spectral_lower" => campaign_bounds[1],
        "campaign_spectral_upper" => campaign_bounds[2],
        "itensors_tol" => params.itensors_tol, "itensors_maxdim" => maxdim,
        "campaign_itensors_tol" => spec.params.itensors_tol,
        "purification_steps" => params.purification_steps,
        "idempotency_tolerance" => idempotency_tolerance,
        "initial_rho_max_chi" => maxlinkdim(initialization.value),
        "initial_rho_mean_chi" => mean_bond_dimension(initialization.value),
        "initialization_time_s" => initialization.time,
        "initialization_allocations_bytes" => initialization.bytes,
        "initialization_gc_time_s" => initialization.gctime,
        "converged" => converged, "termination_reason" => string(termination_reason),
        "iterations" => last_iteration, "final_rho_max_chi" => maxlinkdim(rho),
        "final_rho_squared_max_chi" => final.rho_squared_chi,
        "final_rho_mean_chi" => mean_bond_dimension(rho),
        "final_trace" => final.trace, "final_trace_error" => abs(final.trace - Ne),
        "final_idempotency_residual" => final.idempotency,
        "final_hermiticity_residual" => final.hermiticity,
        "final_measurement_time_s" => final_measurement.time,
        "final_measurement_allocations_bytes" => final_measurement.bytes,
        "final_measurement_gc_time_s" => final_measurement.gctime,
    ))
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "started_at" => string(started_at), "finished_at" => string(now(UTC)),
        "julia_version" => string(VERSION), "threads" => Threads.nthreads(),
        "git_revision" => git_revision(repo_root),
        "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
        "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
        "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", "local"),
    ))
end

println("Fixed-Hamiltonian McWeeny diagnostic: label=$(spec.label) maxdim=$maxdim mu=$chemical_potential")
println("converged=$converged termination=$termination_reason iterations=$last_iteration")
println("final trace error=$(abs(final.trace - Ne)) idempotency=$(final.idempotency) chi=$(maxlinkdim(rho))")
println("Result directory: $output")
