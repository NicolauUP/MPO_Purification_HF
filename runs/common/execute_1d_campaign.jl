#!/usr/bin/env julia

"""
Execute one explicit 1D campaign point.

Usage:
    julia --project=. runs/common/execute_1d_campaign.jl CAMPAIGN_FILE TASK_INDEX

`CAMPAIGN_FILE` must define `campaign`, a named tuple with `name` and `runs`.
Every entry of `runs` must provide `label`, `params`, `spectral_bounds`, and
`purification_method`.
"""

using Dates
using LinearAlgebra
using Printf
using SHA
using TOML
using ITensorMPS
using MPO_MeanField

length(ARGS) == 2 || error("usage: execute_1d_campaign.jl CAMPAIGN_FILE TASK_INDEX")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
isnothing(task_index) && error("TASK_INDEX must be an integer, got $(ARGS[2])")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
hasproperty(campaign, :name) || error("campaign must provide `name`")
hasproperty(campaign, :runs) || error("campaign must provide `runs`")
1 <= task_index <= length(campaign.runs) || error(
    "TASK_INDEX=$task_index is outside 1:$(length(campaign.runs))",
)

spec = campaign.runs[task_index]
for field in (:label, :params, :spectral_bounds, :purification_method)
    hasproperty(spec, field) || error("campaign run $task_index is missing `$field`")
end
spec.params isa Parameters1D || error("execute_1d_campaign.jl requires Parameters1D")
spec.spectral_bounds isa Tuple && length(spec.spectral_bounds) == 2 || error(
    "spectral_bounds must be a two-element tuple",
)

function csv_escape(value)
    text = string(value)
    return '"' * replace(text, '"' => "\"\"") * '"'
end

function write_csv_row(io, values)
    println(io, join(csv_escape.(values), ','))
end

function safe_component(text)
    component = replace(string(text), r"[^A-Za-z0-9._-]+" => "_")
    isempty(component) && error("run label must contain at least one safe character")
    return component
end

function git_revision(repo_root)
    try
        return readchomp(`git -C $repo_root rev-parse HEAD`)
    catch
        return "unavailable"
    end
end

function sha1_file(path)
    return bytes2hex(sha1(read(path)))
end

results_root = get(ENV, "MPO_RESULTS_ROOT", "")
isempty(results_root) && error("set MPO_RESULTS_ROOT to an external result directory")
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
run_name = @sprintf("task_%04d_%s", task_index, safe_component(spec.label))
run_dir = joinpath(results_root, safe_component(campaign.name), run_name)
ispath(run_dir) && error("refusing to overwrite existing result directory: $run_dir")
mkpath(run_dir)

cp(campaign_file, joinpath(run_dir, "input.jl"))
selection = Dict(
    "campaign" => string(campaign.name),
    "task_index" => task_index,
    "label" => string(spec.label),
    "purification_method" => string(spec.purification_method),
    "spectral_bounds" => collect(Float64.(spec.spectral_bounds)),
)
open(joinpath(run_dir, "selection.toml"), "w") do io
    TOML.print(io, selection)
end

metadata = Dict(
    "campaign" => string(campaign.name),
    "task_index" => task_index,
    "label" => string(spec.label),
    "started_at" => string(now(UTC)),
    "julia_version" => string(VERSION),
    "git_revision" => git_revision(repo_root),
    "project_sha1" => sha1_file(joinpath(repo_root, "Project.toml")),
    "manifest_sha1" => sha1_file(joinpath(repo_root, "Manifest.toml")),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", string(task_index)),
    "threads" => Threads.nthreads(),
)

sys = System(spec.params)
converged = run_scf!(sys, Float64(spec.spectral_bounds[1]), Float64(spec.spectral_bounds[2]);
    purification_method=spec.purification_method,
    verify_spectral_bounds=get(spec, :verify_spectral_bounds, false),
    spectral_safety_margin=get(spec, :spectral_safety_margin, 0.0),
    chemical_potential=get(spec, :chemical_potential, nothing),
    record_energy=true,
    verbose=get(spec, :verbose, :nothing),
)

diagnostics = scf_diagnostics(sys)
obs = observables_1d(sys)

open(joinpath(run_dir, "scf_history.csv"), "w") do io
    write_csv_row(io, (
        "iteration", "trace", "vh_residual", "vf_residual", "rho_residual",
        "commutator_residual", "two_cycle_residual", "purification_converged",
        "purification_termination_reason", "purification_iterations", "energy_total",
    ))
    for record in diagnostics.history
        write_csv_row(io, (
            record.iteration, record.trace, record.vh_residual, record.vf_residual,
            record.rho_residual, record.commutator_residual, record.two_cycle_residual,
            record.purification_converged, record.purification_termination_reason,
            record.purification_iterations, record.energy_total,
        ))
    end
end

open(joinpath(run_dir, "site_density.csv"), "w") do io
    write_csv_row(io, ("site", "density"))
    for (site, density) in enumerate(obs.site_density)
        write_csv_row(io, (site, density))
    end
end

open(joinpath(run_dir, "bond_order.csv"), "w") do io
    write_csv_row(io, ("site_left", "site_right", "real", "imag"))
    for (site, value) in enumerate(obs.bond_order)
        write_csv_row(io, (site, site + 1, real(value), imag(value)))
    end
end

observables_summary = Dict(
    "particle_number" => obs.particle_number,
    "energy_kinetic" => obs.energy.kinetic,
    "energy_hartree" => obs.energy.hartree,
    "energy_fock" => obs.energy.fock,
    "energy_interaction" => obs.energy.interaction,
    "energy_total" => obs.energy.total,
    "hermiticity_residual" => obs.hermiticity_residual,
    "idempotency_residual" => obs.idempotency_residual,
    "stationarity_residual" => obs.stationarity_residual,
)
open(joinpath(run_dir, "observables.toml"), "w") do io
    TOML.print(io, observables_summary)
end

metadata["finished_at"] = string(now(UTC))
metadata["scf_converged"] = converged
metadata["scf_termination_reason"] = string(diagnostics.termination_reason)
metadata["scf_iterations"] = length(diagnostics.history)
metadata["final_bond_dimension"] = maxlinkdim(sys.ρ)
open(joinpath(run_dir, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

println("Result directory: $run_dir")
println("SCF: converged=$converged termination=$(diagnostics.termination_reason)")
exit(converged ? 0 : 2)
