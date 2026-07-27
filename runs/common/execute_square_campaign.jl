#!/usr/bin/env julia

"""Execute one explicit square-lattice MPO Hartree--Fock campaign point."""

using Dates
using Printf
using SHA
using TOML
using ITensorMPS
using MPO_MeanField

length(ARGS) == 2 || error("usage: execute_square_campaign.jl CAMPAIGN_FILE TASK_INDEX")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
isnothing(task_index) && error("TASK_INDEX must be an integer")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
for field in (:label, :params, :spectral_bounds, :purification_method)
    hasproperty(spec, field) || error("campaign run $task_index is missing `$field`")
end
params = spec.params
params isa ParametersSquare || error("square campaign requires ParametersSquare")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
safe_component(value) = begin
    result = replace(string(value), r"[^A-Za-z0-9._-]+" => "_")
    isempty(result) && error("unsafe empty path component")
    result
end
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end

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
    "mcweeny_form" => string(get(spec, :mcweeny_form, :standard)),
    "allow_unconverged_purification" => get(
        spec, :allow_unconverged_purification, false,
    ),
)
for field in (:chemical_potential, :mcweeny_trace_target, :mcweeny_trace_tolerance)
    value = get(spec, field, nothing)
    isnothing(value) || (selection[string(field)] = value)
end
open(joinpath(run_dir, "selection.toml"), "w") do io
    TOML.print(io, selection)
end
metadata = Dict(
    "campaign" => string(campaign.name), "task_index" => task_index,
    "label" => string(spec.label), "started_at" => string(now(UTC)),
    "julia_version" => string(VERSION), "git_revision" => git_revision(repo_root),
    "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
    "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", string(task_index)),
    "threads" => Threads.nthreads(),
)

sys = System(params)
progress_path = joinpath(run_dir, "progress.txt")
converged = open(progress_path, "w") do progress
    run_scf!(sys, Float64(spec.spectral_bounds[1]), Float64(spec.spectral_bounds[2]);
        purification_method=spec.purification_method,
        chemical_potential=get(spec, :chemical_potential, nothing),
        mcweeny_form=get(spec, :mcweeny_form, :standard),
        mcweeny_trace_target=get(spec, :mcweeny_trace_target, nothing),
        mcweeny_trace_tolerance=get(spec, :mcweeny_trace_tolerance, nothing),
        allow_unconverged_purification=get(
            spec, :allow_unconverged_purification, false,
        ),
        verify_spectral_bounds=get(spec, :verify_spectral_bounds, false),
        record_energy=true, verbose=get(spec, :verbose, :all), io=progress,
        overwrite_progress=false,
    )
end
diagnostics = scf_diagnostics(sys)
obs = observables_square(sys)

open(joinpath(run_dir, "scf_history.csv"), "w") do io
    write_csv_row(io, ("iteration", "trace", "vh_residual", "vf_residual", "rho_residual", "commutator_residual", "two_cycle_residual", "purification_converged", "purification_termination_reason", "purification_iterations", "purification_selected_iteration", "rho_bond_dimension", "hartree_bond_dimension", "fock_bond_dimension", "effective_hamiltonian_bond_dimension", "energy_total"))
    for record in diagnostics.history
        write_csv_row(io, (record.iteration, record.trace, record.vh_residual, record.vf_residual, record.rho_residual, record.commutator_residual, record.two_cycle_residual, record.purification_converged, record.purification_termination_reason, record.purification_iterations, record.purification_selected_iteration, record.rho_bond_dimension, record.hartree_bond_dimension, record.fock_bond_dimension, record.effective_hamiltonian_bond_dimension, record.energy_total))
    end
end
open(joinpath(run_dir, "site_density.csv"), "w") do io
    write_csv_row(io, ("site", "density"))
    for (site, density) in enumerate(obs.site_density)
        write_csv_row(io, (site, density))
    end
end
open(joinpath(run_dir, "bond_order.csv"), "w") do io
    write_csv_row(io, ("site_left", "site_right", "orientation", "real", "imag"))
    for (bond, value) in zip(obs.horizontal_bonds, obs.horizontal_bond_order)
        write_csv_row(io, (bond[1], bond[2], "horizontal", real(value), imag(value)))
    end
    for (bond, value) in zip(obs.vertical_bonds, obs.vertical_bond_order)
        write_csv_row(io, (bond[1], bond[2], "vertical", real(value), imag(value)))
    end
end
open(joinpath(run_dir, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "particle_number" => obs.particle_number,
        "energy_kinetic" => obs.energy.kinetic, "energy_hartree" => obs.energy.hartree,
        "energy_fock" => obs.energy.fock, "energy_interaction" => obs.energy.interaction,
        "energy_total" => obs.energy.total, "hermiticity_residual" => obs.hermiticity_residual,
        "idempotency_residual" => obs.idempotency_residual,
        "stationarity_residual" => obs.stationarity_residual,
    ))
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
