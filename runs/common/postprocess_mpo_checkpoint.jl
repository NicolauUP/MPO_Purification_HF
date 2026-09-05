#!/usr/bin/env julia

"""Compute full observables from a saved public MPO checkpoint.

Usage:
  postprocess_mpo_checkpoint.jl CAMPAIGN.jl TASK SOURCE_ROOT OUTPUT_ROOT

The source result is resolved from SOURCE_ROOT/campaign/task_XXXX_label. SCF
is not rerun. This postprocessor intentionally computes the complete site and
nearest-neighbour bond arrays, energy, Hermiticity, idempotency, and optional
stationarity diagnostics in a separate non-overwriting output directory.
"""

using Dates
using HDF5
using ITensors
using ITensorMPS
using MPO_MeanField
using TOML

length(ARGS) == 4 || error("usage: CAMPAIGN.jl TASK SOURCE_ROOT OUTPUT_ROOT")
campaign_path = abspath(ARGS[1])
task = parse(Int, ARGS[2])
source_root = abspath(ARGS[3])
output_root = abspath(ARGS[4])
isfile(campaign_path) || error("campaign does not exist: $campaign_path")

include(campaign_path)
@isdefined(campaign) || error("campaign source must define `campaign`")
case = campaign_case(campaign, task)
case.model isa SquareModel || error("postprocessor currently requires SquareModel")

safe_component(text) = begin
    value = replace(String(text), r"[^A-Za-z0-9._-]+" => "_")
    isempty(value) ? error("empty path component") : value
end
source_directory = joinpath(source_root, safe_component(campaign.name),
    "task_" * lpad(string(task), 4, '0') * "_" * safe_component(case.label))
checkpoint = joinpath(source_directory, "converged_state.h5")
isfile(checkpoint) || error("checkpoint does not exist: $checkpoint")

params = legacy_parameters(case.model, case.representation, case.solver)
checkpoint_maxdim = h5open(checkpoint, "r") do file
    Int(read(file, "itensors_maxdim"))
end
if checkpoint_maxdim != params.itensors_maxdim
    params = ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
        itensors_maxdim=checkpoint_maxdim, density=params.density,
        purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end
loaded = read_mpo_checkpoint(checkpoint, params)
system = loaded.system
println("postprocessing task=$task label=$(case.label) side=$(2^div(params.L,2))")
started = time()
values = observables_square(system; measure_stationarity=false)
elapsed = time() - started

output = joinpath(output_root, safe_component(campaign.name),
    "task_" * lpad(string(task), 4, '0') * "_" * safe_component(case.label))
ispath(output) && error("refusing to overwrite existing output: $output")
mkpath(output)
open(joinpath(output, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "available" => true,
        "source_directory" => source_directory,
        "particle_number" => values.particle_number,
        "energy_kinetic" => values.energy.kinetic,
        "energy_hartree" => values.energy.hartree,
        "energy_fock" => values.energy.fock,
        "energy_interaction" => values.energy.interaction,
        "energy_total" => values.energy.total,
        "hermiticity_residual" => values.hermiticity_residual,
        "idempotency_residual" => values.idempotency_residual,
        "stationarity_residual" => values.stationarity_residual,
        "measurement_time_s" => elapsed,
    ))
end
open(joinpath(output, "site_density.csv"), "w") do io
    println(io, "site,density")
    for (site, value) in enumerate(values.site_density)
        println(io, site, ",", value)
    end
end
open(joinpath(output, "bond_order.csv"), "w") do io
    println(io, "site_left,site_right,orientation,real,imag")
    for (bond, value) in zip(values.horizontal_bonds, values.horizontal_bond_order)
        println(io, bond[1], ",", bond[2], ",horizontal,", real(value), ",", imag(value))
    end
    for (bond, value) in zip(values.vertical_bonds, values.vertical_bond_order)
        println(io, bond[1], ",", bond[2], ",vertical,", real(value), ",", imag(value))
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "format_version" => 1,
        "created_at_utc" => string(now(UTC)),
        "campaign" => campaign.name,
        "campaign_file" => campaign_path,
        "task" => task,
        "label" => case.label,
        "source_directory" => source_directory,
        "checkpoint" => checkpoint,
        "maxdim" => checkpoint_maxdim,
        "postprocessor" => "observables_square_from_mpo_checkpoint",
        "scf_converged" => loaded.converged,
        "scf_termination_reason" => string(loaded.termination_reason),
        "scf_iterations" => loaded.scf_iterations,
    ))
end
cp(joinpath(source_directory, "scf_history.csv"), joinpath(output, "scf_history.csv"))
println("postprocessed result: $output")
