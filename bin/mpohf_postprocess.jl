#!/usr/bin/env julia

"""Measure final host-side observables from a saved MPO checkpoint.

Usage:
    julia --project=. bin/mpohf_postprocess.jl CAMPAIGN.jl TASK RESULT_DIRECTORY
"""

using Dates
using MPO_MeanField
using TOML

length(ARGS) == 3 || error(
    "usage: mpohf_postprocess.jl CAMPAIGN.jl TASK RESULT_DIRECTORY",
)
campaign_path = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
result_directory = abspath(ARGS[3])
isdir(result_directory) || error("result directory does not exist: $result_directory")
checkpoint_path = joinpath(result_directory, "converged_state.h5")
isfile(checkpoint_path) || error("checkpoint does not exist: $checkpoint_path")
output = joinpath(result_directory, "postprocess")
ispath(output) && error("refusing to overwrite postprocess output: $output")

namespace = Module(gensym(:MPOHFPostprocessCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
isdefined(namespace, :campaign) || error("campaign source must define campaign")
campaign = Base.invokelatest(getfield, namespace, :campaign)
case = campaign_case(campaign, task_index)
params = Base.invokelatest(
    legacy_parameters, case.model, case.representation, case.solver,
)
loaded = read_mpo_checkpoint(checkpoint_path, params)

println("Measuring host-side observables from $checkpoint_path")
flush(stdout)
measurement = @timed observables(loaded.system)
measured = measurement.value

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
csv(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_row(io, values) = println(io, join(csv.(values), ','))
try
    open(joinpath(temporary, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "created_at" => string(Dates.now()),
            "campaign" => campaign_path,
            "task_index" => task_index,
            "label" => case.label,
            "source_checkpoint" => checkpoint_path,
            "checkpoint_format_version" => loaded.format_version,
            "checkpoint_scf_converged" => loaded.converged,
            "checkpoint_scf_termination_reason" => string(loaded.termination_reason),
            "checkpoint_scf_iterations" => loaded.scf_iterations,
            "measurement_time_s" => measurement.time,
            "measurement_allocations_bytes" => measurement.bytes,
            "measurement_gc_time_s" => measurement.gctime,
        ))
    end
    open(joinpath(temporary, "observables.toml"), "w") do io
        TOML.print(io, Dict(
            "particle_number" => measured.particle_number,
            "energy_kinetic" => measured.energy.kinetic,
            "energy_hartree" => measured.energy.hartree,
            "energy_fock" => measured.energy.fock,
            "energy_interaction" => measured.energy.interaction,
            "energy_total" => measured.energy.total,
            "hermiticity_residual" => measured.hermiticity_residual,
            "idempotency_residual" => measured.idempotency_residual,
            "stationarity_residual" => measured.stationarity_residual,
        ))
    end
    open(joinpath(temporary, "site_density.csv"), "w") do io
        write_row(io, ("site", "density"))
        for (site, density) in enumerate(measured.site_density)
            write_row(io, (site, density))
        end
    end
    open(joinpath(temporary, "bond_order.csv"), "w") do io
        write_row(io, ("site_left", "site_right", "orientation", "real", "imag"))
        if case.model isa ChainModel
            for (site, value) in enumerate(measured.bond_order)
                write_row(io, (site, site + 1, "chain", real(value), imag(value)))
            end
        else
            for (bond, value) in zip(measured.horizontal_bonds, measured.horizontal_bond_order)
                write_row(io, (bond[1], bond[2], "horizontal", real(value), imag(value)))
            end
            for (bond, value) in zip(measured.vertical_bonds, measured.vertical_bond_order)
                write_row(io, (bond[1], bond[2], "vertical", real(value), imag(value)))
            end
        end
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end

println("Post-processing complete: $output")
println("measurement_time_s=$(measurement.time)")
flush(stdout)
