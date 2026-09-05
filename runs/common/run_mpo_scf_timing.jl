#!/usr/bin/env julia

"""Run one public MPO-SCF campaign case and record only SCF-loop timings.

Usage:
    julia --project=. runs/common/run_mpo_scf_timing.jl CAMPAIGN.jl TASK OUTPUT cpu|cuda

The final site-density and bond-order measurements are intentionally skipped.
Set `MPO_SCFTIMING_CHECKPOINT=true` to save only the converged QTT/MPO
checkpoint; this does not evaluate real-space observables.
"""

using TOML

length(ARGS) == 4 || error("usage: run_mpo_scf_timing.jl CAMPAIGN.jl TASK OUTPUT cpu|cuda")
campaign_path = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
backend = Symbol(ARGS[4])
backend in (:cpu, :cuda) || error("backend must be cpu or cuda")
isfile(campaign_path) || error("campaign source does not exist: $campaign_path")
ispath(output) && error("refusing to overwrite existing output directory: $output")

if backend == :cuda
    using CUDA
end
using MPO_MeanField

function load_campaign(path::AbstractString)
    namespace = Module(gensym(:MPOHFSCFTimingCampaign))
    Core.eval(namespace, :(using MPO_MeanField))
    Base.include(namespace, path)
    campaign = Base.invokelatest(getfield, namespace, :campaign)
    campaign isa CampaignSpec || error("campaign source must bind `campaign` to a CampaignSpec")
    return campaign
end

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

campaign = load_campaign(campaign_path)
case = campaign_case(campaign, task)
runtime = RuntimeSettings(backend=backend)
timings = NamedTuple[]
save_checkpoint = get(ENV, "MPO_SCFTIMING_CHECKPOINT", "false") == "true"
checkpoint_path = save_checkpoint ? joinpath(output, "converged_state.h5") : nothing

result = Base.invokelatest(
    solve,
    case.model,
    case.representation,
    case.solver;
    runtime=runtime,
    method=:mpo,
    spectral_bounds=case.spectral_bounds,
    verbose=:all,
    checkpoint_path=checkpoint_path,
    measure_observables=false,
    phase_callback=event -> push!(timings, event),
)

mkpath(output)
cp(campaign_path, joinpath(output, "campaign.jl"))

open(joinpath(output, "scf_timing.csv"), "w") do io
    write_csv_row(io, (
        "iteration", "initialization_time_s", "purification_time_s",
        "density_to_host_time_s", "mean_field_time_s", "fields_to_device_time_s",
        "device_diagnostics_time_s", "residuals_time_s", "energy_time_s",
        "mixing_time_s", "scf_iteration_time_s", "purification_iterations",
        "rho_bond_dimension", "hartree_bond_dimension", "fock_bond_dimension",
        "effective_hamiltonian_bond_dimension", "mixing_method",
    ))
    for event in timings
        write_csv_row(io, (
            event.iteration, event.initialization_time_s, event.purification_time_s,
            event.density_to_host_time_s, event.mean_field_time_s,
            event.fields_to_device_time_s, event.device_diagnostics_time_s,
            event.residuals_time_s, event.energy_time_s, event.mixing_time_s,
            event.measured_iteration_time_s, event.purification_iterations,
            event.rho_bond_dimension, event.hartree_bond_dimension,
            event.fock_bond_dimension, event.effective_hamiltonian_bond_dimension,
            event.mixing_method,
        ))
    end
end

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => campaign.name,
        "label" => case.label,
        "task" => task,
        "backend" => string(backend),
        "scf_converged" => result.converged,
        "scf_termination_reason" => string(result.termination_reason),
        "scf_iterations" => length(result.diagnostics.history),
        "solve_elapsed_time_s" => result.elapsed_time_s,
        "post_scf_observables_measured" => false,
        "checkpoint_written" => (save_checkpoint && result.converged),
    ))
end

println("Timing output: $output")
println("SCF: converged=$(result.converged) termination=$(result.termination_reason)")
exit(result.converged ? 0 : 2)
