#!/usr/bin/env julia

"""Targeted reciprocal-space survey of an MPO/SP2 checkpoint.

The calculation makes one QTT Fourier transform and then queries selected
Fourier coefficients. It never materializes the N-point density or spectrum.
Predicted quasiperiodic satellites are compared with nearby same-line,
off-line, and deterministic generic controls, so a visible satellite is not
mistaken for an open-boundary/background artifact.
"""

using Dates
using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using TOML

length(ARGS) == 4 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY",
)
campaign_path = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
result_directory = abspath(ARGS[3])
output = abspath(ARGS[4])
checkpoint_path = let candidates = (
    joinpath(result_directory, "converged_state.h5"),
    joinpath(result_directory, "final_fixed_sp2_state.h5"),
)
    found = findfirst(isfile, candidates)
    isnothing(found) && error(
        "checkpoint does not exist: expected converged_state.h5 or final_fixed_sp2_state.h5 in $result_directory",
    )
    candidates[found]
end
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:QTTReciprocalCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
campaign = Base.invokelatest(getfield, namespace, :campaign)
case, params = if hasproperty(campaign, :runs)
    legacy_case = campaign.runs[task_index]
    (legacy_case, legacy_case.params)
else
    public_case = campaign_case(campaign, task_index)
    (public_case, Base.invokelatest(
        legacy_parameters, public_case.model, public_case.representation, public_case.solver,
    ))
end
params isa ParametersSquare || error("reciprocal survey requires ParametersSquare")
checkpoint_maxdim = h5open(checkpoint_path, "r") do file
    Int(read(file, "itensors_maxdim"))
end
params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
    itensors_maxdim=checkpoint_maxdim, density=params.density,
    purification_steps=params.purification_steps,
    scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
    scf_max_iterations=params.scf_max_iterations,
)
system = read_mpo_checkpoint(checkpoint_path, params).system
side_bits = params.L ÷ 2
side = 1 << side_bits
carrier = side ÷ 2
tau = sqrt(2.0) - 5.0 / 6.0
cutoff = params.itensors_tol
maxdim = params.itensors_maxdim

csv(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_row(io, values) = println(io, join(csv.(values), ','))
wrap(index) = mod(index, side)

println("Extracting QTT density and Fourier transform for $(side)x$(side)")
density = density_diagonal_mps(system)
centered, mean_density, trace_value = centered_density_mps(
    density, system.sites; cutoff=cutoff, maxdim=maxdim,
)
fourier_timing = @timed qtt_fourier_square(
    centered, system.sites; sign=-1, cutoff_MPO=cutoff, cutoff=cutoff, maxdim=maxdim,
)
fourier = fourier_timing.value

rows = NamedTuple[]
function sample!(family, label, kx, ky; m=missing)
    value = qtt_fourier_amplitude(fourier, wrap(kx), wrap(ky))
    push!(rows, (
        family=family, label=label, m=m, kx=wrap(kx), ky=wrap(ky),
        momentum_x=2pi * wrap(kx) / side, momentum_y=2pi * wrap(ky) / side,
        real=real(value), imag=imag(value), power=abs2(value),
    ))
end

sample!("carrier", "checkerboard", carrier, carrier; m=0)
satellite_indices = Int[]
for m in -5:5
    index = wrap(round(Int, side * mod(pi + 2pi * tau * m, 2pi) / (2pi)))
    push!(satellite_indices, index)
    m == 0 && continue
    # Predicted modulation wave vectors along the two separable axes.
    sample!("satellite", "horizontal", index, carrier; m=m)
    sample!("satellite", "vertical", carrier, index; m=m)
    # Matched controls on the same reciprocal line but displaced by 17 bins.
    # The shift is much larger than a single-bin peak width and is checked not
    # to overlap any predicted sample for this finite lattice.
    control = wrap(index + 17)
    control in satellite_indices && (control = wrap(index + 29))
    sample!("line_control", "horizontal_offset", control, carrier; m=m)
    sample!("line_control", "vertical_offset", carrier, control; m=m)
    # Same k component, displaced off the expected modulation line.
    sample!("off_line_control", "horizontal", index, carrier + 17; m=m)
    sample!("off_line_control", "vertical", carrier + 17, index; m=m)
end
for r in 1:32
    # Deterministic generic reciprocal controls, avoiding the carrier lines.
    kx = wrap(97r + 11)
    ky = wrap(193r + 29)
    (kx == carrier || ky == carrier) && continue
    sample!("generic_control", "deterministic_$r", kx, ky)
end

power_values(family) = [row.power for row in rows if row.family == family]
carrier_power = only(power_values("carrier"))
summary = Dict(
    "created_at" => string(now()),
    "diagnostic" => "qtt_targeted_reciprocal_space_survey",
    "source_checkpoint" => checkpoint_path,
    "side" => side,
    "matrix_dimension" => side^2,
    "cutoff" => cutoff,
    "maxdim" => maxdim,
    "mean_density" => real(mean_density),
    "trace" => real(trace_value),
    "carrier_power" => carrier_power,
    "satellite_power_max" => maximum(power_values("satellite")),
    "satellite_power_mean" => sum(power_values("satellite")) / length(power_values("satellite")),
    "line_control_power_max" => maximum(power_values("line_control")),
    "line_control_power_mean" => sum(power_values("line_control")) / length(power_values("line_control")),
    "off_line_control_power_max" => maximum(power_values("off_line_control")),
    "off_line_control_power_mean" => sum(power_values("off_line_control")) / length(power_values("off_line_control")),
    "generic_control_power_max" => maximum(power_values("generic_control")),
    "generic_control_power_mean" => sum(power_values("generic_control")) / length(power_values("generic_control")),
    "satellite_to_carrier_power_ratio" => maximum(power_values("satellite")) / carrier_power,
    "satellite_to_line_control_ratio" => maximum(power_values("satellite")) / maximum(power_values("line_control")),
    "satellite_to_off_line_control_ratio" => maximum(power_values("satellite")) / maximum(power_values("off_line_control")),
    "fourier_max_chi" => maxlinkdim(fourier),
    "fourier_time_s" => fourier_timing.time,
)

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
try
    open(joinpath(temporary, "summary.toml"), "w") do io
        TOML.print(io, summary)
    end
    open(joinpath(temporary, "samples.csv"), "w") do io
        write_row(io, propertynames(first(rows)))
        for row in rows
            write_row(io, Tuple(row))
        end
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end
println("QTT reciprocal survey complete: $output")
