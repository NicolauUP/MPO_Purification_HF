#!/usr/bin/env julia

"""Analyze a converged square SP2/MPO checkpoint with QTT FFT and D2.

This is intentionally independent of the exhaustive observables postprocessor:
it never enumerates all sites or bonds and never squares the full projector.
Two analysis compression policies provide an internal convergence check.
"""

using Dates
using ITensors, ITensorMPS
using HDF5
using LinearAlgebra
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
isfile(checkpoint_path) || error("checkpoint does not exist: $checkpoint_path")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:QTTChargeCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
isdefined(namespace, :campaign) || error("campaign source must define campaign")
campaign = Base.invokelatest(getfield, namespace, :campaign)
# Both the public CampaignSpec and the established fixed-H diagnostic
# campaigns are accepted. The latter is deliberately reported as a fixed
# initial-field state below, never as an SCF solution.
case, params, label = if hasproperty(campaign, :runs)
    legacy_case = campaign.runs[task_index]
    (legacy_case, legacy_case.params, string(legacy_case.label))
else
    public_case = campaign_case(campaign, task_index)
    (public_case, Base.invokelatest(
        legacy_parameters, public_case.model, public_case.representation, public_case.solver,
    ), string(public_case.label))
end
params isa ParametersSquare || error("QTT charge analysis currently requires a square model")
checkpoint_maxdim = h5open(checkpoint_path, "r") do file
    Int(read(file, "itensors_maxdim"))
end
if checkpoint_maxdim != params.itensors_maxdim
    # Fixed-H cap ladders share a physical campaign but intentionally alter
    # only the representation cap. Rebuild that exact representation here so
    # checkpoint validation remains strict for every rung.
    params = ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
        itensors_maxdim=checkpoint_maxdim, density=params.density,
        purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end
loaded = read_mpo_checkpoint(checkpoint_path, params)
system = loaded.system
state_kind = basename(checkpoint_path) == "final_fixed_sp2_state.h5" ?
    "fixed_initial_field_sp2" : "self_consistent_hf"
levels = length(system.sites)
side_bits = levels ÷ 2
side = 1 << side_bits
N = side^2

policies = [
    (name="loose", cutoff=1e-8, maxdim=min(256, params.itensors_maxdim)),
    (name="tight", cutoff=params.itensors_tol, maxdim=params.itensors_maxdim),
]
fit_scales = side_bits >= 5 ? collect(2:(side_bits - 1)) : collect(1:side_bits)

csv(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_row(io, values) = println(io, join(csv.(values), ','))

function hann_window_mps(sites, tolerance)
    window(z) = begin
        x, y = square_lattice_decoder(Int(z), levels)
        wx = 0.5 * (1.0 - cos(2pi * x / (side - 1)))
        wy = 0.5 * (1.0 - cos(2pi * y / (side - 1)))
        wx * wy
    end
    _, _, state = MPO_MeanField.Quantics_TCI(window, Float64, sites, tolerance)
    return state
end

function quasiperiodic_samples(transformed, tau)
    rows = NamedTuple[]
    checkerboard = side ÷ 2
    for m in -5:5
        momentum = mod(pi + 2pi * tau * m, 2pi)
        index = mod(round(Int, side * momentum / (2pi)), side)
        for (orientation, kx, ky) in (
            ("horizontal_satellite", index, checkerboard),
            ("vertical_satellite", checkerboard, index),
        )
            value = qtt_fourier_amplitude(transformed, kx, ky)
            push!(rows, (; m, orientation, kx, ky,
                momentum_x=2pi * kx / side,
                momentum_y=2pi * ky / side,
                real=real(value), imag=imag(value), power=abs2(value)))
        end
    end
    return rows
end

println("Extracting density diagonal QTT from $checkpoint_path")
flush(stdout)
density_extraction = @timed density_diagonal_mps(system)
density_source = density_extraction.value
window = hann_window_mps(system.sites, min(params.tci_tol, 1e-12))

results = NamedTuple[]
scale_rows = NamedTuple[]
sample_rows = NamedTuple[]
for policy in policies
    println("QTT charge analysis policy=$(policy.name) cutoff=$(policy.cutoff) maxdim=$(policy.maxdim)")
    flush(stdout)
    density = copy(density_source)
    truncate!(density; cutoff=policy.cutoff, maxdim=policy.maxdim)
    centered, mean_density, trace_value = centered_density_mps(
        density, system.sites; cutoff=policy.cutoff, maxdim=policy.maxdim,
    )

    real_timing = @timed qtt_multiscale_d2(
        centered; cutoff=policy.cutoff, maxdim=policy.maxdim,
        keep=:prefix, fit_scales=fit_scales,
    )
    raw_fft_timing = @timed qtt_fourier_square(
        centered, system.sites;
        sign=-1, cutoff_MPO=policy.cutoff,
        cutoff=policy.cutoff, maxdim=policy.maxdim,
    )
    raw_fft = raw_fft_timing.value
    raw_fourier_timing = @timed qtt_multiscale_d2(
        raw_fft; cutoff=policy.cutoff, maxdim=policy.maxdim,
        keep=:suffix, fit_scales=fit_scales,
    )

    windowed_timing = @timed MPO_MeanField.Quantics.automul(
        centered, window; cutoff=policy.cutoff, maxdim=policy.maxdim,
    )
    windowed = windowed_timing.value
    hann_fft_timing = @timed qtt_fourier_square(
        windowed, system.sites;
        sign=-1, cutoff_MPO=policy.cutoff,
        cutoff=policy.cutoff, maxdim=policy.maxdim,
    )
    hann_fft = hann_fft_timing.value
    hann_fourier_timing = @timed qtt_multiscale_d2(
        hann_fft; cutoff=policy.cutoff, maxdim=policy.maxdim,
        keep=:suffix, fit_scales=fit_scales,
    )

    real_result = real_timing.value
    raw_result = raw_fourier_timing.value
    hann_result = hann_fourier_timing.value
    real_norm = real(inner(centered, centered))
    raw_norm = real(inner(raw_fft, raw_fft))
    hann_real_norm = real(inner(windowed, windowed))
    hann_fourier_norm = real(inner(hann_fft, hann_fft))
    push!(results, (
        policy=policy.name, cutoff=policy.cutoff, maxdim=policy.maxdim,
        density_max_chi=maxlinkdim(density), centered_max_chi=maxlinkdim(centered),
        trace=real(trace_value), trace_error=abs(real(trace_value) - params.density * N),
        mean_density=real(mean_density),
        real_d2=real_result.d2, real_d2_r_squared=real_result.r_squared,
        real_ipr=real_result.finest_ipr, real_participation=real_result.participation,
        raw_fourier_d2=raw_result.d2, raw_fourier_d2_r_squared=raw_result.r_squared,
        raw_fourier_ipr=raw_result.finest_ipr,
        raw_fourier_participation=raw_result.participation,
        hann_fourier_d2=hann_result.d2, hann_fourier_d2_r_squared=hann_result.r_squared,
        hann_fourier_ipr=hann_result.finest_ipr,
        hann_fourier_participation=hann_result.participation,
        raw_fourier_max_chi=maxlinkdim(raw_fft),
        hann_windowed_max_chi=maxlinkdim(windowed), hann_fourier_max_chi=maxlinkdim(hann_fft),
        raw_parseval_relative_error=abs(raw_norm - real_norm) / real_norm,
        hann_parseval_relative_error=abs(hann_fourier_norm - hann_real_norm) / hann_real_norm,
        real_d2_time_s=real_timing.time, raw_fft_time_s=raw_fft_timing.time,
        raw_fourier_d2_time_s=raw_fourier_timing.time,
        hann_window_time_s=windowed_timing.time, hann_fft_time_s=hann_fft_timing.time,
        hann_fourier_d2_time_s=hann_fourier_timing.time,
    ))
    for (space, analysis) in (
        ("real", real_result),
        ("fourier_raw", raw_result),
        ("fourier_hann", hann_result),
    )
        for row in analysis.scales
            local_d2 = row.scale == 0 ? NaN : analysis.local_slopes[row.scale].d2
            push!(scale_rows, merge((policy=policy.name, space=space), row, (local_d2=local_d2,)))
        end
    end
    if policy.name == "tight"
        tau = sqrt(2.0) - 5.0 / 6.0
        append!(sample_rows, [merge((spectrum="raw",), row) for row in quasiperiodic_samples(raw_fft, tau)])
        append!(sample_rows, [merge((spectrum="hann",), row) for row in quasiperiodic_samples(hann_fft, tau)])
    end
end

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
try
    open(joinpath(temporary, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "created_at" => string(now()),
            "diagnostic" => "qtt_charge_fft_d2_from_sp2_checkpoint",
            "campaign" => campaign_path,
            "task_index" => task_index,
            "label" => label,
            "source_checkpoint" => checkpoint_path,
            "state_kind" => state_kind,
            "checkpoint_scf_converged" => loaded.converged,
            "checkpoint_scf_termination_reason" => string(loaded.termination_reason),
            "side" => side,
            "matrix_dimension" => N,
            "fit_scales" => fit_scales,
            "density_measure" => "|n(x,y)-mean(n)|^2 normalized to unit mass",
            "fourier_normalization" => "1/sqrt(N)",
            "momentum_interval" => "[0, 2pi)",
            "density_source_max_chi" => maxlinkdim(density_source),
            "density_extraction_time_s" => density_extraction.time,
            "hann_window_max_chi" => maxlinkdim(window),
        ))
    end
    open(joinpath(temporary, "configurations.csv"), "w") do io
        header = propertynames(first(results))
        write_row(io, header)
        for row in results
            write_row(io, Tuple(row))
        end
    end
    open(joinpath(temporary, "scales.csv"), "w") do io
        header = propertynames(first(scale_rows))
        write_row(io, header)
        for row in scale_rows
            write_row(io, Tuple(row))
        end
    end
    open(joinpath(temporary, "fourier_samples.csv"), "w") do io
        header = propertynames(first(sample_rows))
        write_row(io, header)
        for row in sample_rows
            write_row(io, Tuple(row))
        end
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end

println("SP2 checkpoint QTT charge analysis complete: $output")
for row in results
    println("$(row.policy): D2_real=$(row.real_d2) D2_k_raw=$(row.raw_fourier_d2) " *
            "D2_k_hann=$(row.hann_fourier_d2) parseval=$(row.raw_parseval_relative_error)")
end
