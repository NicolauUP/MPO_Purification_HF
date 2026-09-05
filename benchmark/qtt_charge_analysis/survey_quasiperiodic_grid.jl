#!/usr/bin/env julia

"""Evaluate QTT Fourier power on a quasiperiodic momentum grid.

Usage: survey_quasiperiodic_grid.jl CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY MMIN MMAX

For m,n in [MMIN,MMAX], evaluates the nearest available Fourier bins to
k_m = pi + 2*pi*tau*m and writes P_mn = |n_tilde(k_m,k_n)|^2.
The full density raster is never materialized.
"""

using Dates
using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using TOML

length(ARGS) == 6 || error("usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY MMIN MMAX")
campaign_path = abspath(ARGS[1])
task_text = ARGS[2]
result_text = abspath(ARGS[3])
output = abspath(ARGS[4])
task = parse(Int, task_text)
mmin, mmax = parse.(Int, ARGS[5:6])
mmin <= mmax || error("MMIN must not exceed MMAX")
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
ispath(output) && error("refusing to overwrite existing output: $output")

result_directory = result_text
checkpoint_path = let candidates = (joinpath(result_directory, "converged_state.h5"), joinpath(result_directory, "final_fixed_sp2_state.h5"))
    found = findfirst(isfile, candidates)
    isnothing(found) && error("no converged checkpoint in $result_directory")
    candidates[found]
end

namespace = Module(gensym(:QTTGridCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
campaign = Base.invokelatest(getfield, namespace, :campaign)
public_case = campaign_case(campaign, task)
params0 = Base.invokelatest(legacy_parameters, public_case.model, public_case.representation, public_case.solver)
checkpoint_maxdim = h5open(checkpoint_path, "r") do file Int(read(file, "itensors_maxdim")) end
params = ParametersSquare(L=params0.L, t=params0.t, U=params0.U, W=params0.W, S=params0.S,
    tci_tol=params0.tci_tol, itensors_tol=params0.itensors_tol, itensors_maxdim=checkpoint_maxdim,
    density=params0.density, purification_steps=params0.purification_steps, scf_mixing=params0.scf_mixing,
    scf_tol=params0.scf_tol, scf_max_iterations=params0.scf_max_iterations)
system = read_mpo_checkpoint(checkpoint_path, params).system
levels = length(system.sites)
levels % 2 == 0 || error("square QTT requires an even number of levels")
side = 1 << (levels ÷ 2)
tau = sqrt(2.0) - 5.0 / 6.0

println("Building QTT Fourier transform for $(side)x$(side), m,n=$mmin:$mmax")
density = density_diagonal_mps(system)
centered, mean_density, trace_value = centered_density_mps(density, system.sites;
    cutoff=params.itensors_tol, maxdim=checkpoint_maxdim)
function checkerboard(linear)
    x, y = square_lattice_decoder(Int(linear), levels)
    iseven(x + y) ? 1.0 : -1.0
end
_, _, checkerboard_state = MPO_MeanField.Quantics_TCI(checkerboard, Float64, system.sites, 1e-12)
demodulated = MPO_MeanField.Quantics.automul(centered, checkerboard_state;
    cutoff=params.itensors_tol, maxdim=checkpoint_maxdim)
_, _, ones_state = MPO_MeanField.Quantics_TCI(_ -> 1.0, Float64, system.sites, 1e-12)
carrier_mean = real(inner(ones_state, demodulated)) / side^2
demodulated = +(demodulated, -carrier_mean * ones_state;
    cutoff=params.itensors_tol, maxdim=checkpoint_maxdim)
fft_timing = @timed qtt_fourier_square(demodulated, system.sites; sign=-1,
    cutoff_MPO=params.itensors_tol, cutoff=params.itensors_tol, maxdim=checkpoint_maxdim)
fourier = fft_timing.value
wrap(k) = mod(k, side)
bin(k) = wrap(round(Int, side * mod(k, 2pi) / (2pi)))
indices = Dict(m => bin(pi + 2pi * tau * m) for m in mmin:mmax)
physical = Dict(m => mod(pi + 2pi * tau * m, 2pi) for m in mmin:mmax)
power = Matrix{Float64}(undef, mmax - mmin + 1, mmax - mmin + 1)
for (ix, m) in enumerate(mmin:mmax), (iy, n) in enumerate(mmin:mmax)
    power[ix, iy] = abs2(qtt_fourier_amplitude(fourier, indices[m], indices[n]))
end
normalization = maximum(power)
mkpath(output)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict("created_at" => string(now()), "diagnostic" => "qtt_quasiperiodic_fourier_power_grid",
        "source_checkpoint" => checkpoint_path, "side" => side, "matrix_dimension" => side^2,
        "tau" => tau, "m_min" => mmin, "m_max" => mmax, "fourier_normalization" => "1/sqrt(N)",
        "mean_density" => real(mean_density), "trace" => real(trace_value),
        "fourier_max_chi" => maxlinkdim(fourier), "fourier_time_s" => fft_timing.time,
        "momentum_definition" => "demodulated field: k_m = mod(2*pi*tau*m, 2*pi); raw density carrier removed first",
        "checkerboard_carrier_mean_removed" => carrier_mean,
        "bin_definition" => "nearest integer Fourier bin"))
end
open(joinpath(output, "momentum_bins.csv"), "w") do io
    println(io, "m,k_physical,k_bin")
    for m in mmin:mmax println(io, m, ",", physical[m], ",", indices[m]) end
end
open(joinpath(output, "power.csv"), "w") do io
    println(io, "m,n,power,relative_power")
    for (ix, m) in enumerate(mmin:mmax), (iy, n) in enumerate(mmin:mmax)
        println(io, m, ",", n, ",", power[ix, iy], ",", power[ix, iy] / normalization)
    end
end
println("QTT quasiperiodic Fourier grid complete: $output")
