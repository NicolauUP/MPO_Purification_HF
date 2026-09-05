#!/usr/bin/env julia

"""Sampled audit of a converged square MPO checkpoint.

Usage:
  julia audit_mpo_checkpoint_sampled.jl CAMPAIGN.jl TASK CHECKPOINT OUTPUT [STRIDE] [CROP_SIDE] [CROP_STRIDE]

The audit intentionally avoids `observables_square`, which performs a full
site-and-bond sweep and evaluates the energy. It computes global MPO invariants
once, then samples a regular coarse grid and a centred full-resolution crop.
"""

using LinearAlgebra
using Statistics
using TOML
using HDF5
using ITensors
using ITensorMPS
using MPO_MeanField

4 <= length(ARGS) <= 7 || error("usage: CAMPAIGN.jl TASK CHECKPOINT OUTPUT [STRIDE] [CROP_SIDE] [CROP_STRIDE]")
campaign_path = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
checkpoint = abspath(ARGS[3])
output = abspath(ARGS[4])
stride = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 16
crop_side = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 256
crop_stride = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 1
stride > 0 || error("stride must be positive")
crop_stride > 0 || error("crop_stride must be positive")
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
isfile(checkpoint) || error("checkpoint does not exist: $checkpoint")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:CheckpointAuditCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
campaign = getfield(namespace, :campaign)
case = campaign_case(campaign, task_index)
params = legacy_parameters(case.model, case.representation, case.solver)
loaded = read_mpo_checkpoint(checkpoint, params)
system = loaded.system
side = 2^div(params.L, 2)
iszero(side % stride) || error("stride=$stride must divide side=$side")
0 < crop_side <= side || error("crop_side must lie in 1:$side")
iszero(crop_side % crop_stride) || error("crop_side=$crop_side must be divisible by crop_stride=$crop_stride")

println("Sampled checkpoint audit: side=$side stride=$stride crop_side=$crop_side")
rho = system.ρ
rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
trace_value = real(tr(rho))
target_particles = round(Int, side^2 * params.density)
trace_error = abs(trace_value - target_particles)
idempotency = MPO_MeanField.idempotency_residual(rho, rho_squared)
hermiticity = MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
println("global trace=$trace_value trace_error=$trace_error idempotency=$idempotency hermiticity=$hermiticity")

function sampled_density(system, x_values, y_values, side)
    result = Matrix{Float64}(undef, length(x_values), length(y_values))
    for (ix, x) in enumerate(x_values), (iy, y) in enumerate(y_values)
        site = square_lattice_index(x, y, system.params.L)
        result[ix, iy] = real(MatrixChecker(
            system.ρ, system.sites, site, site, system.bra_states, system.ket_states,
        ))
    end
    return result
end

coarse_values = 0:stride:(side - 1)
crop_start = (side - crop_side) ÷ 2
crop_values = crop_start:crop_stride:(crop_start + crop_side - 1)
coarse = sampled_density(system, coarse_values, coarse_values, side)
crop = sampled_density(system, crop_values, crop_values, side)

function write_field(path, field, xs, ys)
    open(path, "w") do io
        println(io, "x,y,density,raw,checkerboard_demodulated")
        for (ix, x) in enumerate(xs), (iy, y) in enumerate(ys)
            value = field[ix, iy]
            raw = value - 0.5
            demodulated = (-1.0)^(x + y) * raw
            println(io, "$x,$y,$value,$raw,$demodulated")
        end
    end
end

function fft_peaks(field, xs, ys)
    nx, ny = size(field)
    physical_step = xs[2] - xs[1]
    demodulated = [(-1.0)^(xs[ix] + ys[iy]) * (field[ix, iy] - 0.5)
                   for ix in 1:nx, iy in 1:ny]
    demodulated .-= mean(demodulated)
    window = [0.5 * (1 - cos(2pi * (ix - 1) / (nx - 1))) *
              0.5 * (1 - cos(2pi * (iy - 1) / (ny - 1)))
              for ix in 1:nx, iy in 1:ny]
    signal = demodulated .* window
    fourier = [cis(-2pi * k * x / nx) for k in 0:(nx - 1), x in 0:(nx - 1)]
    amplitudes = fourier * ComplexF64.(signal) * transpose(fourier)
    power = zeros(Float64, nx, ny)
    for ky in 0:(ny - 1), kx in 0:(nx - 1)
        power[mod(kx + nx ÷ 2, nx) + 1, mod(ky + ny ÷ 2, ny) + 1] = abs2(amplitudes[kx + 1, ky + 1])
    end
    power ./= maximum(power)
    peaks = Tuple{Float64, Int, Int}[]
    for iy in 1:ny, ix in 1:nx
        abs(ix - (nx ÷ 2 + 1)) <= 1 && abs(iy - (ny ÷ 2 + 1)) <= 1 && continue
        push!(peaks, (power[ix, iy], ix - (nx ÷ 2 + 1), iy - (ny ÷ 2 + 1)))
    end
    sort!(peaks; by=first, rev=true)
    return peaks, physical_step
end

tmp = output * ".tmp.$(getpid())"
mkpath(tmp)
try
    write_field(joinpath(tmp, "density_coarse.csv"), coarse, coarse_values, coarse_values)
    write_field(joinpath(tmp, "density_center_crop.csv"), crop, crop_values, crop_values)
    peaks, physical_step = fft_peaks(coarse, coarse_values, coarse_values)
    open(joinpath(tmp, "coarse_fft_peaks.csv"), "w") do io
        println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
        for (rank, (weight, kx, ky)) in enumerate(peaks[1:min(32, length(peaks))])
            # The coarse sample has physical spacing `stride`, so report
            # wavevectors in the original lattice units (and expose aliasing
            # honestly rather than labelling coarse-grid frequencies as k).
            println(io, "$rank,$kx,$ky,$(2pi * kx / side),$(2pi * ky / side),$weight")
        end
    end
    open(joinpath(tmp, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "campaign" => campaign_path, "task_index" => task_index,
            "source_checkpoint" => checkpoint, "side" => side,
            "matrix_dimension" => side^2, "sample_stride" => stride,
            "coarse_side" => length(coarse_values), "crop_side" => crop_side,
            "crop_stride" => crop_stride,
            "crop_start" => crop_start, "target_particles" => target_particles,
            "trace" => trace_value,
            "trace_error" => trace_error, "idempotency_residual" => idempotency,
            "hermiticity_residual" => hermiticity,
            "audit_scope" => "global invariants plus coarse grid and centred crop; no energy or full bond sweep",
            "fft_scope" => "Hann-windowed checkerboard-demodulated coarse grid",
        ))
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote sampled audit to $output")
