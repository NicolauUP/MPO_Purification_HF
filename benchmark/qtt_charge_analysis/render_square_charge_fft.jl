#!/usr/bin/env julia

"""Render real-space charge-density and Fourier maps from a square MPO checkpoint.

Usage:
    julia --project=. render_square_charge_fft.jl CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY

The script is intentionally a lightweight figure postprocessor. It never
restarts SCF, rebuilds mean fields, or measures the full energy. For the
small reference lattice used here, it contracts the QTT density MPS once and
uses a direct, Hann-windowed DFT. This keeps the visual diagnostic independent
of an optional FFT package and makes its normalization explicit.
"""

using Dates
using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using LinearAlgebra
using Statistics
using TOML

length(ARGS) == 4 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY",
)
campaign_path = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
result_directory = abspath(ARGS[3])
output = abspath(ARGS[4])
checkpoint_path = joinpath(result_directory, "converged_state.h5")
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
isfile(checkpoint_path) || error("checkpoint does not exist: $checkpoint_path")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:SquareChargeFigureCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
campaign = getfield(namespace, :campaign)
case = campaign_case(campaign, task_index)
params = legacy_parameters(case.model, case.representation, case.solver)
checkpoint_cap = h5open(checkpoint_path, "r") do file
    Int(read(file, "itensors_maxdim"))
end
checkpoint_cap == params.itensors_maxdim || error(
    "checkpoint maxdim=$checkpoint_cap does not match campaign maxdim=$(params.itensors_maxdim)",
)
loaded = read_mpo_checkpoint(checkpoint_path, params)
system = loaded.system
levels = length(system.sites)
iszero(levels % 2) || error("square QTT requires an even number of levels")
side = 1 << (levels >>> 1)

println("Extracting $(side)x$(side) density from $checkpoint_path")
density_mps = density_diagonal_mps(system)
qtt_dense = vec(Array(reduce(*, density_mps), system.sites...))
density = Matrix{Float64}(undef, side, side)
for linear in 0:(side^2 - 1)
    x, y = square_lattice_decoder(linear, levels)
    reversed = Int(bitreverse(UInt(linear)) >>> (8sizeof(UInt) - levels))
    density[x + 1, y + 1] = real(qtt_dense[reversed + 1])
end

# Separate the checkerboard carrier from the slowly varying quasiperiodic
# envelope before Fourier analysis. The separable Hann window suppresses the
# open-boundary edge background; it broadens, but does not shift, peaks.
raw = density .- 0.5
envelope = [(-1.0)^(x + y) * raw[x, y] for x in 1:side, y in 1:side]
envelope .-= mean(envelope)
window = [0.5 * (1 - cos(2pi * (x - 1) / (side - 1))) *
          0.5 * (1 - cos(2pi * (y - 1) / (side - 1)))
          for x in 1:side, y in 1:side]
signal = envelope .* window

fourier_matrix = [cis(-2pi * k * x / side) for k in 0:(side - 1), x in 0:(side - 1)]
amplitudes = fourier_matrix * ComplexF64.(signal) * transpose(fourier_matrix)
power = zeros(Float64, side, side)
for ky in 0:(side - 1), kx in 0:(side - 1)
    power[mod(kx + side ÷ 2, side) + 1, mod(ky + side ÷ 2, side) + 1] =
        abs2(amplitudes[kx + 1, ky + 1])
end
power ./= maximum(power)

function diverging(value, scale)
    u = clamp(value / scale, -1.0, 1.0)
    u < 0 && return (round(UInt8, 255 * (1 + u)), round(UInt8, 255 * (1 + u)), UInt8(255))
    return (UInt8(255), round(UInt8, 255 * (1 - u)), round(UInt8, 255 * (1 - u)))
end

function spectrum_colour(value)
    u = clamp((value + 8) / 8, 0.0, 1.0)
    u < 0.55 && return (round(UInt8, 255 * (1 - u / 0.55)),
                        round(UInt8, 255 * (1 - 0.35u / 0.55)), UInt8(255))
    v = (u - 0.55) / 0.45
    return (UInt8(0), round(UInt8, 166 * (1 - v)), round(UInt8, 255 * (1 - v)))
end

function write_ppm(path, image)
    height, width, _ = size(image)
    open(path, "w") do io
        write(io, "P6\n$width $height\n255\n")
        write(io, reinterpret(UInt8, vec(permutedims(image, (3, 2, 1)))))
    end
end

cell = 16
gap = 48
panel = side * cell
real_image = fill(UInt8(255), panel, 2panel + gap, 3)
for (panel_index, (field, scale)) in enumerate(((raw, maximum(abs, raw)), (envelope, maximum(abs, envelope))))
    xoffset = (panel_index - 1) * (panel + gap)
    for y in 1:side, x in 1:side
        colour = diverging(field[x, side - y + 1], scale)
        xs = (xoffset + (x - 1) * cell + 1):(xoffset + x * cell)
        ys = ((y - 1) * cell + 1):(y * cell)
        for yy in ys, xx in xs
            real_image[yy, xx, :] .= UInt8[colour...]
        end
    end
end

# A separate bulk-only panel avoids interpreting the open-boundary response as
# quasiperiodic structure. The FFT above still uses the full Hann-windowed
# field; this crop is only for real-space presentation.
crop_margin = max(8, side ÷ 16)
crop_side = side - 2crop_margin
bulk_image = fill(UInt8(255), crop_side * cell, crop_side * cell, 3)
bulk_scale = maximum(abs, envelope[(crop_margin + 1):(side - crop_margin),
                                   (crop_margin + 1):(side - crop_margin)])
for y in 1:crop_side, x in 1:crop_side
    colour = diverging(envelope[x + crop_margin, side - (y + crop_margin) + 1], bulk_scale)
    xs = ((x - 1) * cell + 1):(x * cell)
    ys = ((y - 1) * cell + 1):(y * cell)
    for yy in ys, xx in xs
        bulk_image[yy, xx, :] .= UInt8[colour...]
    end
end

fft_image = fill(UInt8(255), panel, panel, 3)
log_power = log10.(max.(power, 1e-8))
for y in 1:side, x in 1:side
    colour = spectrum_colour(log_power[x, side - y + 1])
    xs = ((x - 1) * cell + 1):(x * cell)
    ys = ((y - 1) * cell + 1):(y * cell)
    for yy in ys, xx in xs
        fft_image[yy, xx, :] .= UInt8[colour...]
    end
end

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
try
    write_ppm(joinpath(temporary, "charge_density_raw_and_demodulated.ppm"), real_image)
    write_ppm(joinpath(temporary, "charge_density_demodulated_bulk.ppm"), bulk_image)
    write_ppm(joinpath(temporary, "charge_density_demodulated_hann_fft.ppm"), fft_image)
    open(joinpath(temporary, "density.csv"), "w") do io
        println(io, "x,y,density,raw,checkerboard_demodulated")
        for x in 1:side, y in 1:side
            println(io, "$(x - 1),$(y - 1),$(density[x, y]),$(raw[x, y]),$(envelope[x, y])")
        end
    end
    open(joinpath(temporary, "fft_power.csv"), "w") do io
        println(io, "kx_index,ky_index,kx,ky,relative_power")
        for x in 1:side, y in 1:side
            kx = x - (side ÷ 2 + 1)
            ky = y - (side ÷ 2 + 1)
            println(io, "$kx,$ky,$(2pi * kx / side),$(2pi * ky / side),$(power[x, y])")
        end
    end
    peaks = Tuple{Float64,Int,Int}[]
    for y in 1:side, x in 1:side
        abs(x - (side ÷ 2 + 1)) <= 1 && abs(y - (side ÷ 2 + 1)) <= 1 && continue
        push!(peaks, (power[x, y], x - (side ÷ 2 + 1), y - (side ÷ 2 + 1)))
    end
    sort!(peaks; by=first, rev=true)
    open(joinpath(temporary, "strongest_fft_peaks.csv"), "w") do io
        println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
        for (rank, (weight, kx, ky)) in enumerate(peaks[1:min(32, length(peaks))])
            println(io, "$rank,$kx,$ky,$(2pi * kx / side),$(2pi * ky / side),$weight")
        end
    end
    open(joinpath(temporary, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "created_at" => string(now()), "campaign" => campaign_path,
            "task_index" => task_index, "source_checkpoint" => checkpoint_path,
            "side" => side, "matrix_dimension" => side^2,
            "density_qtt_max_chi" => maxlinkdim(density_mps),
            "real_space_definition" => "raw=n(x,y)-1/2; envelope=(-1)^(x+y) raw - mean",
            "fft_definition" => "direct DFT of Hann-windowed checkerboard-demodulated envelope",
            "fft_power_normalization" => "power/max(power)",
            "demodulated_bulk_crop_margin" => crop_margin,
            "demodulated_bulk_crop_side" => crop_side,
            "demodulated_bulk_scale" => bulk_scale,
        ))
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end
println("Wrote charge and FFT figures to $output")
