#!/usr/bin/env julia

"""Render a square charge-density CSV and its Hann-windowed demodulated DFT.

Usage:
  julia render_square_site_density_fft.jl RESULT_DIRECTORY OUTPUT_DIRECTORY

This postprocessor deliberately depends only on Julia standard libraries.  It
works for both dense-HF and MPO-HF outputs because both write `site_density.csv`.
The site ordering is the public square/QTT ordering used by `site_density.csv`.
"""

using TOML
using Statistics

length(ARGS) == 2 || error("usage: $(PROGRAM_FILE) RESULT_DIRECTORY OUTPUT_DIRECTORY")
result = abspath(ARGS[1])
output = abspath(ARGS[2])
density_path = joinpath(result, "site_density.csv")
metadata_path = joinpath(result, "metadata.toml")
isfile(density_path) || error("missing $density_path")
isfile(metadata_path) || error("missing $metadata_path")
ispath(output) && error("refusing to overwrite existing output: $output")

metadata = TOML.parsefile(metadata_path)
size_xy = metadata["physical_size"]
length(size_xy) == 2 && size_xy[1] == size_xy[2] || error("expected a square physical_size")
side = Int(size_xy[1])
side > 1 || error("side must exceed one")

# Interleaved QTT order: bits of x and y alternate, most-significant first.
function decode_square_site(linear::Int, side::Int)
    bits = trailing_zeros(side)
    x = 0
    y = 0
    for bit in 0:(bits - 1)
        # Keep exactly the convention of `square_lattice_decoder`: even
        # interleaved bits encode y and odd bits encode x.
        y |= ((linear >> (2bit)) & 1) << bit
        x |= ((linear >> (2bit + 1)) & 1) << bit
    end
    return x, y
end

density = zeros(Float64, side, side)
lines = readlines(density_path)
length(lines) == side^2 + 1 || error("unexpected density row count")
for line in @view lines[2:end]
    site, value = split(line, ',')
    x, y = decode_square_site(parse(Int, replace(site, '"' => "")) - 1, side)
    density[x + 1, y + 1] = parse(Float64, replace(value, '"' => ""))
end

raw = density .- 0.5
envelope = [(-1.0)^(x + y) * raw[x, y] for x in 1:side, y in 1:side]
envelope .-= mean(envelope)
window = [0.5 * (1 - cos(2pi * (x - 1) / (side - 1))) *
          0.5 * (1 - cos(2pi * (y - 1) / (side - 1)))
          for x in 1:side, y in 1:side]
signal = envelope .* window
fourier = [cis(-2pi * k * x / side) for k in 0:(side - 1), x in 0:(side - 1)]
amplitudes = fourier * ComplexF64.(signal) * transpose(fourier)
power = zeros(Float64, side, side)
for ky in 0:(side - 1), kx in 0:(side - 1)
    power[mod(kx + side ÷ 2, side) + 1, mod(ky + side ÷ 2, side) + 1] = abs2(amplitudes[kx + 1, ky + 1])
end
power ./= maximum(power)

function diverging(value, scale)
    u = clamp(value / scale, -1.0, 1.0)
    u < 0 && return (round(UInt8, 255 * (1 + u)), round(UInt8, 255 * (1 + u)), UInt8(255))
    return (UInt8(255), round(UInt8, 255 * (1 - u)), round(UInt8, 255 * (1 - u)))
end
function spectrum_colour(value)
    u = clamp((value + 8) / 8, 0.0, 1.0)
    u < 0.55 && return (round(UInt8, 255 * (1 - u / 0.55)), round(UInt8, 255 * (1 - 0.35u / 0.55)), UInt8(255))
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
function field_image(field, colour; cell=16)
    field_side = size(field, 1)
    size(field, 2) == field_side || error("field image requires a square matrix")
    image = fill(UInt8(255), field_side * cell, field_side * cell, 3)
    scale = maximum(abs, field)
    for y in 1:field_side, x in 1:field_side
        rgb = colour(field[x, field_side - y + 1], scale)
        image[((y - 1) * cell + 1):(y * cell), ((x - 1) * cell + 1):(x * cell), :] .=
            reshape(UInt8[rgb...], 1, 1, 3)
    end
    return image
end

tmp = output * ".tmp.$(getpid())"
mkpath(tmp)
try
    write_ppm(joinpath(tmp, "charge_density_demodulated.ppm"), field_image(envelope, diverging))
    # The bulk image is the presentation diagnostic: exclude open-boundary
    # response before setting its colour scale, while leaving the full field
    # and its Hann-windowed FFT available for audit.
    crop_margin = max(4, side >>> 3)
    bulk = envelope[(crop_margin + 1):(side - crop_margin),
                    (crop_margin + 1):(side - crop_margin)]
    write_ppm(joinpath(tmp, "charge_density_demodulated_bulk.ppm"), field_image(bulk, diverging))
    write_ppm(
        joinpath(tmp, "charge_density_demodulated_hann_fft.ppm"),
        field_image(log10.(max.(power, 1e-8)), (value, _) -> spectrum_colour(value)),
    )
    open(joinpath(tmp, "strongest_fft_peaks.csv"), "w") do io
        println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
        peaks = Tuple{Float64,Int,Int}[]
        for y in 1:side, x in 1:side
            abs(x - (side ÷ 2 + 1)) <= 1 && abs(y - (side ÷ 2 + 1)) <= 1 && continue
            push!(peaks, (power[x, y], x - (side ÷ 2 + 1), y - (side ÷ 2 + 1)))
        end
        sort!(peaks; by=first, rev=true)
        for (rank, (weight, kx, ky)) in enumerate(peaks[1:min(32, length(peaks))])
            println(io, "$rank,$kx,$ky,$(2pi * kx / side),$(2pi * ky / side),$weight")
        end
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote $output")
