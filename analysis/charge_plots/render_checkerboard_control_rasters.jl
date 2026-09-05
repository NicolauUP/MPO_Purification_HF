#!/usr/bin/env julia

"""Render lightweight conference rasters for the V2=0, U=1 checkerboard control.

Usage:
  julia --project=analysis/charge_plots render_checkerboard_control_rasters.jl INPUT.h5 OUTPUT_DIR

The Fourier raster is centered on the checkerboard ordering vector by
multiplying n(x,y)-1/2 by (-1)^(x+y) before a Hann-windowed FFT. Thus its
central peak represents the original charge peak at (pi, pi).
"""

using FFTW
using HDF5

length(ARGS) == 2 || error("usage: render_checkerboard_control_rasters.jl INPUT.h5 OUTPUT_DIR")
input, output = abspath.(ARGS)
isfile(input) || error("input does not exist: $input")
ispath(output) && error("refusing to overwrite output: $output")

density = h5open(input, "r") do file
    Float64.(read(file, "density"))
end
side = size(density, 1)
size(density, 2) == side || error("density must be square")
raw = density .- 0.5
edge = max(10, side >>> 7)
interior = raw[(edge + 1):(side - edge), (edge + 1):(side - edge)]
scale = max(maximum(abs, interior), eps(Float64))

function diverging(value, scale)
    u = clamp(value / scale, -1.0, 1.0)
    if u < 0
        return (UInt8(round(255 * (1 + u))), UInt8(round(255 * (1 + u))), UInt8(255))
    end
    return (UInt8(255), UInt8(round(255 * (1 - u))), UInt8(round(255 * (1 - u))))
end

function spectrum_colour(value)
    u = clamp((value + 12) / 12, 0.0, 1.0)
    if u < 0.5
        v = u / 0.5
        return (UInt8(round(20 * v)), UInt8(round(30 + 160v)), UInt8(round(80 + 150v)))
    end
    v = (u - 0.5) / 0.5
    return (UInt8(round(20 + 235v)), UInt8(round(190 + 60v)), UInt8(round(230 * (1 - v))))
end

function write_ppm(path, image)
    height, width, _ = size(image)
    open(path, "w") do io
        write(io, "P6\n$width $height\n255\n")
        for y in 1:height, x in 1:width
            write(io, image[y, x, 1]); write(io, image[y, x, 2]); write(io, image[y, x, 3])
        end
    end
end

function density_raster(field, scale; magnification=1)
    n, m = size(field)
    image = Array{UInt8}(undef, m * magnification, n * magnification, 3)
    for y in 1:m, x in 1:n
        colour = diverging(field[x, m - y + 1], scale)
        ys = ((y - 1) * magnification + 1):(y * magnification)
        xs = ((x - 1) * magnification + 1):(x * magnification)
        for yy in ys, xx in xs
            image[yy, xx, 1] = colour[1]
            image[yy, xx, 2] = colour[2]
            image[yy, xx, 3] = colour[3]
        end
    end
    return image
end

carrier = [iseven(x + y) ? raw[x, y] : -raw[x, y] for x in 1:side, y in 1:side]
window = [0.5 * (1 - cos(2pi * (x - 1) / (side - 1))) *
          0.5 * (1 - cos(2pi * (y - 1) / (side - 1))) for x in 1:side, y in 1:side]
power = abs2.(circshift(fft(carrier .* window), (side >>> 1, side >>> 1)))
power ./= maximum(power)
log_power = log10.(max.(power, 1e-12))

fft_image = Array{UInt8}(undef, side, side, 3)
for y in 1:side, x in 1:side
    colour = spectrum_colour(log_power[x, side - y + 1])
    fft_image[y, x, 1] = colour[1]
    fft_image[y, x, 2] = colour[2]
    fft_image[y, x, 3] = colour[3]
end

crop_side = min(256, side - 2edge)
start = (side - crop_side) ÷ 2 + 1
crop = raw[start:(start + crop_side - 1), start:(start + crop_side - 1)]

tmp = output * ".tmp.$(getpid())"
mkpath(tmp)
try
    write_ppm(joinpath(tmp, "charge_density_full_2048.ppm"), density_raster(raw, scale))
    write_ppm(joinpath(tmp, "charge_density_center_256.ppm"), density_raster(crop, scale; magnification=4))
    write_ppm(joinpath(tmp, "charge_fft_checkerboard_centered.ppm"), fft_image)
    open(joinpath(tmp, "README.txt"), "w") do io
        println(io, "Charge-density control: V2=0, U=1, 2048x2048, checkerboard + seed.")
        println(io, "density files show n(x,y)-1/2 with white at zero, blue negative, red positive.")
        println(io, "FFT is Hann-windowed and carrier-shifted: its centre is the original (pi,pi) charge-order peak.")
        println(io, "edge margin excluded only when selecting the colour scale: $edge sites.")
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote control rasters: $output")
