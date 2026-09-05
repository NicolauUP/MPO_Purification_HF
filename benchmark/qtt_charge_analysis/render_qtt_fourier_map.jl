#!/usr/bin/env julia

"""Render a Hann-windowed QTT density structure-factor map from an MPO checkpoint.

The transform is QTT-native. Only its final rank-compressed MPS is contracted
once to a dense raster for visualisation; no individual-site oracle loop or
coarse real-space averaging is used.
"""

using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using Statistics

length(ARGS) in (4, 5, 6, 7) || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY [ROTATION_RADIANS] [markers|none] [raw|demodulate]",
)
campaign_path, result_directory, output = abspath.(ARGS[[1, 3, 4]])
task = parse(Int, ARGS[2])
rotation = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.0
marker_mode = length(ARGS) >= 6 ? ARGS[6] : "markers"
field_mode = length(ARGS) >= 7 ? ARGS[7] : "raw"
marker_mode in ("markers", "none") || error("marker mode must be markers or none")
field_mode in ("raw", "demodulate") || error("field mode must be raw or demodulate")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:QTTFourierCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
case = campaign_case(getfield(namespace, :campaign), task)
params = legacy_parameters(case.model, case.representation, case.solver)
checkpoint = joinpath(result_directory, "converged_state.h5")
cap = h5open(checkpoint, "r") do file
    Int(read(file, "itensors_maxdim"))
end
params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=params.itensors_tol, itensors_maxdim=cap,
    density=params.density, purification_steps=params.purification_steps,
    scf_mixing=params.scf_mixing, scf_tol=params.scf_tol, scf_max_iterations=params.scf_max_iterations,
)
system = read_mpo_checkpoint(checkpoint, params).system
levels = length(system.sites)
side = 1 << (levels ÷ 2)
iseven(levels) || error("QTT Fourier map requires an even number of QTT levels")

function hann_window_mps(sites, tolerance)
    function window(linear)
        x, y = square_lattice_decoder(Int(linear), levels)
        return (0.5 * (1 - cos(2pi * x / (side - 1)))) *
               (0.5 * (1 - cos(2pi * y / (side - 1))))
    end
    _, _, state = MPO_MeanField.Quantics_TCI(window, Float64, sites, tolerance)
    return state
end

function checkerboard_mps(sites, tolerance)
    function checkerboard(linear)
        x, y = square_lattice_decoder(Int(linear), levels)
        return iseven(x + y) ? 1.0 : -1.0
    end
    _, _, state = MPO_MeanField.Quantics_TCI(checkerboard, Float64, sites, tolerance)
    return state
end

println("Extracting, windowing, and QTT-transforming the saved charge density...")
flush(stdout)
density = density_diagonal_mps(system)
centered, _, _ = centered_density_mps(density, system.sites; cutoff=1e-10, maxdim=cap)
analysis_field = field_mode == "demodulate" ? MPO_MeanField.Quantics.automul(
    centered, checkerboard_mps(system.sites, 1e-12); cutoff=1e-10, maxdim=cap,
) : centered
windowed = MPO_MeanField.Quantics.automul(
    analysis_field, hann_window_mps(system.sites, 1e-12); cutoff=1e-10, maxdim=cap,
)
transformed = qtt_fourier_square(
    windowed, system.sites; cutoff_MPO=1e-10, cutoff=1e-10, maxdim=cap,
)
println("QTT ranks: density=$(maxlinkdim(density)), field=$(maxlinkdim(analysis_field)), Fourier=$(maxlinkdim(transformed)); rendering raster...")
flush(stdout)

# Contract once. Quantics stores the first bit as the most-significant bit,
# while Julia's dense tensor layout uses it as the fastest array index.
raw = vec(Array(reduce(*, transformed), system.sites...))
power = Matrix{Float64}(undef, side, side)
for linear in 0:(side^2 - 1)
    x, y = square_lattice_decoder(linear, levels)
    reversed = Int(bitreverse(UInt(linear)) >> (8sizeof(UInt) - levels))
    power[x + 1, y + 1] = abs2(raw[reversed + 1])
end

function write_ppm(path, rgb)
    open(path, "w") do io
        write(io, "P6\n$(size(rgb, 2)) $(size(rgb, 1))\n255\n")
        write(io, reinterpret(UInt8, vec(permutedims(rgb, (3, 2, 1)))))
    end
end
function magma(value)
    value = clamp(value, 0.0, 1.0)
    value < 0.35 && return (round(UInt8, 8 + 67value / 0.35), round(UInt8, 5 + 16value / 0.35), round(UInt8, 32 + 74value / 0.35))
    value < 0.7 && return (round(UInt8, 75 + 133(value - 0.35) / 0.35), round(UInt8, 21 + 43(value - 0.35) / 0.35), round(UInt8, 106 + 32(value - 0.35) / 0.35))
    value = (value - 0.7) / 0.3
    return (round(UInt8, 208 + 47value), round(UInt8, 64 + 172value), round(UInt8, 138 - 107value))
end

# Retain the physical convention k_x,k_y in [0,2pi). The displayed scale is
# deliberately fixed, so maps from future cap ladders can be compared directly.
logpower = log10.(power .+ eps(Float64))
log_floor, log_ceiling = -13.0, -1.0
rgb = Array{UInt8}(undef, side, side, 3)
for y in 1:side, x in 1:side
    c = magma((logpower[x, y] - log_floor) / (log_ceiling - log_floor))
    rgb[y, x, 1] = c[1]; rgb[y, x, 2] = c[2]; rgb[y, x, 3] = c[3]
end

# Crosses mark the reciprocal axes of the hopping modulation, displaced from
# the checkerboard wave vector. theta=0 is the existing separable convention.
function mark_cross!(image, px, py)
    for offset in -5:5, (x, y) in ((px + offset, py), (px, py + offset))
        1 <= x <= side && 1 <= y <= side && (image[y, x, :] .= UInt8[255, 255, 255])
    end
end
centre = side ÷ 2 + 1
tau = 0.5808802290397618
rotation_cos, rotation_sin = cos(rotation), sin(rotation)
momentum_pixel(k) = mod(round(Int, mod(k, 2pi) / (2pi) * side), side) + 1
if marker_mode == "markers"
    for m in -5:5
        displacement = 2pi * tau * m
        mark_cross!(rgb,
            momentum_pixel(pi + displacement * rotation_cos),
            momentum_pixel(pi + displacement * rotation_sin),
        )
        mark_cross!(rgb,
            momentum_pixel(pi - displacement * rotation_sin),
            momentum_pixel(pi + displacement * rotation_cos),
        )
    end
end

mkpath(output)
write_ppm(joinpath(output, "density_qtt_structure_factor.ppm"), rgb)
bar = Array{UInt8}(undef, side, 32, 3)
for y in 1:side
    c = magma((side - y) / (side - 1))
    bar[y, :, 1] .= c[1]; bar[y, :, 2] .= c[2]; bar[y, :, 3] .= c[3]
end
write_ppm(joinpath(output, "density_qtt_structure_factor_colorbar.ppm"), bar)

# Quantify the marked quasiperiodic points against a local leakage background.
# The 9x9 stencil excludes the central 3x3 pixels, which prevents a genuine
# narrow satellite (or its immediate Hann broadening) from defining its own
# reference background.
function local_background(array, x, y; radius=4, exclusion=1)
    values = Float64[]
    for dx in -radius:radius, dy in -radius:radius
        max(abs(dx), abs(dy)) <= exclusion && continue
        xx = mod1(x + dx, side); yy = mod1(y + dy, side)
        push!(values, array[xx, yy])
    end
    return median(values)
end
open(joinpath(output, "quasiperiodic_samples.csv"), "w") do io
    println(io, "m,orientation,kx,ky,power,local_background,peak_to_background")
    for m in -5:5
        displacement = 2pi * tau * m
        for (orientation, x, y) in (
            ("rotated_x", momentum_pixel(pi + displacement * rotation_cos),
                momentum_pixel(pi + displacement * rotation_sin)),
            ("rotated_y", momentum_pixel(pi - displacement * rotation_sin),
                momentum_pixel(pi + displacement * rotation_cos)),
        )
            background = local_background(power, x, y)
            println(io, join((m, orientation, 2pi * (x - 1) / side, 2pi * (y - 1) / side,
                power[x, y], background, power[x, y] / background), ','))
        end
    end
end
open(joinpath(output, "peaks.txt"), "w") do io
    for index in partialsortperm(vec(power), 1:20; rev=true)
        x, y = Tuple(CartesianIndices(power)[index])
        println(io, "kx=$(round(2pi * (x - 1) / side; digits=6)) ky=$(round(2pi * (y - 1) / side; digits=6)) power=$(power[x, y])")
    end
end
println("Wrote QTT Fourier map to $output")
