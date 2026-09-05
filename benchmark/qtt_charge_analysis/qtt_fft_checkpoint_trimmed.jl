#!/usr/bin/env julia

"""Whole-lattice QTT Fourier analysis of an MPO checkpoint.

Usage:
  julia qtt_fft_checkpoint_trimmed.jl CAMPAIGN.jl TASK CHECKPOINT OUTPUT [EDGE_MARGIN]

The density is represented as a QTT/MPS, multiplied by a checkerboard
demodulation and a hard edge mask, and transformed with the QTT Fourier
transform.  Only the final Fourier raster and presentation rasters are
materialized; no site-by-site MatrixChecker loop is used.
"""

using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using LinearAlgebra
using Statistics
using TOML

length(ARGS) in (4, 5) || error("usage: CAMPAIGN.jl TASK CHECKPOINT OUTPUT [EDGE_MARGIN]")
campaign_path = abspath(ARGS[1])
task = parse(Int, ARGS[2])
checkpoint = abspath(ARGS[3])
output = abspath(ARGS[4])
edge_margin = length(ARGS) == 5 ? parse(Int, ARGS[5]) : 10
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
isfile(checkpoint) || error("checkpoint does not exist: $checkpoint")
edge_margin >= 0 || error("edge margin must be nonnegative")
ispath(output) && error("refusing to overwrite existing output: $output")

namespace = Module(gensym(:TrimmedQTTFFT))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
case = campaign_case(getfield(namespace, :campaign), task)
params0 = legacy_parameters(case.model, case.representation, case.solver)
cap = h5open(checkpoint, "r") do file; Int(read(file, "itensors_maxdim")); end
params = ParametersSquare(
    L=params0.L, t=params0.t, U=params0.U, W=params0.W, S=params0.S,
    tci_tol=params0.tci_tol, itensors_tol=params0.itensors_tol,
    itensors_maxdim=cap, density=params0.density,
    purification_steps=params0.purification_steps, scf_mixing=params0.scf_mixing,
    scf_tol=params0.scf_tol, scf_max_iterations=params0.scf_max_iterations,
)
system = read_mpo_checkpoint(checkpoint, params).system
levels = length(system.sites)
iseven(levels) || error("square QTT requires an even number of levels")
side = 1 << (levels ÷ 2)
edge_margin < side ÷ 2 || error("edge margin must leave a nonempty interior")

function mask_mps(sites)
    function mask(linear)
        x, y = square_lattice_decoder(Int(linear), levels)
        return (edge_margin <= x < side - edge_margin && edge_margin <= y < side - edge_margin) ? 1.0 : 0.0
    end
    _, _, state = MPO_MeanField.Quantics_TCI(mask, Float64, sites, min(params.tci_tol, 1e-12))
    return state
end
function carrier_mps(sites)
    function carrier(linear)
        x, y = square_lattice_decoder(Int(linear), levels)
        return iseven(x + y) ? 1.0 : -1.0
    end
    _, _, state = MPO_MeanField.Quantics_TCI(carrier, Float64, sites, min(params.tci_tol, 1e-12))
    return state
end
function hann_mps(sites)
    function window(linear)
        x, y = square_lattice_decoder(Int(linear), levels)
        0.5 * (1 - cos(2pi * x / (side - 1))) *
        0.5 * (1 - cos(2pi * y / (side - 1)))
    end
    _, _, state = MPO_MeanField.Quantics_TCI(window, Float64, sites, min(params.tci_tol, 1e-12))
    return state
end

println("Whole-lattice QTT FFT: side=$side edge_margin=$edge_margin cap=$cap")
flush(stdout)
density = density_diagonal_mps(system)
centered, mean_density, trace_value = centered_density_mps(
    density, system.sites; cutoff=params.itensors_tol, maxdim=cap,
)
window = mask_mps(system.sites)
carrier = carrier_mps(system.sites)
hann = hann_mps(system.sites)
demodulated = MPO_MeanField.Quantics.automul(
    centered, carrier; cutoff=params.itensors_tol, maxdim=cap,
)
# Remove the residual carrier average.  Without this step the checkerboard
# demodulation leaves a large k=0 background that hides the quasiperiodic
# envelope in both the real-space colour scale and the FFT.
ones_state = MPO_MeanField._ones_mps(system.sites, Float64)
carrier_mean = real(inner(ones_state, demodulated)) / side^2
demodulated = +(
    demodulated, -carrier_mean * ones_state;
    cutoff=params.itensors_tol, maxdim=cap,
)
windowed = MPO_MeanField.Quantics.automul(
    demodulated, hann; cutoff=params.itensors_tol, maxdim=cap,
)
transformed = qtt_fourier_square(
    windowed, system.sites; cutoff_MPO=params.itensors_tol,
    cutoff=params.itensors_tol, maxdim=cap,
)
println("QTT ranks: density=$(maxlinkdim(density)) windowed=$(maxlinkdim(windowed)) Fourier=$(maxlinkdim(transformed))")
flush(stdout)

function dense_raster(state)
    raw = vec(Array(reduce(*, state), system.sites...))
    raster = Matrix{Float64}(undef, side, side)
    for linear in 0:(side^2 - 1)
        x, y = square_lattice_decoder(linear, levels)
        reversed = Int(bitreverse(UInt(linear)) >> (8sizeof(UInt) - levels))
        raster[x + 1, y + 1] = real(raw[reversed + 1])
    end
    return raster
end
fft_amplitudes = dense_raster(transformed)
power = abs2.(fft_amplitudes)
maximum(power) > 0 || error("Fourier transform has zero norm")
power ./= maximum(power)
power = circshift(power, (side ÷ 2, side ÷ 2))
real_field = dense_raster(demodulated)
mask = dense_raster(window)
trimmed_field = real_field .* mask

function density_rgb(v, scale)
    u = clamp(v / scale, -1.0, 1.0)
    if u < 0; q = UInt8(round(255 * (1 + u))); return (q, q, UInt8(255)); end
    q = UInt8(round(255 * (1 - u))); return (UInt8(255), q, q)
end
const fft_log_floor = -8.0
function spectrum_rgb(v)
    # The full dynamic range is enormous; clipping the displayed log-power
    # below -8 makes weak C8 satellites visible without changing the data.
    u = clamp((v - fft_log_floor) / (0.0 - fft_log_floor), 0.0, 1.0)
    if u < 0.5
        q = u / 0.5; return (UInt8(round(15q)), UInt8(round(55 + 180q)), UInt8(round(110 + 145q)))
    end
    q = (u - 0.5) / 0.5; return (UInt8(round(15 + 240q)), UInt8(round(235 - 35q)), UInt8(round(255 - 210q)))
end
function write_ppm(path, image)
    h, w, _ = size(image)
    open(path, "w") do io
        write(io, "P6\n$w $h\n255\n")
        for y in 1:h, x in 1:w; write(io, image[y, x, 1], image[y, x, 2], image[y, x, 3]); end
    end
end
function field_image(field, colour, scale)
    h, w = size(field); rgb = Array{UInt8}(undef, h, w, 3)
    for y in 1:h, x in 1:w
        c = colour(field[x, h - y + 1], scale)
        rgb[y, x, 1] = c[1]; rgb[y, x, 2] = c[2]; rgb[y, x, 3] = c[3]
    end
    return rgb
end

function box_blur(field, radius)
    h, w = size(field)
    integral = zeros(Float64, h + 1, w + 1)
    for y in 1:h, x in 1:w
        integral[y + 1, x + 1] = field[y, x] + integral[y, x + 1] +
            integral[y + 1, x] - integral[y, x]
    end
    blurred = Matrix{Float64}(undef, h, w)
    for y in 1:h, x in 1:w
        y0, y1 = max(1, y - radius), min(h, y + radius)
        x0, x1 = max(1, x - radius), min(w, x + radius)
        blurred[y, x] = (integral[y1 + 1, x1 + 1] - integral[y0, x1 + 1] -
            integral[y1 + 1, x0] + integral[y0, x0]) / ((y1 - y0 + 1) * (x1 - x0 + 1))
    end
    return blurred
end

function envelope_rgb(v, scale)
    u = clamp(v / scale, 0.0, 1.0)
    # Dark blue -> cyan -> yellow, suitable for a nonnegative envelope.
    if u < 0.5
        q = 2u
        return (UInt8(round(20 + 10q)), UInt8(round(35 + 130q)), UInt8(round(105 + 120q)))
    end
    q = 2(u - 0.5)
    return (UInt8(round(30 + 225q)), UInt8(round(165 + 70q)), UInt8(round(225 - 185q)))
end

tmp = output * ".tmp.$(getpid())"; mkpath(tmp)
try
    # Use a robust presentation scale: a single exceptional tensor entry
    # must not wash out the modulation over the rest of the lattice.
    scale = quantile(abs.(vec(trimmed_field)), 0.999)
    scale > 0 || error("real-space field has zero interior contrast")
    write_ppm(joinpath(tmp, "density_demodulated_edge10_contrast.ppm"), field_image(trimmed_field, density_rgb, scale))
    envelope = box_blur(abs.(trimmed_field), 8)
    envelope_scale = quantile(vec(envelope), 0.999)
    envelope_scale > 0 || error("real-space envelope has zero contrast")
    write_ppm(joinpath(tmp, "density_envelope_edge10_box8.ppm"), field_image(envelope, envelope_rgb, envelope_scale))
    logpower = clamp.(log10.(max.(power, 1e-12)), fft_log_floor, 0.0)
    write_ppm(joinpath(tmp, "whole_lattice_fft_hann_edge10_contrast.ppm"), field_image(logpower, (v, s) -> spectrum_rgb(v), 1.0))
    peaks = Tuple{Float64, Int, Int}[]
    centre = side ÷ 2 + 1
    for y in 1:side, x in 1:side
        abs(x - centre) <= 1 && abs(y - centre) <= 1 && continue
        push!(peaks, (power[x, y], x - centre, y - centre))
    end
    sort!(peaks; by=first, rev=true)
    open(joinpath(tmp, "whole_lattice_fft_peaks.csv"), "w") do io
        println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
        for (rank, (weight, kx, ky)) in enumerate(peaks[1:min(64, length(peaks))])
            println(io, "$rank,$kx,$ky,$(2pi*kx/side),$(2pi*ky/side),$weight")
        end
    end
    open(joinpath(tmp, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "campaign" => campaign_path, "task_index" => task,
            "source_checkpoint" => checkpoint, "side" => side,
            "matrix_dimension" => side^2, "edge_margin" => edge_margin,
            "trace" => real(trace_value), "mean_density" => real(mean_density),
            "demodulated_mean_removed" => carrier_mean,
            "density_max_chi" => maxlinkdim(density), "fourier_max_chi" => maxlinkdim(transformed),
            "fft_definition" => "whole-lattice QTT FFT of centered density times checkerboard carrier and Hann window",
            "edge_definition" => "10-site mask applied only to real-space presentation raster",
            "power_normalization" => "abs2(F)/maximum(abs2(F))",
            "fft_display_log_floor" => fft_log_floor,
            "real_display_scale" => scale,
            "envelope_definition" => "8-site-radius box average of abs(centered density times checkerboard carrier)",
            "envelope_display_scale" => envelope_scale,
        ))
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote whole-lattice QTT FFT to $output")
