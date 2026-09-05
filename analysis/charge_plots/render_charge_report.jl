#!/usr/bin/env julia

"""Render publication figures from `export_charge_hdf5.jl` output.

Usage:
    julia --project=analysis/charge_plots render_charge_report.jl INPUT.h5 OUTPUT_DIRECTORY [edge_margin] [crop_side]

The renderer is deliberately separate from the numerical solver environment.
It never reads an MPO and never reruns SCF.  Dense HDF5 exports produce full
heatmaps and a Hann-windowed FFT; sampled exports produce a sampled-density
figure and a clear notice that a full FFT requires a dense raster.
"""

if !isempty(ARGS) && ARGS[1] in ("-h", "--help")
    println("usage: $(PROGRAM_FILE) INPUT.h5 OUTPUT_DIRECTORY [edge_margin] [crop_side]")
    exit(0)
end
length(ARGS) >= 2 || error("usage: $(PROGRAM_FILE) INPUT.h5 OUTPUT_DIRECTORY [edge_margin] [crop_side]")

using GLMakie
using FFTW
using HDF5
using Statistics

GLMakie.activate!(visible=false)

input = abspath(ARGS[1])
output = abspath(ARGS[2])
edge_margin = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
crop_side_requested = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 128
isfile(input) || error("input HDF5 does not exist: $input")
edge_margin >= 0 || error("edge_margin must be nonnegative")
crop_side_requested > 0 || error("crop_side must be positive")
mkpath(output)

attrs_as_dict(file) = Dict(String(key) => attrs(file)[key] for key in keys(attrs(file)))
metadata = h5open(input, "r") do file
    attrs_as_dict(file)
end

has_density = h5open(input, "r") do file
    haskey(file, "density")
end

if !has_density
    # A large sampled export intentionally cannot support a full heatmap or
    # FFT. Render the samples so the file remains immediately inspectable.
    samples = h5open(input, "r") do file
        (x=read(file, "sample_x"), y=read(file, "sample_y"), density=read(file, "sample_density"))
    end
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1]; xlabel="x", ylabel="y", title="Sampled charge density")
    scatter!(ax, samples.x, samples.y; color=samples.density, colormap=:viridis,
        markersize=5)
    Colorbar(fig[1, 2]; limits=extrema(samples.density), label="n(x,y)")
    save(joinpath(output, "sampled_density.png"), fig)
    save(joinpath(output, "sampled_density.pdf"), fig)
    open(joinpath(output, "README.txt"), "w") do io
        println(io, "This HDF5 file contains sampled density values only; no full FFT was rendered.")
        println(io, "Use export_charge_hdf5.jl with mode=dense for a publication raster.")
    end
    println("Rendered sampled density only: $output")
    exit(0)
end

density = h5open(input, "r") do file
    Float64.(read(file, "density"))
end
ndims(density) == 2 && size(density, 1) == size(density, 2) ||
    error("density dataset must be a square two-dimensional raster")
side = size(density, 1)
edge_margin * 2 < side || error("edge margin leaves no bulk region")
crop_side = min(crop_side_requested, side)
crop_start = (side - crop_side) ÷ 2 + 1
crop_stop = crop_start + crop_side - 1

raw = density .- 0.5
demodulated = similar(raw)
for y in 1:side, x in 1:side
    demodulated[x, y] = iseven((x - 1) + (y - 1)) ? raw[x, y] : -raw[x, y]
end
demodulated .-= mean(demodulated)

bulk_scale_field = copy(demodulated)
if edge_margin > 0
    bulk_scale_field[1:edge_margin, :] .= NaN
    bulk_scale_field[(side - edge_margin + 1):side, :] .= NaN
    bulk_scale_field[:, 1:edge_margin] .= NaN
    bulk_scale_field[:, (side - edge_margin + 1):side] .= NaN
end
bulk_interior = demodulated[(edge_margin + 1):(side - edge_margin),
                             (edge_margin + 1):(side - edge_margin)]
bulk_scale = maximum(abs, bulk_interior)
bulk_scale = bulk_scale == 0 ? 1.0 : bulk_scale

function box_average(field::AbstractMatrix{<:Real}, radius::Int)
    radius >= 0 || error("box radius must be nonnegative")
    n, m = size(field)
    prefix = zeros(Float64, n + 1, m + 1)
    prefix[2:end, 2:end] .= cumsum(cumsum(Float64.(field); dims=1); dims=2)
    result = similar(Float64.(field))
    for y in 1:m, x in 1:n
        x0, x1 = max(1, x - radius), min(n, x + radius)
        y0, y1 = max(1, y - radius), min(m, y + radius)
        result[x, y] = (prefix[x1 + 1, y1 + 1] - prefix[x0, y1 + 1] -
                        prefix[x1 + 1, y0] + prefix[x0, y0]) / ((x1 - x0 + 1) * (y1 - y0 + 1))
    end
    return result
end

envelope = box_average(abs.(demodulated), 8)
window = [0.5 * (1 - cos(2pi * (x - 1) / max(side - 1, 1))) *
          0.5 * (1 - cos(2pi * (y - 1) / max(side - 1, 1)))
          for x in 1:side, y in 1:side]
windowed = demodulated .* window
fft_shifted = circshift(fft(windowed), (side ÷ 2, side ÷ 2))
power = abs2.(fft_shifted)
power ./= max(maximum(power), eps(Float64))
log_power = log10.(max.(power, 1e-12))
k = 2pi .* ((-side ÷ 2):(side ÷ 2 - 1)) ./ side

function heatmap_panel!(fig, slot, field, title, colormap, colorrange)
    ax = Axis(fig[slot[1], slot[2]]; aspect=DataAspect(), title=title, xlabel="x", ylabel="y")
    heatmap!(ax, 1:size(field, 1), 1:size(field, 2), field';
        colormap=colormap, colorrange=colorrange, nan_color=:white)
    Colorbar(fig[slot[1], slot[2] + 1]; limits=colorrange, colormap=colormap)
    return ax
end

fig = Figure(size=(1800, 1250), figure_padding=18)
heatmap_panel!(fig, (1, 1), raw, "Full density deviation", :RdBu,
    (-maximum(abs, raw), maximum(abs, raw)))
heatmap_panel!(fig, (1, 3), bulk_scale_field, "Demodulated density (edge masked; full lattice)", :RdBu,
    (-bulk_scale, bulk_scale))
heatmap_panel!(fig, (1, 5), envelope, "Checkerboard envelope", :viridis,
    (minimum(envelope), maximum(envelope)))
heatmap_panel!(fig, (2, 1), raw[crop_start:crop_stop, crop_start:crop_stop],
    "Central $(crop_side)×$(crop_side) density", :RdBu,
    (-maximum(abs, raw), maximum(abs, raw)))
heatmap_panel!(fig, (2, 3), demodulated[crop_start:crop_stop, crop_start:crop_stop],
    "Central demodulated density", :RdBu, (-bulk_scale, bulk_scale))
heatmap_panel!(fig, (2, 5), log_power, "Hann FFT: log₁₀ relative power", :magma,
    (-12.0, 0.0))

save(joinpath(output, "charge_report.png"), fig; px_per_unit=2)

# The full Brillouin-zone map is useful for checking all spectral weight, but
# finite open-boundary leakage can conceal the physically relevant satellites.
# Save a contrast-enhanced central view of the *same* demodulated, windowed
# transform for presentation use. No Fourier coefficients are discarded.
zoom = min(1.5, pi)
zoomfig = Figure(size=(1050, 900), figure_padding=18)
zoomax = Axis(zoomfig[1, 1], aspect=DataAspect(),
    title="Hann FFT: central reciprocal-space detail", xlabel=L"k_x", ylabel=L"k_y",
    xticks=([-1, 0, 1], [L"-1", L"0", L"1"]),
    yticks=([-1, 0, 1], [L"-1", L"0", L"1"]))
zoommap = heatmap!(zoomax, k, k, log_power'; colormap=:magma,
    colorrange=(-6.0, 0.0))
xlims!(zoomax, -zoom, zoom)
ylims!(zoomax, -zoom, zoom)
Colorbar(zoomfig[1, 2], zoommap; label=L"\log_{10} S_n(\mathbf{k})")
save(joinpath(output, "fourier_hann_central_zoom.png"), zoomfig; px_per_unit=2)

peaks = Tuple{Float64, Int, Int}[]
center = side ÷ 2 + 1
for ky in 1:side, kx in 1:side
    abs(kx - center) <= 1 && abs(ky - center) <= 1 && continue
    push!(peaks, (power[kx, ky], kx - center, ky - center))
end
sort!(peaks; by=first, rev=true)
open(joinpath(output, "fft_peaks.csv"), "w") do io
    println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
    for (rank, (value, kx, ky)) in enumerate(peaks[1:min(64, length(peaks))])
        println(io, "$rank,$kx,$ky,$(2pi * kx / side),$(2pi * ky / side),$value")
    end
end

open(joinpath(output, "metadata.txt"), "w") do io
    println(io, "source=$input")
    println(io, "side=$side")
    println(io, "edge_margin=$edge_margin")
    println(io, "crop_side=$crop_side")
    println(io, "fft_window=separable Hann window on full lattice")
    println(io, "demodulation=((-1)^(x+y))*(n(x,y)-0.5), then mean subtraction")
    for (key, value) in sort!(collect(metadata); by=first)
        println(io, "source_$key=$value")
    end
end
println("Rendered charge report: $output")
