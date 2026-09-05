#!/usr/bin/env julia

# Render edge-trimmed pi-demodulated density and clean Hann FFT figures.
using CairoMakie
using FFTW
using HDF5
using Statistics

length(ARGS) == 4 || error("usage: render_conference_figures.jl C8_H5 SEPARABLE_H5 OUTPUT_ROOT EDGE")
c8_path, sep_path, output_root = abspath.(ARGS[1:3])
edge = parse(Int, ARGS[4])

function read_density(path)
    h5open(path, "r") do f
        haskey(f, "density") || error("missing dense density dataset in $path")
        Float64.(read(f["density"]))
    end
end

function render_case(name, path, output_root, edge; robust_density=false, fft_log_floor=-12.0, fft_zoom=1.6, satellite_enhanced=false)
    density = read_density(path)
    nx, ny = size(density)
    nx == ny || error("expected square density in $path")
    2edge < nx || error("edge margin too large")
    out = joinpath(output_root, name)
    mkpath(out)
    raw = density .- 0.5
    stagger = [iseven((x - 1) + (y - 1)) ? 1.0 : -1.0 for x in 1:nx, y in 1:ny]
    bulk_x = (edge + 1):(nx - edge)
    bulk_y = (edge + 1):(ny - edge)
    demod_raw = stagger .* raw
    checkerboard_carrier = mean(@view demod_raw[bulk_x, bulk_y])
    demod = demod_raw .- checkerboard_carrier
    interior = demod[bulk_x, bulk_y]
    scale = robust_density ? quantile(abs.(vec(interior)), 0.995) : maximum(abs, interior)
    scale = max(scale, eps(Float64))
    masked = copy(demod)
    masked[1:edge, :] .= NaN; masked[(nx-edge+1):nx, :] .= NaN
    masked[:, 1:edge] .= NaN; masked[:, (ny-edge+1):ny] .= NaN
    crop = min(256, nx - 2edge)
    x0 = (nx - crop) ÷ 2 + 1; y0 = (ny - crop) ÷ 2 + 1
    xr, yr = x0:(x0 + crop - 1), y0:(y0 + crop - 1)
    fig = Figure(size=(1800, 820), fontsize=24)
    ax1 = Axis(fig[1, 1], title="$name: carrier-removed π-demodulated density", xlabel="x", ylabel="y", aspect=DataAspect())
    hm = heatmap!(ax1, 0:nx-1, 0:ny-1, masked'; colormap=:RdBu, colorrange=(-scale, scale), nan_color=:white)
    ax2 = Axis(fig[1, 2], title="$name: central $crop×$crop crop", xlabel="x", ylabel="y", aspect=DataAspect())
    heatmap!(ax2, xr .- 1, yr .- 1, demod[xr, yr]'; colormap=:RdBu, colorrange=(-scale, scale))
    Colorbar(fig[1, 3], hm, label="δd(x,y): demodulated density minus bulk mean")
    save(joinpath(out, "pi_demodulated_fluctuation_edge$(edge)_crop$(crop).png"), fig; px_per_unit=2)

    win = 0.5 .* (1 .- cos.(2π .* (0:nx-1) ./ max(nx - 1, 1)))
    window = win * win'
    weights = zeros(Float64, nx, ny)
    weights[bulk_x, bulk_y] .= window[bulk_x, bulk_y]
    weighted_mean = sum(demod .* weights) / sum(weights)
    fft_signal = (demod .- weighted_mean) .* weights
    power = abs2.(circshift(fft(fft_signal), (nx ÷ 2, ny ÷ 2)))
    power ./= max(maximum(power), eps(Float64))
    log_power = log10.(max.(power, 1e-12))
    k = 2π .* ((-nx ÷ 2):(nx ÷ 2 - 1)) ./ nx
    figfft = Figure(size=(1100, 950), fontsize=24)
    ax = Axis(figfft[1, 1], title="$name: carrier-removed, edge-trimmed Hann FFT", xlabel="kₓ", ylabel="kᵧ", aspect=DataAspect(), xticks=([-π, 0, π], ["-π", "0", "π"]), yticks=([-π, 0, π], ["-π", "0", "π"]))
    hmfft = heatmap!(ax, k, k, log_power'; colormap=:magma, colorrange=(fft_log_floor, 0))
    Colorbar(figfft[1, 2], hmfft, label="log₁₀ relative power")
    save(joinpath(out, "pi_demodulated_fft_carrier_removed_edge$(edge).png"), figfft; px_per_unit=2)
    zoom = min(fft_zoom, maximum(abs, k))
    xlims!(ax, -zoom, zoom); ylims!(ax, -zoom, zoom)
    save(joinpath(out, "pi_demodulated_fft_carrier_removed_edge$(edge)_zoom.png"), figfft; px_per_unit=2)
    if satellite_enhanced
        # Keep the data and normalization unchanged, but saturate peaks above 10^-2
        # so that weak quasiperiodic satellites remain visible in presentations.
        figsat = Figure(size=(1100, 950), fontsize=24)
        axsat = Axis(
            figsat[1, 1],
            title="$name: satellite-enhanced FFT (primary peaks saturated)",
            xlabel="kₓ",
            ylabel="kᵧ",
            aspect=DataAspect(),
            xticks=([-π, 0, π], ["-π", "0", "π"]),
            yticks=([-π, 0, π], ["-π", "0", "π"]),
        )
        hmsat = heatmap!(axsat, k, k, log_power'; colormap=:magma, colorrange=(-6, -2))
        # Overlay well-separated strongest Fourier pixels. Marker size is a
        # presentation aid only; coordinates and colors retain the FFT data.
        selected = CartesianIndex{2}[]
        for linear_index in sortperm(vec(power); rev=true)
            candidate = CartesianIndices(power)[linear_index]
            log_power[candidate] < -5 && break
            all(maximum(abs.(Tuple(candidate) .- Tuple(previous))) > 3 for previous in selected) || continue
            push!(selected, candidate)
            length(selected) == 256 && break
        end
        scatter!(
            axsat,
            [k[index[1]] for index in selected],
            [k[index[2]] for index in selected];
            color=[log_power[index] for index in selected],
            colormap=:magma,
            colorrange=(-6, -2),
            markersize=7,
            strokewidth=0,
        )
        Colorbar(figsat[1, 2], hmsat, label="log₁₀ relative power (clipped above 10⁻²)")
        save(joinpath(out, "pi_demodulated_fft_satellites_edge$(edge).png"), figsat; px_per_unit=2)
    end
    println("$name: carrier=$(checkerboard_carrier), centered_bulk_mean=$(mean(interior)), weighted_fft_mean=$(weighted_mean)")
    println("$name: density and FFT figures written to $out")
end

render_case("C8_2048x2048", c8_path, output_root, edge)
render_case(
    "separable_2048x2048",
    sep_path,
    output_root,
    edge;
    robust_density=true,
    fft_log_floor=-6.0,
    fft_zoom=π,
    satellite_enhanced=true,
)
