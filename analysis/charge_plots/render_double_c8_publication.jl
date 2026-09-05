#!/usr/bin/env julia

"""Create a publication-style GLMakie figure for a double-C8 density export.

Usage:
  julia --project=analysis/charge_plots render_double_c8_publication.jl INPUT.h5 OUTPUT_DIR [edge_margin]

The edge-masked panel retains the complete lattice; only the outer margin is
excluded from the color-range estimate. Hopping modulation is evaluated with
coordinates centered at the lattice midpoint.
"""

length(ARGS) in (2, 3) || error("usage: $(PROGRAM_FILE) INPUT.h5 OUTPUT_DIR [edge_margin]")
input = abspath(ARGS[1])
output = abspath(ARGS[2])
edge = length(ARGS) == 3 ? parse(Int, ARGS[3]) : 10
isfile(input) || error("input does not exist: $input")
ispath(output) && error("refusing to overwrite existing output: $output")
edge >= 0 || error("edge margin must be nonnegative")

using FFTW
using GLMakie
using HDF5
using Statistics
using TOML

# Render off-screen so the script works from a terminal or batch job while
# still using the GLMakie backend.
GLMakie.activate!(visible=false)

density = h5open(input, "r") do file
    Float64.(read(file["density"]))
end
ndims(density) == 2 && size(density, 1) == size(density, 2) ||
    error("expected square /density dataset")
side = size(density, 1)
2edge < side || error("edge margin leaves no interior")

raw = density .- 0.5
demod = [(-1.0)^(x + y) * raw[x, y] for x in 1:side, y in 1:side]
demod .-= mean(demod)

masked = copy(demod)
if edge > 0
    masked[1:edge, :] .= NaN
    masked[(side - edge + 1):side, :] .= NaN
    masked[:, 1:edge] .= NaN
    masked[:, (side - edge + 1):side] .= NaN
end
interior = demod[(edge + 1):(side - edge), (edge + 1):(side - edge)]
demod_scale = max(maximum(abs, interior), eps(Float64))

tau_fast = sqrt(2.0) / 8.0
tau_slow = (41.0 / 40.0) * tau_fast
angles = (0.0, pi / 4, pi / 2, 3pi / 4)
c8(x, y, tau) = sum(cos(2pi * tau * (x * cos(a) + y * sin(a))) for a in angles) / 4
modulation(x, y) = 0.05 * c8(x, y, tau_fast) + 0.05 * c8(x, y, tau_slow)
center = (side + 1) / 2
delta_hopping = [0.5 * (modulation(x - center + 0.5, y - center) +
                         modulation(x - center, y - center + 0.5))
                  for x in 1:side, y in 1:side]

window = [0.5 * (1 - cos(2pi * (x - 1) / max(side - 1, 1))) *
          0.5 * (1 - cos(2pi * (y - 1) / max(side - 1, 1)))
          for x in 1:side, y in 1:side]
fft_shifted = circshift(fft(demod .* window), (side ÷ 2, side ÷ 2))
power = abs2.(fft_shifted)
power ./= max(maximum(power), eps(Float64))
log_power = log10.(max.(power, 1e-12))

function panel!(fig, row, col, field, title, cmap, limits; label="")
    ax = Axis(fig[row, col]; aspect=DataAspect(), title=title,
              xlabel="x", ylabel="y")
    heatmap!(ax, 1:side, 1:side, field'; colormap=cmap,
             colorrange=limits, nan_color=:white)
    Colorbar(fig[row, col + 1]; limits=limits, colormap=cmap, label=label)
    return ax
end

fig = Figure(size=(1900, 1000), figure_padding=24, fontsize=22)
panel!(fig, 1, 1, density, "Charge density  n(x,y)", :viridis,
       (minimum(density), maximum(density)); label="n")
panel!(fig, 1, 3, masked, "Demodulated density (full lattice; edges masked)", :RdBu,
       (-demod_scale, demod_scale); label="(-1)^(x+y)(n-1/2)")
panel!(fig, 1, 5, delta_hopping, "Centered hopping modulation  δt(x,y)", :RdBu,
       (-maximum(abs, delta_hopping), maximum(abs, delta_hopping)); label="δt")
panel!(fig, 2, 1, log_power, "Hann-windowed FFT of demodulated density", :magma,
       (-12.0, 0.0); label="log₁₀ relative power")
ax = Axis(fig[2, 3]; aspect=DataAspect(), title="Centered coordinates",
          xlabel="x - (L+1)/2", ylabel="y - (L+1)/2")
heatmap!(ax, (1:side) .- center, (1:side) .- center, delta_hopping';
         colormap=:RdBu, colorrange=(-maximum(abs, delta_hopping), maximum(abs, delta_hopping)),
         nan_color=:white)
Colorbar(fig[2, 4]; limits=(-maximum(abs, delta_hopping), maximum(abs, delta_hopping)),
         colormap=:RdBu, label="δt")
Label(fig[2, 5], "τ₁=$(round(tau_fast, sigdigits=6))\nτ₂=$(round(tau_slow, sigdigits=6))\nV₁=V₂=0.05";
      tellheight=false, tellwidth=false, halign=:left, valign=:top)

tmp = output * ".tmp.$(getpid())"
mkpath(tmp)
try
    save(joinpath(tmp, "double_c8_publication.png"), fig; px_per_unit=2)
    center_k = side ÷ 2 + 1
    peaks = Tuple{Float64, Int, Int}[]
    for ky in 1:side, kx in 1:side
        abs(kx - center_k) <= 1 && abs(ky - center_k) <= 1 && continue
        push!(peaks, (power[kx, ky], kx - center_k, ky - center_k))
    end
    sort!(peaks; by=first, rev=true)
    open(joinpath(tmp, "fft_peaks.csv"), "w") do io
        println(io, "rank,kx_index,ky_index,kx,ky,relative_power")
        for (rank, (p, kx, ky)) in enumerate(peaks[1:min(64, length(peaks))])
            println(io, "$rank,$kx,$ky,$(2pi*kx/side),$(2pi*ky/side),$p")
        end
    end
    open(joinpath(tmp, "metadata.toml"), "w") do io
        TOML.print(io, Dict("source_hdf5" => input, "side" => side,
            "edge_margin" => edge, "fft_window" => "separable Hann",
            "hopping_coordinates" => "centered at (L+1)/2 with bond offsets",
            "tau_fast" => tau_fast, "tau_slow" => tau_slow,
            "V_fast" => 0.05, "V_slow" => 0.05))
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote GLMakie publication figure to $output")
