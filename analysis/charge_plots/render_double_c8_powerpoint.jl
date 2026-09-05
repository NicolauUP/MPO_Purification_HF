#!/usr/bin/env julia

"""Render transparent, PowerPoint-oriented views of a double-C8 density.

Usage:
  julia --project=analysis/charge_plots render_double_c8_powerpoint.jl INPUT.h5 OUTPUT_DIR [edge_margin]

The input contains the dense site density exported from the MPO checkpoint.
The bare bond fields are reconstructed from the *same* 2048-site campaign:
coordinates are zero based and the two directions are kept separate. The
figures use Makie's raster `image!` primitive rather than `heatmap!`.
"""

length(ARGS) in (2, 3) || error("usage: $(PROGRAM_FILE) INPUT.h5 OUTPUT_DIR [edge_margin]")
input = abspath(ARGS[1])
output = abspath(ARGS[2])
edge = length(ARGS) == 3 ? parse(Int, ARGS[3]) : 10
isfile(input) || error("input does not exist: $input")
ispath(output) && error("refusing to overwrite existing output: $output")
edge >= 0 || error("edge margin must be nonnegative")

using GLMakie
using HDF5
using Statistics
using TOML

GLMakie.activate!(visible=false)

density = h5open(input, "r") do file
    Float64.(read(file["density"]))
end
ndims(density) == 2 && size(density, 1) == size(density, 2) ||
    error("expected square /density dataset")
side = size(density, 1)
edge * 2 < side || error("edge margin leaves no interior")

# This is the exact double-C8 field used by the 2048² campaign. Each cosine
# supplies the opposite wavevector automatically, so four angles make an
# eightfold star. The bond arguments below follow SquareModel's zero-based
# (x,y) convention and are evaluated at bond midpoints.
tau_fast = sqrt(2.0) / 8.0 + 1.0 / 1024.0
tau_slow = sqrt(2.0) / 8.0 - 1.0 / 1024.0
angles = (0.0, pi / 4, pi / 2, 3pi / 4)
c8(x, y, tau) = sum(cos(2pi * tau * (x * cos(a) + y * sin(a))) for a in angles) / 4
modulation(x, y) = 0.05 * c8(x, y, tau_fast) + 0.05 * c8(x, y, tau_slow)

tx_value(x, y) = -1.0 - modulation(x + 0.5, y)
ty_value(x, y) = -1.0 - modulation(x, y + 0.5)

delta = density .- 0.5
interior = delta[(edge + 1):(side - edge), (edge + 1):(side - edge)]
density_scale = max(maximum(abs, interior), eps(Float64))
demod = [(-1.0)^(x + y) * delta[x, y] for x in 1:side, y in 1:side]
demod .-= mean(demod)
demod_interior = demod[(edge + 1):(side - edge), (edge + 1):(side - edge)]
demod_scale = max(maximum(abs, demod_interior), eps(Float64))
crop_side = min(128, side)
center_start = 1 + (side - crop_side) ÷ 2
quarter_start = clamp(1 + side ÷ 4 - crop_side ÷ 2, 1, side - crop_side + 1)
starts = (center_start, quarter_start)
labels = ("center", "quarter")

transparent_axis!(fig, slot; title, xlabel="x", ylabel="y") = Axis(fig[slot...];
    aspect=DataAspect(), title=title, xlabel=xlabel, ylabel=ylabel,
    backgroundcolor=:transparent, xgridvisible=false, ygridvisible=false,
    topspinevisible=false, rightspinevisible=false)

function save_single(path, field, title, cmap, limits; xlabel="x", ylabel="y", xvalues=1:side, yvalues=1:side)
    fig = Figure(size=(1100, 950), backgroundcolor=:white, figure_padding=24, fontsize=26)
    ax = transparent_axis!(fig, (1, 1); title, xlabel, ylabel)
    image!(ax, (first(xvalues), last(xvalues)), (first(yvalues), last(yvalues)), field';
           colormap=cmap, colorrange=limits, interpolate=false)
    Colorbar(fig[1, 2]; limits=limits, colormap=cmap, label=title, vertical=true)
    save(path, fig; px_per_unit=2)
end

function save_two(path, fields, titles, cmap, limits; ranges=nothing)
    fig = Figure(size=(1900, 900), backgroundcolor=:white, figure_padding=24, fontsize=24)
    for j in 1:2
        xr = isnothing(ranges) ? (1:side, 1:side) : ranges[j]
        ax = transparent_axis!(fig, (1, 2j - 1); title=titles[j])
        image!(ax, (first(xr[1]), last(xr[1])), (first(xr[2]), last(xr[2])), fields[j]';
               colormap=cmap, colorrange=limits, interpolate=false)
        Colorbar(fig[1, 2j]; limits=limits, colormap=cmap, label=titles[j], vertical=true)
    end
    save(path, fig; px_per_unit=2)
end

mkpath(output * ".tmp.$(getpid())")
tmp = output * ".tmp.$(getpid())"
try
    # Full site density, diverging around n=0.5 (white at half filling).
    save_single(joinpath(tmp, "density_full_2048.png"), density,
        "n(x,y)", :RdBu, (0.5 - density_scale, 0.5 + density_scale))
    save_single(joinpath(tmp, "density_demodulated_full_2048.png"), demod,
        "(-1)^(x+y)(n-1/2)", :RdBu, (-demod_scale, demod_scale))

    # Two bulk 128x128 views with identical density color limits.
    crops = [density[s:s + crop_side - 1, s:s + crop_side - 1] for s in starts]
    fig = Figure(size=(1900, 900), backgroundcolor=:transparent, figure_padding=24, fontsize=24)
    for j in 1:2
        s = starts[j]
        ax = transparent_axis!(fig, (1, 2j - 1); title="Density $(labels[j]) ($(crop_side)×$(crop_side))")
        image!(ax, (s, s + crop_side - 1), (s, s + crop_side - 1), collect(crops[j]');
               colormap=:RdBu, colorrange=(0.5 - density_scale, 0.5 + density_scale), interpolate=false)
        Colorbar(fig[1, 2j]; limits=(0.5 - density_scale, 0.5 + density_scale),
                 colormap=:RdBu, label="n(x,y)", vertical=true)
    end
    save(joinpath(tmp, "density_crops_128_center_quarter.png"), fig; px_per_unit=2)
    demod_crops = [demod[s:s + crop_side - 1, s:s + crop_side - 1] for s in starts]
    save_two(joinpath(tmp, "density_demodulated_crops_128_center_quarter.png"), demod_crops,
        ("center demodulated density", "quarter demodulated density"), :RdBu,
        (-demod_scale, demod_scale),
        ranges=((starts[1]:(starts[1] + crop_side - 1), starts[1]:(starts[1] + crop_side - 1)),
                (starts[2]:(starts[2] + crop_side - 1), starts[2]:(starts[2] + crop_side - 1))))

    # Separate bare x/y bond-modulation fields, not their misleading average.
    # Only the requested 128x128 bulk views are materialized; a full pair of
    # 2048² bond arrays would duplicate hundreds of MB for no presentation gain.
    bond_crops = Vector{Tuple{Matrix{Float64},Matrix{Float64}}}(undef, 2)
    for j in 1:2
        s = starts[j]
        bond_crops[j] = (
            [tx_value(x - 1, y - 1) + 1.0 for x in s:(s + crop_side - 1), y in s:(s + crop_side - 1)],
            [ty_value(x - 1, y - 1) + 1.0 for x in s:(s + crop_side - 1), y in s:(s + crop_side - 1)],
        )
    end
    bond_scale = max(maximum(abs, bond_crops[1][1]), maximum(abs, bond_crops[1][2]),
                      maximum(abs, bond_crops[2][1]), maximum(abs, bond_crops[2][2]))
    fig = Figure(size=(1900, 1700), backgroundcolor=:white, figure_padding=24, fontsize=24)
    for j in 1:2, k in 1:2
        field = bond_crops[j][k]
        ax = transparent_axis!(fig, (j, 2k - 1); title="$(labels[j]) δt$(k == 1 ? "ₓ" : "ᵧ") ($(crop_side)×$(crop_side))")
        s = starts[j]
        image!(ax, (s, s + crop_side - 1), (s, s + crop_side - 1), field';
               colormap=:RdBu, colorrange=(-bond_scale, bond_scale), interpolate=false)
        Colorbar(fig[j, 2k]; limits=(-bond_scale, bond_scale), colormap=:RdBu,
                 label="δt", vertical=true)
    end
    save(joinpath(tmp, "bond_modulation_crops_128_center_quarter.png"), fig; px_per_unit=2)

    # Correlations between the bare hopping and density at each bond endpoint.
    # Do not use only (n_i+n_j)/2: for a staggered state that average can cancel
    # the physically relevant response. Statistics use every bond; plotting is
    # deterministically subsampled to keep the PNG light.
    mutable struct CorrelationStats
        count::Int
        sx::Float64
        sy::Float64
        sxx::Float64
        syy::Float64
        sxy::Float64
    end
    CorrelationStats() = CorrelationStats(0, 0.0, 0.0, 0.0, 0.0, 0.0)
    function add!(a::CorrelationStats, x, y)
        a.count += 1
        a.sx += x; a.sy += y; a.sxx += x*x; a.syy += y*y; a.sxy += x*y
        return nothing
    end
    function correlation(a::CorrelationStats)
        mx, my = a.sx / a.count, a.sy / a.count
        cov = a.sxy / a.count - mx * my
        vx = a.sxx / a.count - mx^2
        vy = a.syy / a.count - my^2
        return cov / sqrt(max(vx * vy, eps(Float64))), cov / max(vx, eps(Float64))
    end
    total_bonds = 2 * side * (side - 1)
    plot_cap = min(400_000, total_bonds)
    stride = max(1, cld(total_bonds, plot_cap))
    xs = Float64[]
    y_origin = Float64[]; y_neighbour = Float64[]
    y_demod_origin = Float64[]; y_demod_neighbour = Float64[]
    sizehint!(xs, cld(total_bonds, stride))
    for ys in (y_origin, y_neighbour, y_demod_origin, y_demod_neighbour)
        sizehint!(ys, cld(total_bonds, stride))
    end
    origin_stats = CorrelationStats(); neighbour_stats = CorrelationStats()
    demod_origin_stats = CorrelationStats(); demod_neighbour_stats = CorrelationStats()
    average_stats = CorrelationStats()
    count = 0
    sampled = 0
    for x in 0:(side - 2), y in 0:(side - 1)
        t = tx_value(x, y)
        no, nn = density[x + 1, y + 1], density[x + 2, y + 1]
        do_, dn = demod[x + 1, y + 1], demod[x + 2, y + 1]
        add!(origin_stats, t, no); add!(neighbour_stats, t, nn)
        add!(demod_origin_stats, t, do_); add!(demod_neighbour_stats, t, dn)
        add!(average_stats, t, 0.5 * (no + nn))
        count += 1
        if (count - 1) % stride == 0
            push!(xs, t); push!(y_origin, no); push!(y_neighbour, nn)
            push!(y_demod_origin, do_); push!(y_demod_neighbour, dn); sampled += 1
        end
    end
    for x in 0:(side - 1), y in 0:(side - 2)
        t = ty_value(x, y)
        no, nn = density[x + 1, y + 1], density[x + 1, y + 2]
        do_, dn = demod[x + 1, y + 1], demod[x + 1, y + 2]
        add!(origin_stats, t, no); add!(neighbour_stats, t, nn)
        add!(demod_origin_stats, t, do_); add!(demod_neighbour_stats, t, dn)
        add!(average_stats, t, 0.5 * (no + nn))
        count += 1
        if (count - 1) % stride == 0
            push!(xs, t); push!(y_origin, no); push!(y_neighbour, nn)
            push!(y_demod_origin, do_); push!(y_demod_neighbour, dn); sampled += 1
        end
    end
    stats = (origin=correlation(origin_stats), neighbour=correlation(neighbour_stats),
             demod_origin=correlation(demod_origin_stats), demod_neighbour=correlation(demod_neighbour_stats),
             average=correlation(average_stats))
    fig = Figure(size=(1900, 1700), backgroundcolor=:white, figure_padding=24, fontsize=24)
    endpoint_data = ((y_origin, "origin nᵢ", stats.origin),
                     (y_neighbour, "neighbour nⱼ", stats.neighbour),
                     (y_demod_origin, "origin demodulated density", stats.demod_origin),
                     (y_demod_neighbour, "neighbour demodulated density", stats.demod_neighbour))
    for j in 1:4
        row, col = j <= 2 ? (1, 2j - 1) : (2, 2(j - 2) - 1)
        ys, title, (r, sl) = endpoint_data[j]
        ylabel = j <= 2 ? "density" : "demodulated density"
        ax = transparent_axis!(fig, (row, col); title, xlabel="initial tᵢⱼ", ylabel)
        scatter!(ax, xs, ys; markersize=2, color=(:steelblue, 0.18), strokewidth=0)
        Label(fig[row, col], "all bonds: $(count)   plotted: $(sampled)\nPearson r=$(round(r, sigdigits=5)), slope=$(round(sl, sigdigits=5))";
            tellheight=false, tellwidth=false, halign=:left, valign=:top, padding=(70, 0, 0, 70))
    end
    save(joinpath(tmp, "density_vs_initial_hopping_scatter.png"), fig; px_per_unit=2)

    open(joinpath(tmp, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "source_hdf5" => input, "side" => side, "edge_margin" => edge,
            "crop_side" => crop_side, "crop_starts_1based" => collect(starts),
            "tau_fast" => tau_fast, "tau_slow" => tau_slow,
            "bond_coordinates" => "SquareModel zero-based coordinates; midpoint x/y bonds",
            "density_bond_value" => "endpoint densities reported separately; endpoint average retained only as a reference",
            "correlation_bonds" => count, "scatter_points" => sampled,
            "pearson_r_origin" => stats.origin[1], "pearson_r_neighbour" => stats.neighbour[1],
            "pearson_r_demod_origin" => stats.demod_origin[1],
            "pearson_r_demod_neighbour" => stats.demod_neighbour[1],
            "pearson_r_endpoint_average_reference" => stats.average[1],
            "image_primitive" => "Makie image!", "background" => "white canvas; transparent axes"))
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end
println("Wrote PowerPoint figures to $output")
