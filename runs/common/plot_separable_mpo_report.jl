#!/usr/bin/env julia

"""Render presentation figures from `collect_separable_mpo_report.jl` output.

Usage:
    julia --project=. runs/common/plot_separable_mpo_report.jl REPORT_DIRECTORY

This is deliberately an analysis-side utility. It requires CairoMakie to be
available in the Julia environment used for plotting, but does not alter the
solver project or any simulation result directory.
"""

using TOML
using Statistics

try
    @eval using CairoMakie
catch error_message
    error("CairoMakie is required only for plotting. Activate the Julia analysis environment that contains it. Original error: $error_message")
end

length(ARGS) == 1 || error("usage: plot_separable_mpo_report.jl REPORT_DIRECTORY")
report = abspath(ARGS[1])
isdir(report) || error("report directory does not exist: $report")
metadata = TOML.parsefile(joinpath(report, "metadata.toml"))
figure_directory = joinpath(report, "figures")
mkpath(figure_directory)

function csv_fields(line::AbstractString)
    fields, buffer, quoted = String[], IOBuffer(), false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            next = nextind(line, i)
            if quoted && next <= lastindex(line) && line[next] == '"'
                print(buffer, '"'); i = nextind(line, next); continue
            end
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(fields, String(take!(buffer)))
        elseif c != '\r' && c != '\n'
            print(buffer, c)
        end
        i = nextind(line, i)
    end
    push!(fields, String(take!(buffer)))
    return fields
end

function read_csv(path)
    isfile(path) || return Dict{String,String}[]
    open(path) do io
        eof(io) && return Dict{String,String}[]
        header = csv_fields(readline(io))
        rows = Dict{String,String}[]
        for line in eachline(io)
            isempty(strip(line)) && continue
            fields = csv_fields(line)
            length(fields) == length(header) || error("malformed row in $path")
            push!(rows, Dict(zip(header, fields)))
        end
        rows
    end
end

asfloat(row, field) = begin
    value = tryparse(Float64, get(row, field, ""))
    isnothing(value) ? NaN : value
end
finite_rows(rows, fields) = [row for row in rows if all(isfinite(asfloat(row, field)) for field in fields)]

set_theme!(Theme(
    fontsize=18,
    Axis=(backgroundcolor=:white, xgridcolor=(:gray, 0.18), ygridcolor=(:gray, 0.18)),
))

compression = read_csv(joinpath(report, "compression_ladder.csv"))
sp2 = read_csv(joinpath(report, "sp2_summary.csv"))
if !isempty(compression) || !isempty(sp2)
    figure = Figure(size=(1180, 500))
    ax_density = Axis(figure[1, 1], xlabel="maximum bond dimension χ", ylabel="density RMS error", xscale=log10, yscale=log10)
    ax_idem = Axis(figure[1, 2], xlabel="maximum bond dimension χ", ylabel="trace / idempotency error", xscale=log10, yscale=log10)
    direct = finite_rows(compression, ["max_chi", "density_rms_error"])
    if !isempty(direct)
        sort!(direct, by=row -> asfloat(row, "max_chi"))
        lines!(ax_density, asfloat.(direct, "max_chi"), asfloat.(direct, "density_rms_error"), color=:dodgerblue, marker=:circle, label="exact projector → MPO")
        lines!(ax_idem, asfloat.(direct, "max_chi"), max.(asfloat.(direct, "trace_error"), eps()), color=:dodgerblue, marker=:circle, label="compression trace error")
    end
    fixed_h = finite_rows(sp2, ["requested_maxdim", "density_rms_error"])
    if !isempty(fixed_h)
        sort!(fixed_h, by=row -> asfloat(row, "requested_maxdim"))
        lines!(ax_density, asfloat.(fixed_h, "requested_maxdim"), asfloat.(fixed_h, "density_rms_error"), color=:darkorange, marker=:rect, label="SP2 → exact projector")
    end
    fixed_h_quality = finite_rows(sp2, ["requested_maxdim", "trace_error", "idempotency_residual"])
    if !isempty(fixed_h_quality)
        sort!(fixed_h_quality, by=row -> asfloat(row, "requested_maxdim"))
        χ = asfloat.(fixed_h_quality, "requested_maxdim")
        lines!(ax_idem, χ, max.(asfloat.(fixed_h_quality, "trace_error"), eps()), color=:darkorange, marker=:rect, label="SP2 trace error")
        lines!(ax_idem, χ, max.(asfloat.(fixed_h_quality, "idempotency_residual"), eps()), color=:seagreen, marker=:utriangle, label="SP2 idempotency")
    end
    axislegend(ax_density, position=:lb); axislegend(ax_idem, position=:lb)
    Label(figure[0, :], "Deterministic MPO validation: rank ladder", fontsize=24, font=:bold)
    save(joinpath(figure_directory, "projector_accuracy_vs_rank.pdf"), figure)
end

trajectory = read_csv(joinpath(report, "sp2_trajectory.csv"))
if !isempty(trajectory)
    figure = Figure(size=(1180, 510))
    ax_error = Axis(figure[1, 1], xlabel="SP2 iteration", ylabel="residual", yscale=log10)
    ax_chi = Axis(figure[1, 2], xlabel="SP2 iteration", ylabel="MPO bond dimension χ")
    labels = unique(get.(trajectory, "label", "fixed-H SP2"))
    palette = [:dodgerblue, :darkorange, :seagreen, :purple]
    for (index, label) in enumerate(labels)
        rows = filter(row -> get(row, "label", "") == label, trajectory)
        sort!(rows, by=row -> asfloat(row, "iteration"))
        color = palette[mod1(index, length(palette))]
        iteration = asfloat.(rows, "iteration")
        trace_error = max.(asfloat.(rows, "trace_error"), eps())
        idempotency = max.(asfloat.(rows, "idempotency_residual"), eps())
        lines!(ax_error, iteration, trace_error, color=color, label="$label: trace")
        lines!(ax_error, iteration, idempotency, color=(color, 0.55), linestyle=:dash, label="$label: P²−P")
        lines!(ax_chi, iteration, asfloat.(rows, "rho_max_chi"), color=color, label=label)
    end
    axislegend(ax_error, position=:rt, labelsize=11); axislegend(ax_chi, position=:rb)
    Label(figure[0, :], "Fixed-H SP2 trajectory and rank growth", fontsize=24, font=:bold)
    save(joinpath(figure_directory, "sp2_trajectory.pdf"), figure)
end

scf = read_csv(joinpath(report, "scf_trajectory.csv"))
if !isempty(scf)
    figure = Figure(size=(1180, 510))
    ax_residual = Axis(figure[1, 1], xlabel="SCF iteration", ylabel="field / density residual", yscale=log10)
    ax_energy = Axis(figure[1, 2], xlabel="SCF iteration", ylabel="total energy")
    labels = unique(["$(get(row, "solver", "SCF")): $(get(row, "label", ""))" for row in scf])
    palette = [:mediumpurple, :teal, :firebrick, :goldenrod]
    for (index, label) in enumerate(labels)
        rows = filter(scf) do row
            "$(get(row, "solver", "SCF")): $(get(row, "label", ""))" == label
        end
        sort!(rows, by=row -> asfloat(row, "iteration"))
        color = palette[mod1(index, length(palette))]
        iteration = asfloat.(rows, "iteration")
        residual = max.(asfloat.(rows, "rho_residual"), eps())
        lines!(ax_residual, iteration, residual, color=color, label=label)
        lines!(ax_energy, iteration, asfloat.(rows, "energy_total"), color=color, label=label)
    end
    axislegend(ax_residual, position=:rt); axislegend(ax_energy, position=:rb)
    Label(figure[0, :], "Self-consistent field convergence", fontsize=24, font=:bold)
    save(joinpath(figure_directory, "scf_trajectory.pdf"), figure)
end

# Optional exact/SP2 density comparison heatmaps. The comparison writer stores
# x and y already, therefore this does not rely on an implicit QTT ordering.
density_directory = string(get(metadata, "density_comparison_directory", ""))
density_path = joinpath(density_directory, "site_density_comparison.csv")
if !isempty(density_directory) && isfile(density_path)
    rows = read_csv(density_path)
    nx = maximum(round(Int, asfloat(row, "x")) for row in rows) + 1
    ny = maximum(round(Int, asfloat(row, "y")) for row in rows) + 1
    exact = fill(NaN, nx, ny); sp2_density = fill(NaN, nx, ny); difference = fill(NaN, nx, ny)
    for row in rows
        x, y = round(Int, asfloat(row, "x")) + 1, round(Int, asfloat(row, "y")) + 1
        exact[x, y] = asfloat(row, "exact")
        sp2_density[x, y] = asfloat(row, "sp2")
        difference[x, y] = asfloat(row, "error")
    end
    difference_limit = maximum(abs, difference)
    density_limit = maximum(abs, vcat(vec(exact .- 0.5), vec(sp2_density .- 0.5)))
    figure = Figure(size=(1350, 440))
    for (column, data, title, range, map) in (
        (1, exact, "dense reference density", (0.5 - density_limit, 0.5 + density_limit), :balance),
        (2, sp2_density, "SP2 MPO density", (0.5 - density_limit, 0.5 + density_limit), :balance),
        (3, difference, "SP2 − dense", (-difference_limit, difference_limit), :balance),
    )
        ax = Axis(figure[1, column], title=title, xlabel="x", ylabel="y", aspect=DataAspect())
        image!(ax, 0:nx-1, 0:ny-1, data; colormap=map, colorrange=range)
    end
    Label(figure[0, :], "Density comparison at fixed Hamiltonian", fontsize=24, font=:bold)
    save(joinpath(figure_directory, "fixed_h_density_comparison.pdf"), figure)
end

println("Figures written to $figure_directory")
