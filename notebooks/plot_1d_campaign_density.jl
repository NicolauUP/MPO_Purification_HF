### A Pluto.jl notebook ###
# v0.20.13

# This notebook is intentionally an analysis-only tool. It does not activate
# or modify the solver environment. Run it locally after downloading a whole
# campaign result directory from the cluster. Install CairoMakie in a separate
# local analysis environment if it is not already available:
#
#   mkdir -p analysis_env
#   julia --project=analysis_env -e 'using Pkg; Pkg.add("CairoMakie")'
#   julia --project=analysis_env -e 'using Pluto; Pluto.run()'
#
# Then open this file in Pluto and set `result_root` below. Only Julia stdlib
# modules plus CairoMakie are used; CSV parsing is deliberately dependency-free.

using Markdown
using InteractiveUtils

# ╔═╡ 0f72f5fe-2f6d-4c97-b46b-c8d92407cc70
begin
    using CairoMakie
    using Statistics
    using TOML
end

# ╔═╡ 3c3e762b-417e-4a16-bd2b-2e36de18344b
md"""
# 1D charge-density comparison

This notebook compares every `site_density.csv` found directly below one
downloaded campaign result directory. It also reports particle number, final
energy, numerical residuals, and bulk staggered density order.

The current reference campaign is expected at a directory named
`aubry_andre_nn_l10_seed0p1`, containing task directories such as
`task_0001_v2_0_u_1`.
"""

# ╔═╡ 1434ab81-6fec-4ab1-97bb-1c5adacb788b
result_root = raw"/CHANGE/THIS/TO/aubry_andre_nn_l10_seed0p1"

# ╔═╡ f8ccf7d5-75b2-4b7d-a975-60ef3ad53830
begin
    function _unquote(value::AbstractString)
        strip(strip(value), '"')
    end

    function read_site_density(path::AbstractString)
        lines = readlines(path)
        isempty(lines) && error("empty density CSV: $path")
        strip(lines[1]) == "\"site\",\"density\"" || error(
            "unexpected site-density header in $path: $(lines[1])",
        )
        sites = Int[]
        density = Float64[]
        for line in Iterators.drop(lines, 1)
            isempty(strip(line)) && continue
            columns = split(strip(line), ','; limit=2)
            length(columns) == 2 || error("malformed density row in $path: $line")
            push!(sites, parse(Int, _unquote(columns[1])))
            push!(density, parse(Float64, _unquote(columns[2])))
        end
        issorted(sites) || error("site indices are not sorted in $path")
        sites == collect(1:length(sites)) || error("site indices are not contiguous in $path")
        return sites, density
    end

    function bulk_staggered_order(sites::Vector{Int}, density::Vector{Float64}; edge_fraction=0.1)
        length(sites) == length(density) || error("site and density lengths differ")
        count = length(sites)
        first_bulk = floor(Int, edge_fraction * count) + 1
        last_bulk = count - floor(Int, edge_fraction * count)
        bulk = first_bulk:last_bulk
        return mean((isodd(sites[index]) ? 1.0 : -1.0) * (density[index] - 0.5) for index in bulk)
    end

    function load_campaign(root::AbstractString)
        isdir(root) || error("result_root is not a directory: $root")
        task_dirs = sort(filter(path -> isdir(path) && isfile(joinpath(path, "site_density.csv")),
            readdir(root; join=true)))
        isempty(task_dirs) && error("no task directories containing site_density.csv found in $root")
        return map(task_dirs) do directory
            sites, density = read_site_density(joinpath(directory, "site_density.csv"))
            observables_path = joinpath(directory, "observables.toml")
            observables = isfile(observables_path) ? TOML.parsefile(observables_path) : Dict{String,Any}()
            (
                label=replace(basename(directory), r"^task_\\d+_" => ""),
                directory=directory,
                sites=sites,
                density=density,
                observables=observables,
                bulk_staggered_order=bulk_staggered_order(sites, density),
            )
        end
    end
end

# ╔═╡ 271d1323-9200-4df1-8fbd-2ffdf932744a
cases = load_campaign(result_root)

# ╔═╡ 9c6179fc-1c90-4491-b1bd-098c5b99a71d
summary = [
    (
        case=entry.label,
        N=length(entry.sites),
        particle_number=get(entry.observables, "particle_number", missing),
        energy_total=get(entry.observables, "energy_total", missing),
        idempotency=get(entry.observables, "idempotency_residual", missing),
        stationarity=get(entry.observables, "stationarity_residual", missing),
        bulk_staggered_order=entry.bulk_staggered_order,
    )
    for entry in cases
]

# ╔═╡ d1ad9bf0-3191-42ae-8300-221e8b6495f5
md"""
## Full density and central-window view

The upper panel exposes boundary response. The lower panel magnifies a central
window; persistent even/odd splitting there is a bulk Hartree--Fock CDW
branch, while an amplitude decaying toward the centre is boundary Friedel
response.
"""

# ╔═╡ 5b9b7c4d-ebfb-4a7d-bf24-0e61c0fde728
begin
    CairoMakie.activate!()
    figure = Figure(size=(1100, 760))
    full_axis = Axis(
        figure[1, 1], xlabel="site", ylabel="charge density",
        title="Full 1D density profile",
    )
    center_axis = Axis(
        figure[2, 1], xlabel="site", ylabel="charge density",
        title="Central density window",
    )
    for entry in cases
        lines!(full_axis, entry.sites, entry.density; label=entry.label, linewidth=1.8)
        center = cld(length(entry.sites), 2)
        half_width = min(80, center - 1, length(entry.sites) - center)
        window = (center - half_width):(center + half_width)
        lines!(center_axis, entry.sites[window], entry.density[window]; label=entry.label, linewidth=1.8)
    end
    axislegend(full_axis; position=:rb)
    axislegend(center_axis; position=:rb)
    figure
end

# ╔═╡ a6032f2d-1c62-4a2b-ac01-3d762e66727f
md"""
## Optional export

Run the following cell only if you want a PNG next to the downloaded result
directory. The notebook itself does not alter any solver output.
"""

# ╔═╡ 37b86430-39b2-4bc3-9d77-77fab209472f
# save(joinpath(result_root, "charge_density_comparison.png"), figure)

# ╔═╡ Cell order:
# ╟─0f72f5fe-2f6d-4c97-b46b-c8d92407cc70
# ╟─3c3e762b-417e-4a16-bd2b-2e36de18344b
# ╠═1434ab81-6fec-4ab1-97bb-1c5adacb788b
# ╠═f8ccf7d5-75b2-4b7d-a975-60ef3ad53830
# ╠═271d1323-9200-4df1-8fbd-2ffdf932744a
# ╠═9c6179fc-1c90-4491-b1bd-098c5b99a71d
# ╟─d1ad9bf0-3191-42ae-8300-221e8b6495f5
# ╠═5b9b7c4d-ebfb-4a7d-bf24-0e61c0fde728
# ╟─a6032f2d-1c62-4a2b-ac01-3d762e66727f
# ╠═37b86430-39b2-4bc3-9d77-77fab209472f
