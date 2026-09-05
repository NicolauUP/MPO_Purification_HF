#!/usr/bin/env julia

length(ARGS) in (2, 3) || error("usage: $(PROGRAM_FILE) INPUT.h5 OUTPUT.png [crop_side]")
input, output = abspath.(ARGS[1:2])
crop_side = length(ARGS) == 3 ? parse(Int, ARGS[3]) : 128

using GLMakie
using HDF5
GLMakie.activate!(visible=false)
density = h5open(input, "r") do f
    Float64.(read(f["density"]))
end
side = size(density, 1)
crop_side <= side || error("crop is larger than density")
starts = (1 + (side - crop_side) ÷ 2, 1 + side ÷ 4 - crop_side ÷ 2)
scale = maximum(abs, density .- 0.5)
fig = Figure(size=(1900, 900), backgroundcolor=:white, figure_padding=24, fontsize=24)
for (j, s) in enumerate(starts)
    ax = Axis(fig[1, 2j - 1], aspect=DataAspect(),
        title=(j == 1 ? "center" : "quarter") * " density ($(crop_side)×$(crop_side))",
        xlabel="x", ylabel="y", backgroundcolor=:transparent,
        xgridvisible=false, ygridvisible=false, topspinevisible=false,
        rightspinevisible=false)
    crop = density[s:(s + crop_side - 1), s:(s + crop_side - 1)]
    image!(ax, crop; colormap=:RdBu, colorrange=(0.5 - scale, 0.5 + scale), interpolate=false)
    Colorbar(fig[1, 2j]; limits=(0.5 - scale, 0.5 + scale), colormap=:RdBu,
        label="n(x,y)", vertical=true)
end
save(output, fig; px_per_unit=2)
println("Wrote $output")
