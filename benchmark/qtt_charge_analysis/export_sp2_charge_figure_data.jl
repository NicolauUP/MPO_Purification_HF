#!/usr/bin/env julia

"""Export exact QTT charge data for a two-panel real-space figure.

The full panel retains every one of the 1024x1024 raw density values, so an
alternating checkerboard component is never averaged away. The QTT/MPS is
contracted once into a dense 2^20-amplitude tensor, rather than evaluating one
million independent MPS overlaps; the zoom is then a direct slice.
"""

using Dates
using HDF5
using ITensors, ITensorMPS
using MPO_MeanField
using TOML

length(ARGS) == 4 || error("usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT_DIRECTORY")
campaign_path, result_directory, output = abspath(ARGS[1]), abspath(ARGS[3]), abspath(ARGS[4])
task = parse(Int, ARGS[2])
checkpoint_path = joinpath(result_directory, "converged_state.h5")
isfile(checkpoint_path) || error("missing converged checkpoint: $checkpoint_path")
ispath(output) && error("refusing to overwrite output: $output")

namespace = Module(gensym(:FigureCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
campaign = getfield(namespace, :campaign)
case = campaign_case(campaign, task)
params = legacy_parameters(case.model, case.representation, case.solver)
checkpoint_maxdim = h5open(checkpoint_path, "r") do file
    Int(read(file, "itensors_maxdim"))
end
if checkpoint_maxdim != params.itensors_maxdim
    params = ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
        itensors_maxdim=checkpoint_maxdim, density=params.density,
        purification_steps=params.purification_steps, scf_mixing=params.scf_mixing,
        scf_tol=params.scf_tol, scf_max_iterations=params.scf_max_iterations,
    )
end
system = read_mpo_checkpoint(checkpoint_path, params).system
levels = length(system.sites)
iszero(levels % 2) || error("square QTT requires an even number of levels")
side = 1 << (levels ÷ 2)
side == 1024 || error("this figure exporter is configured for a 1024x1024 state")
density = density_diagonal_mps(system)

# Contract the QTT MPS once. Julia stores the first QTT site as the fastest
# array index, while the physical QTT index regards it as the most-significant
# bit, hence the bit reversal when recovering a site value below.
dense_timing = @timed vec(Array(reduce(*, density), system.sites...))
dense_qtt_order = dense_timing.value
full = Matrix{Float64}(undef, side, side)
full_timing = @timed Threads.@threads for linear in 0:(side^2 - 1)
    x, y = square_lattice_decoder(linear, levels)
    reversed = Int(bitreverse(UInt(linear)) >> (8sizeof(UInt) - levels))
    full[x + 1, y + 1] = real(dense_qtt_order[reversed + 1])
end

zoom_side = 128
zoom_start = (side - zoom_side) ÷ 2
zoom = Matrix{Float64}(undef, zoom_side, zoom_side)
zoom_timing = @timed Threads.@threads for linear in 0:(zoom_side^2 - 1)
    dx, dy = square_lattice_decoder(linear, 14)
    x, y = zoom_start + dx, zoom_start + dy
    zoom[dx + 1, dy + 1] = full[x + 1, y + 1]
end

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
try
    open(joinpath(temporary, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "created_at" => string(now()), "diagnostic" => "qtt_charge_figure_export",
            "campaign" => campaign_path, "task" => task, "source_checkpoint" => checkpoint_path,
            "side" => side, "full_panel_side" => side, "full_panel_block_side" => 1,
            "zoom_start_x" => zoom_start, "zoom_start_y" => zoom_start, "zoom_side" => zoom_side,
            "density_qtt_max_chi" => maxlinkdim(density),
            "dense_qtt_contraction_time_s" => dense_timing.time,
            "full_reordering_time_s" => full_timing.time, "zoom_copy_time_s" => zoom_timing.time,
            "threads" => Threads.nthreads(),
        ))
    end
    open(joinpath(temporary, "full_raw_density.csv"), "w") do io
        println(io, "x,y,density")
        for x in 0:(side - 1), y in 0:(side - 1)
            println(io, "$(x),$(y),$(full[x + 1, y + 1])")
        end
    end
    open(joinpath(temporary, "center_128_density.csv"), "w") do io
        println(io, "x,y,density")
        for dx in 0:(zoom_side - 1), dy in 0:(zoom_side - 1)
            println(io, "$(zoom_start + dx),$(zoom_start + dy),$(zoom[dx + 1, dy + 1])")
        end
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end
println("QTT figure data exported: $output")
