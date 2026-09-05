#!/usr/bin/env julia

"""Export a converged square-MPO charge density to an analysis HDF5 file.

This is a post-processing tool: it does not restart SCF, recompute the energy,
or construct the dense density matrix.  The diagonal density is first formed
as a QTT/MPS and is then either materialized as a `side × side` Float64 array
or sampled directly from the QTT.  The latter mode is intended for lattices
whose full raster is too large to store or render.

Usage:
    julia --project=. export_charge_hdf5.jl CAMPAIGN.jl TASK RESULT OUTPUT.h5 [mode] [max_dense_sites] [sample_stride]

`mode` is `auto` (default), `dense`, or `sampled`.  In `auto` mode a dense
raster is written only when `side^2 <= max_dense_sites` (default 16,777,216).
The HDF5 file is never overwritten.
"""

function usage()
    println("usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT.h5 [auto|dense|sampled] [max_dense_sites] [sample_stride]")
end

if !isempty(ARGS) && ARGS[1] in ("-h", "--help")
    usage()
    exit(0)
end

using Dates
using HDF5
using ITensors, ITensorMPS
using MPO_MeanField

length(ARGS) >= 4 || (usage(); error("missing arguments"))
campaign_path = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
result_directory = abspath(ARGS[3])
output = abspath(ARGS[4])
mode_requested = length(ARGS) >= 5 ? lowercase(ARGS[5]) : "auto"
mode_requested in ("auto", "dense", "sampled") || error("mode must be auto, dense, or sampled")
max_dense_sites = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 16_777_216
sample_stride_requested = length(ARGS) >= 7 ? parse(Int, ARGS[7]) : 64
max_dense_sites > 0 || error("max_dense_sites must be positive")
sample_stride_requested > 0 || error("sample_stride must be positive")
isfile(campaign_path) || error("campaign does not exist: $campaign_path")
ispath(output) && error("refusing to overwrite output: $output")

checkpoint_path = let candidates = (
    joinpath(result_directory, "converged_state.h5"),
    joinpath(result_directory, "final_fixed_sp2_state.h5"),
)
    found = findfirst(isfile, candidates)
    isnothing(found) && error("no converged_state.h5 or final_fixed_sp2_state.h5 in $result_directory")
    candidates[found]
end

namespace = Module(gensym(:ChargeExportCampaign))
Core.eval(namespace, :(using MPO_MeanField))
Base.include(namespace, campaign_path)
isdefined(namespace, :campaign) || error("campaign source must define campaign")
campaign = Base.invokelatest(getfield, namespace, :campaign)

case, params, label = if hasproperty(campaign, :runs)
    legacy_case = campaign.runs[task_index]
    (legacy_case, legacy_case.params, string(legacy_case.label))
else
    public_case = Base.invokelatest(campaign_case, campaign, task_index)
    (public_case, Base.invokelatest(
        legacy_parameters, public_case.model, public_case.representation, public_case.solver,
    ), string(public_case.label))
end
params isa ParametersSquare || error("charge HDF5 export currently requires a square model")

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

loaded = read_mpo_checkpoint(checkpoint_path, params)
system = loaded.system
levels = length(system.sites)
iseven(levels) || error("square QTT requires an even number of QTT sites, got $levels")
side_bits = levels ÷ 2
side = 1 << side_bits
site_count = side * side
mode = mode_requested == "auto" ? (site_count <= max_dense_sites ? "dense" : "sampled") : mode_requested
mode == "dense" && site_count > max_dense_sites && error(
    "dense export requested for $site_count sites, above max_dense_sites=$max_dense_sites",
)

println("Extracting $side×$side charge density from $checkpoint_path")
flush(stdout)
density_timing = @timed density_diagonal_mps(system)
density_mps = density_timing.value
sites = collect(siteinds(density_mps))

# The QTT index used by qtt_mps_amplitude is the logical interleaved index.
# Array(reduce(*,...)) uses the opposite physical storage order, hence the
# explicit bit reversal in the dense path below.
function logical_density(dense_qtt, linear::Int, levels::Int)
    reversed = Int(bitreverse(UInt(linear)) >>> (8sizeof(UInt) - levels))
    return real(dense_qtt[reversed + 1])
end

tmp = output * ".tmp.$(getpid())"
mkpath(dirname(output))
try
    h5open(tmp, "w") do file
        attrs(file)["format_version"] = 1
        attrs(file)["diagnostic"] = "qtt_charge_density_export"
        attrs(file)["created_at"] = string(now())
        attrs(file)["campaign"] = campaign_path
        attrs(file)["task_index"] = task_index
        attrs(file)["label"] = label
        attrs(file)["source_result"] = result_directory
        attrs(file)["source_checkpoint"] = checkpoint_path
        attrs(file)["state_kind"] = basename(checkpoint_path) == "final_fixed_sp2_state.h5" ?
            "fixed_initial_field_sp2" : "self_consistent_hf"
        attrs(file)["side"] = side
        attrs(file)["matrix_dimension"] = site_count
        attrs(file)["qtt_levels"] = levels
        attrs(file)["qtt_density_max_chi"] = maxlinkdim(density_mps)
        attrs(file)["requested_mode"] = mode_requested
        attrs(file)["export_mode"] = mode
        attrs(file)["max_dense_sites"] = max_dense_sites
        attrs(file)["density_extraction_time_s"] = density_timing.time
        attrs(file)["density_definition"] = "real(rho[i,i]) in public square/QTT ordering"
        attrs(file)["coordinate_definition"] = "x,y = square_lattice_decoder(linear, qtt_levels)"

        if mode == "dense"
            dense_timing = @timed vec(Array(reduce(*, density_mps), sites...))
            dense_qtt = dense_timing.value
            density = Matrix{Float64}(undef, side, side)
            for linear in 0:(site_count - 1)
                x, y = square_lattice_decoder(linear, levels)
                density[x + 1, y + 1] = logical_density(dense_qtt, linear, levels)
            end
            write(file, "density", density)
            attrs(file)["dense_contraction_time_s"] = dense_timing.time
            attrs(file)["density_min"] = minimum(density)
            attrs(file)["density_max"] = maximum(density)
            attrs(file)["density_mean"] = sum(density) / site_count
            attrs(file)["density_trace"] = sum(density)
        else
            # Keep the sampled output bounded for very large QTTs.  A regular
            # coordinate grid is preferable to random points for later plots.
            stride = sample_stride_requested
            nx = cld(side, stride)
            while nx * nx > 4_000_000
                stride *= 2
                nx = cld(side, stride)
            end
            xs = collect(0:stride:(side - 1))
            ys = collect(0:stride:(side - 1))
            sample_x = Vector{Int32}(undef, length(xs) * length(ys))
            sample_y = similar(sample_x)
            sample_density = Vector{Float64}(undef, length(sample_x))
            cursor = 0
            for x in xs, y in ys
                cursor += 1
                sample_x[cursor] = x
                sample_y[cursor] = y
                linear = square_lattice_index(x, y, levels) - 1
                sample_density[cursor] = real(qtt_mps_amplitude(density_mps, sites, linear))
            end
            write(file, "sample_x", sample_x)
            write(file, "sample_y", sample_y)
            write(file, "sample_density", sample_density)
            attrs(file)["sample_stride"] = stride
            attrs(file)["sample_count"] = length(sample_density)
            attrs(file)["sample_definition"] = "regular grid; values evaluated directly from density QTT"
        end
    end
    mv(tmp, output; force=false)
catch
    ispath(tmp) && rm(tmp; recursive=true, force=true)
    rethrow()
end

println("Charge density HDF5 exported: $output (mode=$mode, side=$side)")
