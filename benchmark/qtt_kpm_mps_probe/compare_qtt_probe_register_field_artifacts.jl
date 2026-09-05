#!/usr/bin/env julia

"""Compare saved spatial QTT local-field MPSs using global Hilbert-space norms."""

using Dates
using LinearAlgebra
using Serialization
using TOML
using ITensors, ITensorMPS

length(ARGS) >= 3 || error("usage: $(PROGRAM_FILE) OUTPUT REFERENCE_FIELDS.jls CANDIDATE_FIELDS.jls...")
output, reference_path = abspath(ARGS[1]), abspath(ARGS[2])
candidate_paths = abspath.(ARGS[3:end])
ispath(output) && error("refusing to overwrite existing output directory: $output")
isfile(reference_path) || error("reference artifact does not exist: $reference_path")
all(isfile, candidate_paths) || error("a candidate artifact does not exist")

reference = deserialize(reference_path)
components = (:density, :horizontal, :vertical, :hartree)
all(hasproperty(reference, component) for component in components) || error("invalid reference artifact")

function align_spatial_sites(candidate::MPS, reference::MPS)
    length(candidate) == length(reference) || error("MPS lengths differ")
    aligned = copy(candidate)
    for position in 1:length(candidate)
        aligned[position] = replaceind(aligned[position],
            siteind(candidate, position) => siteind(reference, position))
    end
    return aligned
end

function relative_difference(candidate::MPS, reference::MPS)
    # Each Slurm process reconstructs QTT sites with fresh Index ids. They have
    # identical physical meaning, but must be aligned before cross-job tensor
    # contractions. Internal link indices remain untouched.
    candidate = align_spatial_sites(candidate, reference)
    maxdim = max(maxlinkdim(candidate) + maxlinkdim(reference), 1)
    difference = +(candidate, -reference; cutoff=1e-14, maxdim=maxdim)
    reference_norm = sqrt(abs(real(inner(reference, reference))))
    difference_norm = sqrt(abs(real(inner(difference, difference))))
    return (; absolute=difference_norm, relative=difference_norm / reference_norm,
        candidate_max_chi=maxlinkdim(candidate), reference_max_chi=maxlinkdim(reference))
end

mkpath(output)
rows = NamedTuple[]
for path in candidate_paths
    candidate = deserialize(path)
    all(hasproperty(candidate, component) for component in components) || error("invalid candidate artifact: $path")
    for component in components
        comparison = relative_difference(getproperty(candidate, component), getproperty(reference, component))
        push!(rows, (; candidate=path, component=String(component), comparison...))
    end
end
open(joinpath(output, "global_mps_differences.csv"), "w") do io
    println(io, "candidate,component,absolute_l2_difference,relative_l2_difference,candidate_max_chi,reference_max_chi")
    for row in rows
        println(io, "$(row.candidate),$(row.component),$(row.absolute),$(row.relative),$(row.candidate_max_chi),$(row.reference_max_chi)")
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "reference" => reference_path,
        "candidates" => candidate_paths,
        "norm" => "exact contraction of saved QTT MPS representations; the only approximation is the rank-sum addition used to form each difference",
    ))
end
