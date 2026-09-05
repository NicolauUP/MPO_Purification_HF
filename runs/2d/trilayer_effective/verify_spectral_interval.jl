#!/usr/bin/env julia

"""Check exact initial-field spectra against the conservative HF interval."""

using LinearAlgebra
using MPO_MeanField
using TOML

length(ARGS) in (2, 3) || error("usage: verify_spectral_interval.jl CAMPAIGN.jl OUTPUT.toml [SIDE]")
campaign_file = abspath(ARGS[1])
output = abspath(ARGS[2])
side = length(ARGS) == 3 ? parse(Int, ARGS[3]) : 32
side >= 2 || error("SIDE must be at least 2")
isfile(campaign_file) || error("campaign does not exist: $campaign_file")

include(campaign_file)
@isdefined(campaign) || error("campaign source must define `campaign`")

function initial_matrix(model::SquareModel, side::Int)
    nx, ny = model.size
    nx == ny || error("spectral oracle requires an equal-square campaign, got $(model.size)")
    side <= nx || error("SIDE=$side exceeds campaign side $nx")
    nx = ny = side
    n = side * side
    matrix = zeros(Float64, n, n)
    index(x, y) = x + nx * y + 1
    for y in 0:(ny - 1), x in 0:(nx - 1)
        i = index(x, y)
        !isnothing(model.potential) && (matrix[i, i] += model.potential(x, y))
        if x + 1 < nx
            j = index(x + 1, y)
            value = model.hopping[1](x + 0.5, y)
            matrix[i, j] = value
            matrix[j, i] = value
        end
        if y + 1 < ny
            j = index(x, y + 1)
            value = model.hopping[2](x, y + 0.5)
            matrix[i, j] = value
            matrix[j, i] = value
        end
    end
    return matrix
end

results = Dict{String,Any}("campaign" => campaign.name, "campaign_file" => campaign_file,
    "side" => side, "oracle" => "exact dense spectrum of an open side-by-side restriction",
    "cases" => Any[])
for (task, case) in enumerate(campaign.cases)
    println("checking task $task ($(case.label))")
    flush(stdout)
    model = case.model
    model isa SquareModel || error("case $(case.label) is not a SquareModel")
    h = initial_matrix(model, side)
    println("  matrix built"); flush(stdout)
    eigenvalues = eigvals(Symmetric(h))
    println("  spectrum done"); flush(stdout)
    nx, ny = model.size
    max_tx = maximum(abs(model.hopping[1](x + 0.5, y))
        for y in 0:(side - 1), x in 0:(side - 2))
    max_ty = maximum(abs(model.hopping[2](x, y + 0.5))
        for y in 0:(side - 2), x in 0:(side - 1))
    potential_values = isnothing(model.potential) ? [0.0] :
        [model.potential(x, y) for y in 0:(side - 1), x in 0:(side - 1)]
    potential_bounds = (minimum(potential_values), maximum(potential_values))
    params = legacy_parameters(model, case.representation, case.solver)
    conservative = square_scf_spectral_bounds(params;
        potential_bounds=potential_bounds, hopping_abs_bounds=(max_tx, max_ty), margin=0.5)
    supplied = case.spectral_bounds
    covers = !isnothing(supplied) && supplied[1] <= conservative[1] && supplied[2] >= conservative[2]
    push!(results["cases"], Dict("task" => task, "label" => case.label,
        "initial_lambda_min" => first(eigenvalues), "initial_lambda_max" => last(eigenvalues),
        "max_abs_tx" => max_tx, "max_abs_ty" => max_ty,
        "potential_min" => potential_bounds[1], "potential_max" => potential_bounds[2],
        "conservative_hf_lower" => conservative[1], "conservative_hf_upper" => conservative[2],
        "supplied_lower" => isnothing(supplied) ? NaN : supplied[1],
        "supplied_upper" => isnothing(supplied) ? NaN : supplied[2],
        "supplied_covers_conservative" => covers))
end
mkpath(dirname(output))
open(output, "w") do io
    TOML.print(io, results)
end
println("spectral interval report: $output")
for case in results["cases"]
    println(case["label"], ": initial=[", case["initial_lambda_min"], ", ", case["initial_lambda_max"],
        "] conservative=[", case["conservative_hf_lower"], ", ", case["conservative_hf_upper"], "]")
end
