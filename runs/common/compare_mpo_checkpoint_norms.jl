#!/usr/bin/env julia

"""Compare two converged MPO density matrices without real-space expansion.

Usage:
    julia --project=. compare_mpo_checkpoint_norms.jl CHI256.h5 CHI512.h5 OUTPUT.toml

The relative error is ||rho_512-rho_256||_F / ||rho_512||_F, evaluated from
MPO inner products. No site-resolved observable is evaluated.
"""

using HDF5
using ITensors
using ITensorMPS
using TOML

length(ARGS) == 3 || error("usage: compare_mpo_checkpoint_norms.jl CHI256.h5 CHI512.h5 OUTPUT.toml")
path256, path512, output = abspath.(ARGS)
isfile(path256) || error("missing chi=256 checkpoint: $path256")
isfile(path512) || error("missing chi=512 checkpoint: $path512")
ispath(output) && error("refusing to overwrite output: $output")

function read_rho(path)
    h5open(path, "r") do file
        return (
            rho=read(file, "rho", MPO),
            target_particles=2^(Int(read(file, "qtt_levels")) - 1),
            maxdim=Int(read(file, "itensors_maxdim")),
        )
    end
end

"""Relabel candidate physical legs to the corresponding reference legs.

Independent QTT builds use distinct ITensor index identities even when their
binary ordering and physical dimensions agree. Link indices remain untouched.
"""
function align_physical_indices!(candidate::MPO, reference::MPO)
    for site in eachindex(candidate)
        candidate_legs = collect(siteinds(candidate, site))
        reference_legs = collect(siteinds(reference, site))
        length(candidate_legs) == length(reference_legs) || error(
            "site $site has incompatible physical-leg counts",
        )
        replacements = Pair{Index,Index}[]
        for leg in candidate_legs
            matches = filter(reference_legs) do reference_leg
                dim(reference_leg) == dim(leg) && plev(reference_leg) == plev(leg)
            end
            length(matches) == 1 || error(
                "site $site has no unique matching physical leg for $leg",
            )
            push!(replacements, leg => only(matches))
        end
        candidate[site] = replaceinds(candidate[site], replacements)
    end
    return candidate
end

state256 = read_rho(path256)
state512 = read_rho(path512)
length(state256.rho) == length(state512.rho) || error("MPO lengths differ")
state256.target_particles == state512.target_particles || error("particle targets differ")
align_physical_indices!(state256.rho, state512.rho)

n256_sq = max(0.0, real(inner(state256.rho, state256.rho)))
n512_sq = max(0.0, real(inner(state512.rho, state512.rho)))
cross = real(inner(state512.rho, state256.rho))
difference_sq = max(0.0, n512_sq + n256_sq - 2 * cross)
difference_norm = sqrt(difference_sq)
relative_error = difference_norm / max(sqrt(n512_sq), sqrt(eps(Float64)))

mkpath(dirname(output))
open(output, "w") do io
    TOML.print(io, Dict(
        "comparison" => "qtt_mpo_density_frobenius_norm",
        "reference_maxdim" => state512.maxdim,
        "candidate_maxdim" => state256.maxdim,
        "qtt_levels" => length(state256.rho),
        "target_particles" => state256.target_particles,
        "rho_256_max_bond_dimension" => maxlinkdim(state256.rho),
        "rho_512_max_bond_dimension" => maxlinkdim(state512.rho),
        "rho_256_trace" => real(tr(state256.rho)),
        "rho_512_trace" => real(tr(state512.rho)),
        "rho_difference_frobenius_norm" => difference_norm,
        "rho_relative_frobenius_error" => relative_error,
    ))
end

println("relative MPO density error = $relative_error")
