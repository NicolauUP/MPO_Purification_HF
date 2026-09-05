"""SCF validation using the reusable GPU block-SpAMM SP2 kernel.

Usage: julia run_block_spamm_scf_1d.jl OUTPUT [N] [V2] [U] [SEED]
"""

using CUDA
using LinearAlgebra
using SparseArrays
using TOML

include(joinpath(@__DIR__, "..", "..", "src", "solvers", "spamm", "block_sp2.jl"))
using .BlockSp2Engine

length(ARGS) in 1:5 || error("usage: run_block_spamm_scf_1d.jl OUTPUT [N] [V2] [U] [SEED]")
const OUT = abspath(ARGS[1])
const NS = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1024
const MODULATION = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 2.0
const INTERACTION = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.3
const SEED_FIELD = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 0.5
const BLOCK_SIZE = 16
const PARTICLES = NS ÷ 2
const SPAMM_TAU = 1e-6
const OUTPUT_CUTOFF = 1e-10
const SP2_STEPS = 24
const MIXING = 0.5
const DENSITY_REL_TOL = 1e-3
const STABLE_REQUIRED = 3
const MAX_SCF = 50
const TAU_AA = sqrt(2.0) - 5.0 / 6.0

isdir(OUT) && error("refusing to overwrite existing output directory: $OUT")
iseven(NS) || error("N must be even")
NS <= 2048 || error("this dense-validation driver is restricted to N <= 2048")
mkpath(OUT)

function effective_hamiltonian(n, bonds)
    diagonal = zeros(Float64, NS)
    for i in 1:NS
        i > 1 && (diagonal[i] += INTERACTION * n[i - 1])
        i < NS && (diagonal[i] += INTERACTION * n[i + 1])
    end
    hopping = [
        -1 - MODULATION * cos(2π * TAU_AA * (i + 0.5)) - INTERACTION * bonds[i]
        for i in 1:(NS - 1)
    ]
    return spdiagm(-1 => hopping, 0 => diagonal, 1 => hopping)
end

function initial_projector_argument(H)
    radius = maximum(abs.(diag(H))) + maximum([
        sum(abs, H[:, i]) - abs(H[i, i]) for i in 1:NS
    ]) + 0.25
    return sparse((radius * I - H) / (2radius)), radius
end

function fields_from_projector(P)
    n = Vector{Float64}(diag(P))
    bonds = [P[i, i + 1] for i in 1:(NS - 1)]
    return n, bonds
end

function solve_scf(method::Symbol)
    density = [iseven(i) ? 0.5 + SEED_FIELD : 0.5 - SEED_FIELD for i in 1:NS]
    bonds = zeros(Float64, NS - 1)
    stable = 0
    history = NamedTuple[]
    for iteration in 1:MAX_SCF
        H = effective_hamiltonian(density, bonds)
        if method == :dense
            eig = eigen(Symmetric(Matrix(H)))
            occupied = view(eig.vectors, :, 1:PARTICLES)
            P = occupied * occupied'
        else
            P0, _ = initial_projector_argument(H)
            step_output = joinpath(OUT, "sp2_iteration_$(lpad(iteration, 3, '0'))")
            mkpath(step_output)
            device = purify_gpu_block_sp2!(
                P0, BLOCK_SIZE, PARTICLES, SPAMM_TAU, SP2_STEPS,
                OUTPUT_CUTOFF; output=step_output,
            )
            P = Matrix(cublock_to_host_csr(device, NS))
            device = nothing
        end
        candidate_density, candidate_bonds = fields_from_projector(P)
        relative_density_change = maximum(abs, candidate_density .- density) /
                                  max(maximum(abs, candidate_density), eps())
        stable = relative_density_change <= DENSITY_REL_TOL ? stable + 1 : 0
        push!(history, (; iteration, relative_density_change,
                         particle_number=sum(candidate_density), stable))
        density .= MIXING .* candidate_density .+ (1 - MIXING) .* density
        bonds .= MIXING .* candidate_bonds .+ (1 - MIXING) .* bonds
        if method == :spamm
            GC.gc(true)
            CUDA.reclaim()
        end
        stable >= STABLE_REQUIRED && return (; converged=true, iteration, density,
                                               bonds, history)
    end
    return (; converged=false, iteration=MAX_SCF, density, bonds, history)
end

spamm = solve_scf(:spamm)
dense = solve_scf(:dense)
difference = spamm.density - dense.density
maximum_error = maximum(abs, difference)
rms_error = sqrt(sum(abs2, difference) / NS)
maximum_error <= DENSITY_REL_TOL || error(
    "block-SpAMM SCF density error exceeds 0.1% absolute occupancy: $maximum_error",
)

open(joinpath(OUT, "history.csv"), "w") do io
    println(io, "iteration,relative_density_change,particle_number,stable_iterations")
    for row in spamm.history
        println(io, "$(row.iteration),$(row.relative_density_change),$(row.particle_number),$(row.stable)")
    end
end
open(joinpath(OUT, "density.csv"), "w") do io
    println(io, "site,spamm_density,dense_density,difference")
    for i in 1:NS
        println(io, "$i,$(spamm.density[i]),$(dense.density[i]),$(difference[i])")
    end
end
open(joinpath(OUT, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "scf_converged" => spamm.converged,
        "scf_iterations" => spamm.iteration,
        "dense_scf_converged" => dense.converged,
        "dense_scf_iterations" => dense.iteration,
        "density_max_abs_error" => maximum_error,
        "density_rms_error" => rms_error,
        "particle_number" => sum(spamm.density),
    ))
end
open(joinpath(OUT, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "method" => "sparse_sp2_block_spamm", "backend" => "cuda",
        "matrix_dimension" => NS, "V2" => MODULATION, "U" => INTERACTION,
        "seed_amplitude" => SEED_FIELD, "filling" => 0.5,
        "block_size" => BLOCK_SIZE, "spamm_tau" => SPAMM_TAU,
        "output_cutoff" => OUTPUT_CUTOFF, "sp2_steps" => SP2_STEPS,
        "scf_mixing" => MIXING, "density_relative_tolerance" => DENSITY_REL_TOL,
        "stable_iterations_required" => STABLE_REQUIRED,
        "cuda_reclaim_boundary" => "after each complete SCF iteration",
    ))
end
println("block-SpAMM SCF validation complete: $OUT")
