#!/usr/bin/env julia

"""Execute one exact-diagonalization reference for an explicit 1D HF campaign.

This is an independent real, open-chain nearest-neighbour Hartree--Fock
solver. It diagonalizes the tridiagonal one-particle effective Hamiltonian at
every SCF iteration and occupies its lowest `Ne` eigenvectors. It is *not*
many-body exact diagonalization; it is an exact solver for the same 1D
mean-field functional implemented by the MPO path.

Usage:
    julia --project=. runs/common/execute_1d_dense_hf_campaign.jl CAMPAIGN_FILE TASK_INDEX
"""

using Dates
using LinearAlgebra
using Printf
using SHA
using TOML
using MPO_MeanField

length(ARGS) == 2 || error("usage: execute_1d_dense_hf_campaign.jl CAMPAIGN_FILE TASK_INDEX")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
isnothing(task_index) && error("TASK_INDEX must be an integer, got $(ARGS[2])")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
hasproperty(campaign, :name) && hasproperty(campaign, :runs) || error(
    "campaign must provide `name` and `runs`",
)
1 <= task_index <= length(campaign.runs) || error(
    "TASK_INDEX=$task_index is outside 1:$(length(campaign.runs))",
)
spec = campaign.runs[task_index]
for field in (:label, :params)
    hasproperty(spec, field) || error("campaign run $task_index is missing `$field`")
end
params = spec.params
params isa Parameters1D || error("dense HF reference requires Parameters1D")
params.t isa Number || params.t isa Function || error("1D hopping must be numeric or functional")
isnothing(params.W) || params.W isa Function || error("1D W must be nothing or a function")
isnothing(params.S) || params.S isa Function || error("1D S must be nothing or a function")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

function safe_component(text)
    component = replace(string(text), r"[^A-Za-z0-9._-]+" => "_")
    isempty(component) && error("run label must contain at least one safe character")
    return component
end

function git_revision(repo_root)
    try
        return readchomp(`git -C $repo_root rev-parse HEAD`)
    catch
        return "unavailable"
    end
end

sha1_file(path) = bytes2hex(sha1(read(path)))
relative_change(candidate, reference) = norm(candidate - reference) / max(norm(reference), sqrt(eps(Float64)))

function _one_body_terms(params::Parameters1D)
    N = 2^params.L
    onsite = isnothing(params.W) ? zeros(Float64, N) : [Float64(params.W(site - 1)) for site in 1:N]
    hopping = [Float64(params.t isa Number ? params.t : params.t(bond - 1)) for bond in 1:(N - 1)]
    seed = isnothing(params.S) ? zeros(Float64, N) : [Float64(params.S(site - 1)) for site in 1:N]
    return onsite, hopping, seed
end

function _mean_fields(rho::Matrix{Float64}, U::Float64)
    N = size(rho, 1)
    hartree = zeros(Float64, N)
    for site in 1:N
        site > 1 && (hartree[site] += U * rho[site - 1, site - 1])
        site < N && (hartree[site] += U * rho[site + 1, site + 1])
    end
    fock = [-U * rho[site, site + 1] for site in 1:(N - 1)]
    return hartree, fock
end

function _effective_hamiltonian(onsite, hopping, hartree, fock)
    SymTridiagonal(onsite .+ hartree, hopping .+ fock)
end

function _occupied_projector(H::SymTridiagonal{Float64,Vector{Float64}}, Ne::Int)
    eigenpairs = eigen(H)
    occupied = @view eigenpairs.vectors[:, 1:Ne]
    return Matrix(occupied * occupied')
end

function _energy(onsite, hopping, rho, U)
    N = size(rho, 1)
    density = diag(rho)
    bonds = [rho[site, site + 1] for site in 1:(N - 1)]
    kinetic = sum(onsite .* density) + 2sum(hopping .* bonds)
    hartree = U * sum(density[site] * density[site + 1] for site in 1:(N - 1))
    fock = -U * sum(abs2(value) for value in bonds)
    return (kinetic=kinetic, hartree=hartree, fock=fock,
        interaction=hartree + fock, total=kinetic + hartree + fock)
end

function dense_hf_scf(params::Parameters1D; io::IO=stdout)
    N = 2^params.L
    Ne = round(Int, params.density * N)
    onsite, hopping, hartree = _one_body_terms(params)
    fock = zeros(Float64, N - 1)
    rho_previous = zeros(Float64, N, N)
    rho_two_steps_ago = nothing
    history = NamedTuple[]
    residual_tolerance = params.scf_tol / 100
    converged = false
    termination_reason = :max_iterations

    println(io, "Dense 1D HF: N=$N Ne=$Ne max_iterations=$(params.scf_max_iterations) mixing=$(params.scf_mixing)")
    for iteration in 1:params.scf_max_iterations
        H = _effective_hamiltonian(onsite, hopping, hartree, fock)
        rho = _occupied_projector(H, Ne)
        new_hartree, new_fock = _mean_fields(rho, Float64(params.U))
        commutator = relative_change(Matrix(H) * rho, rho * Matrix(H))
        vh_residual = iteration == 1 ? Inf : relative_change(new_hartree, hartree)
        vf_residual = iteration == 1 ? Inf : relative_change(new_fock, fock)
        rho_residual = iteration == 1 ? Inf : relative_change(rho, rho_previous)
        two_cycle_residual = iteration < 3 ? Inf : relative_change(rho, rho_two_steps_ago)
        energy = _energy(onsite, hopping, rho, Float64(params.U))
        push!(history, (
            iteration=iteration,
            trace=sum(diag(rho)),
            vh_residual=vh_residual,
            vf_residual=vf_residual,
            rho_residual=rho_residual,
            commutator_residual=commutator,
            two_cycle_residual=two_cycle_residual,
            energy_total=energy.total,
        ))
        @printf(io, "Dense SCF %d/%d | Tr=%.12g | VH=%.3e | VF=%.3e | rho=%.3e | [H,rho]=%.3e\n",
            iteration, params.scf_max_iterations, sum(diag(rho)), vh_residual, vf_residual,
            rho_residual, commutator)

        residuals = (vh_residual, vf_residual, rho_residual, commutator)
        stable = iteration >= 2 && all(isfinite, residuals) && maximum(residuals) <= residual_tolerance
        previous_stable = length(history) >= 2 && begin
            previous = history[end - 1]
            all(isfinite, (previous.vh_residual, previous.vf_residual, previous.rho_residual, previous.commutator_residual)) &&
                maximum((previous.vh_residual, previous.vf_residual, previous.rho_residual, previous.commutator_residual)) <= residual_tolerance
        end
        if stable && previous_stable
            converged = true
            termination_reason = :converged
            println(io, "Dense SCF converged in $iteration iterations.")
            return (; converged, termination_reason, history, rho, onsite, hopping, hartree, fock)
        end
        if iteration >= 3 && isfinite(two_cycle_residual) && two_cycle_residual <= residual_tolerance && rho_residual > residual_tolerance
            termination_reason = :two_cycle_detected
            println(io, "Dense SCF stopped: detected a two-cycle.")
            return (; converged, termination_reason, history, rho, onsite, hopping, hartree, fock)
        end

        rho_two_steps_ago = rho_previous
        if iteration == 1
            hartree = new_hartree
            fock = new_fock
        else
            hartree = params.scf_mixing .* new_hartree .+ (1 - params.scf_mixing) .* hartree
            fock = params.scf_mixing .* new_fock .+ (1 - params.scf_mixing) .* fock
        end
        rho_previous = rho
    end
    return (; converged, termination_reason, history, rho=rho_previous, onsite, hopping, hartree, fock)
end

results_root = get(ENV, "MPO_RESULTS_ROOT", "")
isempty(results_root) && error("set MPO_RESULTS_ROOT to an external result directory")
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
dense_campaign_name = safe_component(campaign.name) * "_dense_hf"
run_name = @sprintf("task_%04d_%s", task_index, safe_component(spec.label))
run_dir = joinpath(results_root, dense_campaign_name, run_name)
ispath(run_dir) && error("refusing to overwrite existing result directory: $run_dir")
mkpath(run_dir)

cp(campaign_file, joinpath(run_dir, "input.jl"))
open(joinpath(run_dir, "selection.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => dense_campaign_name,
        "source_campaign" => string(campaign.name),
        "task_index" => task_index,
        "label" => string(spec.label),
        "solver" => "dense_tridiagonal_hf",
    ))
end

metadata = Dict(
    "campaign" => dense_campaign_name,
    "source_campaign" => string(campaign.name),
    "task_index" => task_index,
    "label" => string(spec.label),
    "solver" => "dense_tridiagonal_hf",
    "started_at" => string(now(UTC)),
    "julia_version" => string(VERSION),
    "git_revision" => git_revision(repo_root),
    "project_sha1" => sha1_file(joinpath(repo_root, "Project.toml")),
    "manifest_sha1" => sha1_file(joinpath(repo_root, "Manifest.toml")),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", string(task_index)),
    "threads" => Threads.nthreads(),
)

result = dense_hf_scf(params)
final_energy = _energy(result.onsite, result.hopping, result.rho, Float64(params.U))
final_H = _effective_hamiltonian(result.onsite, result.hopping, result.hartree, result.fock)
idempotency = norm(result.rho * result.rho - result.rho) / norm(result.rho)
hermiticity = norm(result.rho - result.rho') / norm(result.rho)
stationarity = relative_change(Matrix(final_H) * result.rho, result.rho * Matrix(final_H))

open(joinpath(run_dir, "scf_history.csv"), "w") do io
    write_csv_row(io, ("iteration", "trace", "vh_residual", "vf_residual", "rho_residual", "commutator_residual", "two_cycle_residual", "energy_total"))
    for record in result.history
        write_csv_row(io, (record.iteration, record.trace, record.vh_residual, record.vf_residual,
            record.rho_residual, record.commutator_residual, record.two_cycle_residual, record.energy_total))
    end
end

open(joinpath(run_dir, "site_density.csv"), "w") do io
    write_csv_row(io, ("site", "density"))
    for site in axes(result.rho, 1)
        write_csv_row(io, (site, result.rho[site, site]))
    end
end

open(joinpath(run_dir, "bond_order.csv"), "w") do io
    write_csv_row(io, ("site_left", "site_right", "real", "imag"))
    for site in 1:(size(result.rho, 1) - 1)
        write_csv_row(io, (site, site + 1, result.rho[site, site + 1], 0.0))
    end
end

open(joinpath(run_dir, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "particle_number" => sum(diag(result.rho)),
        "energy_kinetic" => final_energy.kinetic,
        "energy_hartree" => final_energy.hartree,
        "energy_fock" => final_energy.fock,
        "energy_interaction" => final_energy.interaction,
        "energy_total" => final_energy.total,
        "hermiticity_residual" => hermiticity,
        "idempotency_residual" => idempotency,
        "stationarity_residual" => stationarity,
    ))
end

metadata["finished_at"] = string(now(UTC))
metadata["scf_converged"] = result.converged
metadata["scf_termination_reason"] = string(result.termination_reason)
metadata["scf_iterations"] = length(result.history)
metadata["matrix_dimension"] = size(result.rho, 1)
open(joinpath(run_dir, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

println("Result directory: $run_dir")
println("Dense SCF: converged=$(result.converged) termination=$(result.termination_reason)")
exit(result.converged ? 0 : 2)
