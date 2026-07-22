#!/usr/bin/env julia

"""Execute an independent dense one-particle HF reference for a square campaign.

The dense path reconstructs the open-square Hamiltonian directly from the
campaign functions; it does not obtain a dense matrix from an MPO. It is exact
diagonalization of the same *mean-field* equations, not many-body ED.
"""

using Dates
using LinearAlgebra
using Printf
using SHA
using TOML
using MPO_MeanField

length(ARGS) == 2 || error("usage: execute_square_dense_hf_campaign.jl CAMPAIGN_FILE TASK_INDEX")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
isnothing(task_index) && error("TASK_INDEX must be an integer")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
for field in (:label, :params)
    hasproperty(spec, field) || error("campaign run $task_index is missing `$field`")
end
params = spec.params
params isa ParametersSquare || error("dense square HF reference requires ParametersSquare")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
safe_component(value) = begin
    result = replace(string(value), r"[^A-Za-z0-9._-]+" => "_")
    isempty(result) && error("unsafe empty path component")
    result
end
relative_change(candidate, reference) = norm(candidate - reference) / max(norm(reference), sqrt(eps(Float64)))
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end

function direct_square_h0(params::ParametersSquare)
    L = params.L
    side = 2^div(L, 2)
    N = side^2
    H0 = zeros(Float64, N, N)
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) : Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) : Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, L)
        H0[site, site] = isnothing(params.W) ? 0.0 : Float64(params.W(x, y))
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, L)
            H0[site, neighbour] = H0[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, L)
            H0[site, neighbour] = H0[neighbour, site] = ty(x, y)
        end
    end
    return H0
end

function initial_hartree_seed(params::ParametersSquare)
    N = 2^params.L
    seed = zeros(Float64, N, N)
    isnothing(params.S) && return seed
    for site in 1:N
        x, y = square_lattice_decoder(site - 1, params.L)
        seed[site, site] = Float64(params.S(x, y))
    end
    return seed
end

function mean_fields_square(rho::Matrix{Float64}, U::Float64, L::Int)
    N = size(rho, 1)
    VH = zeros(Float64, N, N)
    VF = zeros(Float64, N, N)
    for site in 1:N
        for neighbour in values(square_neighbours(site, L))
            !isnothing(neighbour) && (VH[site, site] += U * rho[neighbour, neighbour])
        end
    end
    for (site, neighbour, _) in square_undirected_bonds(L)
        VF[site, neighbour] = VF[neighbour, site] = -U * rho[site, neighbour]
    end
    return VH, VF
end

function hf_energy_square(H0, rho, U, L)
    density = diag(rho)
    hartree = 0.0
    fock = 0.0
    for (site, neighbour, _) in square_undirected_bonds(L)
        hartree += U * density[site] * density[neighbour]
        fock -= U * abs2(rho[site, neighbour])
    end
    kinetic = real(tr(H0 * rho))
    return (; kinetic, hartree, fock, interaction=hartree + fock, total=kinetic + hartree + fock)
end

function dense_square_hf(params::ParametersSquare; io::IO=stdout)
    N = 2^params.L
    Ne = round(Int, params.density * N)
    H0 = direct_square_h0(params)
    VH = initial_hartree_seed(params)
    VF = zeros(Float64, N, N)
    rho_previous = zeros(Float64, N, N)
    rho_two_steps_ago = nothing
    history = NamedTuple[]
    tolerance = params.scf_tol / 100
    println(io, "Dense square HF: N=$N Ne=$Ne max_iterations=$(params.scf_max_iterations) mixing=$(params.scf_mixing)")
    for iteration in 1:params.scf_max_iterations
        H = H0 + VH + VF
        eigenpairs = eigen(Symmetric(H))
        occupied = @view eigenpairs.vectors[:, 1:Ne]
        rho = Matrix(occupied * occupied')
        new_VH, new_VF = mean_fields_square(rho, Float64(params.U), params.L)
        vh_residual = iteration == 1 ? Inf : relative_change(new_VH, VH)
        vf_residual = iteration == 1 ? Inf : relative_change(new_VF, VF)
        rho_residual = iteration == 1 ? Inf : relative_change(rho, rho_previous)
        two_cycle_residual = iteration < 3 ? Inf : relative_change(rho, rho_two_steps_ago)
        commutator_residual = relative_change(H * rho, rho * H)
        energy = hf_energy_square(H0, rho, Float64(params.U), params.L)
        push!(history, (; iteration, trace=tr(rho), vh_residual, vf_residual, rho_residual, commutator_residual, two_cycle_residual, energy_total=energy.total))
        @printf(io, "Dense SCF %d/%d | Tr=%.12g | VH=%.3e | VF=%.3e | rho=%.3e | [H,rho]=%.3e\n", iteration, params.scf_max_iterations, tr(rho), vh_residual, vf_residual, rho_residual, commutator_residual)
        residuals = (vh_residual, vf_residual, rho_residual, commutator_residual)
        stable = iteration >= 2 && all(isfinite, residuals) && maximum(residuals) <= tolerance
        previous_stable = length(history) >= 2 && begin
            record = history[end - 1]
            previous = (record.vh_residual, record.vf_residual, record.rho_residual, record.commutator_residual)
            all(isfinite, previous) && maximum(previous) <= tolerance
        end
        if stable && previous_stable
            println(io, "Dense square SCF converged in $iteration iterations.")
            return (; converged=true, termination_reason=:converged, history, rho, H0, VH, VF, H)
        elseif iteration >= 3 && isfinite(two_cycle_residual) && two_cycle_residual <= tolerance && rho_residual > tolerance
            println(io, "Dense square SCF stopped: detected a two-cycle.")
            return (; converged=false, termination_reason=:two_cycle_detected, history, rho, H0, VH, VF, H)
        end
        rho_two_steps_ago = rho_previous
        if iteration == 1
            VH, VF = new_VH, new_VF
        else
            VH = params.scf_mixing * new_VH + (1 - params.scf_mixing) * VH
            VF = params.scf_mixing * new_VF + (1 - params.scf_mixing) * VF
        end
        rho_previous = rho
    end
    return (; converged=false, termination_reason=:max_iterations, history, rho=rho_previous, H0, VH, VF, H=H0 + VH + VF)
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
    TOML.print(io, Dict("campaign" => dense_campaign_name, "source_campaign" => string(campaign.name), "task_index" => task_index, "label" => string(spec.label), "solver" => "dense_square_hf"))
end
metadata = Dict(
    "campaign" => dense_campaign_name, "source_campaign" => string(campaign.name),
    "task_index" => task_index, "label" => string(spec.label), "solver" => "dense_square_hf",
    "started_at" => string(now(UTC)), "julia_version" => string(VERSION),
    "git_revision" => git_revision(repo_root),
    "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
    "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", string(task_index)),
    "threads" => Threads.nthreads(),
)
result = dense_square_hf(params)
energy = hf_energy_square(result.H0, result.rho, Float64(params.U), params.L)
idempotency = norm(result.rho * result.rho - result.rho) / norm(result.rho)
hermiticity = norm(result.rho - result.rho') / norm(result.rho)
stationarity = relative_change(result.H * result.rho, result.rho * result.H)
open(joinpath(run_dir, "scf_history.csv"), "w") do io
    write_csv_row(io, ("iteration", "trace", "vh_residual", "vf_residual", "rho_residual", "commutator_residual", "two_cycle_residual", "energy_total"))
    for record in result.history
        write_csv_row(io, (record.iteration, record.trace, record.vh_residual, record.vf_residual, record.rho_residual, record.commutator_residual, record.two_cycle_residual, record.energy_total))
    end
end
open(joinpath(run_dir, "site_density.csv"), "w") do io
    write_csv_row(io, ("site", "density"))
    for site in axes(result.rho, 1)
        write_csv_row(io, (site, result.rho[site, site]))
    end
end
open(joinpath(run_dir, "bond_order.csv"), "w") do io
    write_csv_row(io, ("site_left", "site_right", "orientation", "real", "imag"))
    for (site, neighbour, orientation) in square_undirected_bonds(params.L)
        write_csv_row(io, (site, neighbour, orientation, result.rho[site, neighbour], 0.0))
    end
end
open(joinpath(run_dir, "observables.toml"), "w") do io
    TOML.print(io, Dict("particle_number" => tr(result.rho), "energy_kinetic" => energy.kinetic, "energy_hartree" => energy.hartree, "energy_fock" => energy.fock, "energy_interaction" => energy.interaction, "energy_total" => energy.total, "hermiticity_residual" => hermiticity, "idempotency_residual" => idempotency, "stationarity_residual" => stationarity))
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
