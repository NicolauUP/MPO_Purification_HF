#!/usr/bin/env julia

"""Create an exact dense reference for a square campaign's initial H0 + S.

This deliberately stops before purification and SCF. It gives the exact
spectrum and occupied projector of the Hamiltonian that starts SCF iteration
one, allowing SP2 scaling and truncation to be diagnosed independently of
Hartree/Fock extraction.
"""

using Dates
using LinearAlgebra
using Printf
using SHA
using TOML
using MPO_MeanField

length(ARGS) == 2 || error("usage: initial_square_hamiltonian_reference.jl CAMPAIGN_FILE TASK_INDEX")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
isnothing(task_index) && error("TASK_INDEX must be an integer")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
params = spec.params
params isa ParametersSquare || error("initial square reference requires ParametersSquare")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
safe_component(value) = begin
    result = replace(string(value), r"[^A-Za-z0-9._-]+" => "_")
    isempty(result) && error("unsafe empty path component")
    result
end
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end

function direct_initial_hamiltonian(params::ParametersSquare)
    L = params.L
    side = 2^div(L, 2)
    N = side^2
    H = zeros(Float64, N, N)
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) : Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) : Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, L)
        H[site, site] = (isnothing(params.W) ? 0.0 : Float64(params.W(x, y))) +
                        (isnothing(params.S) ? 0.0 : Float64(params.S(x, y)))
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, L)
            H[site, neighbour] = H[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, L)
            H[site, neighbour] = H[neighbour, site] = ty(x, y)
        end
    end
    return H
end

results_root = get(ENV, "MPO_RESULTS_ROOT", "")
isempty(results_root) && error("set MPO_RESULTS_ROOT to an external result directory")
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
reference_campaign = safe_component(campaign.name) * "_initial_reference"
run_name = @sprintf("task_%04d_%s", task_index, safe_component(spec.label))
run_dir = joinpath(results_root, reference_campaign, run_name)
ispath(run_dir) && error("refusing to overwrite existing result directory: $run_dir")
mkpath(run_dir)

H_initial = direct_initial_hamiltonian(params)
N = size(H_initial, 1)
Ne = round(Int, params.density * N)
eigenpairs = eigen(Symmetric(H_initial))
eigenvalues = eigenpairs.values
occupied_vectors = @view eigenpairs.vectors[:, 1:Ne]
rho = Matrix(occupied_vectors * occupied_vectors')
gap = eigenvalues[Ne + 1] - eigenvalues[Ne]

cp(campaign_file, joinpath(run_dir, "input.jl"))
open(joinpath(run_dir, "selection.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => reference_campaign,
        "source_campaign" => string(campaign.name),
        "task_index" => task_index,
        "label" => string(spec.label),
        "solver" => "dense_initial_square_hamiltonian",
    ))
end
open(joinpath(run_dir, "spectrum.csv"), "w") do io
    write_csv_row(io, ("state", "eigenvalue", "occupation"))
    for state in eachindex(eigenvalues)
        write_csv_row(io, (state, eigenvalues[state], state <= Ne ? 1 : 0))
    end
end
open(joinpath(run_dir, "site_density.csv"), "w") do io
    write_csv_row(io, ("site", "density"))
    for site in 1:N
        write_csv_row(io, (site, rho[site, site]))
    end
end
open(joinpath(run_dir, "bond_order.csv"), "w") do io
    write_csv_row(io, ("site_left", "site_right", "orientation", "real", "imag"))
    for (site, neighbour, orientation) in square_undirected_bonds(params.L)
        write_csv_row(io, (site, neighbour, orientation, rho[site, neighbour], 0.0))
    end
end
open(joinpath(run_dir, "initial_reference.toml"), "w") do io
    TOML.print(io, Dict(
        "matrix_dimension" => N,
        "target_particles" => Ne,
        "spectral_lower" => first(eigenvalues),
        "spectral_upper" => last(eigenvalues),
        "homo_energy" => eigenvalues[Ne],
        "lumo_energy" => eigenvalues[Ne + 1],
        "fermi_gap" => gap,
        "projector_trace" => tr(rho),
        "projector_idempotency_residual" => norm(rho * rho - rho) / norm(rho),
        "hamiltonian_hermiticity_residual" => norm(H_initial - H_initial') / max(norm(H_initial), sqrt(eps(Float64))),
        "supplied_spectral_lower" => Float64(spec.spectral_bounds[1]),
        "supplied_spectral_upper" => Float64(spec.spectral_bounds[2]),
        "lower_padding" => first(eigenvalues) - Float64(spec.spectral_bounds[1]),
        "upper_padding" => Float64(spec.spectral_bounds[2]) - last(eigenvalues),
    ))
end
open(joinpath(run_dir, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => reference_campaign,
        "source_campaign" => string(campaign.name),
        "task_index" => task_index,
        "label" => string(spec.label),
        "solver" => "dense_initial_square_hamiltonian",
        "started_at" => string(now(UTC)),
        "finished_at" => string(now(UTC)),
        "julia_version" => string(VERSION),
        "git_revision" => git_revision(repo_root),
        "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
        "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
        "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", string(task_index)),
        "threads" => Threads.nthreads(),
    ))
end

println("Initial dense reference: N=$N Ne=$Ne")
println("Exact spectrum: [$(first(eigenvalues)), $(last(eigenvalues))]")
println("HOMO=$(eigenvalues[Ne]) LUMO=$(eigenvalues[Ne + 1]) gap=$gap")
println("Result directory: $run_dir")
