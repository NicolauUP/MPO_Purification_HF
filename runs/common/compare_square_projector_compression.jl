#!/usr/bin/env julia

"""Compare exact-projector compression with post-compression of an SP2 MPO.

The dense matrix conversion uses paired `(row bit, column bit)` QTT indices in
the same most-significant-bit-first ordering as `MatrixChecker`. This script
is intended for small dense-reference systems, currently the 32x32 square
lattice (`N=1024`).
"""

using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using SHA
using TOML

length(ARGS) in (6, 7, 9) || error(
    "usage: compare_square_projector_compression.jl CAMPAIGN_FILE TASK_INDEX MODE MAXDIM OUTPUT_DIRECTORY CUTOFFS [ITENSORS_TOL [SP2_IDEMPOTENCY_TOLERANCE SP2_RELATIVE_TRACE_TOLERANCE]]",
)
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
mode = Symbol(ARGS[3])
maxdim = tryparse(Int, ARGS[4])
output = abspath(ARGS[5])
cutoffs = try
    parse.(Float64, split(ARGS[6], ','))
catch
    error("CUTOFFS must be a comma-separated list of numbers")
end
itensors_tol_override = length(ARGS) == 7 ? tryparse(Float64, ARGS[7]) : nothing
if length(ARGS) == 9
    itensors_tol_override = tryparse(Float64, ARGS[7])
end
sp2_idempotency_tolerance = length(ARGS) == 9 ? tryparse(Float64, ARGS[8]) : 1e-3
sp2_relative_trace_tolerance = length(ARGS) == 9 ? tryparse(Float64, ARGS[9]) : nothing
isnothing(task_index) && error("TASK_INDEX must be an integer")
isnothing(maxdim) && error("MAXDIM must be an integer")
mode in (:exact, :sp2) || error("MODE must be `exact` or `sp2`")
maxdim > 0 || error("MAXDIM must be positive")
!isempty(cutoffs) || error("at least one cutoff is required")
all(cutoff -> isfinite(cutoff) && cutoff >= 0, cutoffs) ||
    error("all cutoffs must be finite and nonnegative")
if length(ARGS) >= 7
    isnothing(itensors_tol_override) && error("ITENSORS_TOL must be a number")
    isfinite(itensors_tol_override) && itensors_tol_override > 0 ||
        error("ITENSORS_TOL must be finite and positive")
end
isnothing(sp2_idempotency_tolerance) &&
    error("SP2_IDEMPOTENCY_TOLERANCE must be a number")
isfinite(sp2_idempotency_tolerance) && sp2_idempotency_tolerance > 0 ||
    error("SP2_IDEMPOTENCY_TOLERANCE must be finite and positive")
if length(ARGS) == 9
    isnothing(sp2_relative_trace_tolerance) &&
        error("SP2_RELATIVE_TRACE_TOLERANCE must be a number")
    isfinite(sp2_relative_trace_tolerance) && sp2_relative_trace_tolerance > 0 ||
        error("SP2_RELATIVE_TRACE_TOLERANCE must be finite and positive")
end
length(unique(cutoffs)) == length(cutoffs) || error("cutoffs must be unique")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("comparison requires ParametersSquare")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end
rms(values) = sqrt(sum(abs2, values) / length(values))

function with_numerics(
    params::ParametersSquare,
    maxdim::Int,
    itensors_tol::Float64,
)
    return ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=itensors_tol,
        itensors_maxdim=maxdim, density=params.density,
        purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end

function direct_initial_hamiltonian(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    H = zeros(Float64, N, N)
    tx(x, y) = params.t[1] isa Number ?
        Float64(params.t[1]) : Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ?
        Float64(params.t[2]) : Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, params.L)
        H[site, site] = (isnothing(params.W) ? 0.0 : Float64(params.W(x, y))) +
                        (isnothing(params.S) ? 0.0 : Float64(params.S(x, y)))
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, params.L)
            H[site, neighbour] = H[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, params.L)
            H[site, neighbour] = H[neighbour, site] = ty(x, y)
        end
    end
    return H
end

function matrix_to_interleaved_mpo(
    matrix::AbstractMatrix,
    sites::Vector{<:Index};
    cutoff::Real,
    maxdim::Int,
)
    L = length(sites)
    N = 2^L
    size(matrix) == (N, N) || throw(DimensionMismatch(
        "matrix size $(size(matrix)) is incompatible with $L QTT sites",
    ))
    values = zeros(eltype(matrix), ntuple(_ -> 2, 2L))
    coordinates = Vector{Int}(undef, 2L)
    for row in 1:N, column in 1:N
        for site in 1:L
            shift = L - site
            coordinates[2site - 1] = ((row - 1) >> shift & 1) + 1
            coordinates[2site] = ((column - 1) >> shift & 1) + 1
        end
        values[coordinates...] = matrix[row, column]
    end
    physical_indices = Index[]
    for site in sites
        push!(physical_indices, prime(site), dag(site))
    end
    return MPO(
        ITensor(values, physical_indices...),
        sites;
        cutoff=Float64(cutoff),
        maxdim=maxdim,
    )
end

function interleaved_mpo_to_matrix(mpo::MPO, sites::Vector{<:Index})
    L = length(sites)
    N = 2^L
    physical_indices = Index[]
    for site in sites
        push!(physical_indices, prime(site), dag(site))
    end
    values = Array(contract(mpo), physical_indices...)
    matrix = zeros(eltype(values), N, N)
    coordinates = Vector{Int}(undef, 2L)
    for row in 1:N, column in 1:N
        for site in 1:L
            shift = L - site
            coordinates[2site - 1] = ((row - 1) >> shift & 1) + 1
            coordinates[2site] = ((column - 1) >> shift & 1) + 1
        end
        matrix[row, column] = values[coordinates...]
    end
    return matrix
end

function checkerboard_order(density, L)
    total = sum(eachindex(density)) do site
        x, y = square_lattice_decoder(site - 1, L)
        (iseven(x + y) ? 1.0 : -1.0) * real(density[site])
    end
    return total / length(density)
end

function matrix_diagnostics(candidate, exact_projector, H, bonds, L)
    candidate_norm = max(norm(candidate), sqrt(eps(Float64)))
    exact_norm = max(norm(exact_projector), sqrt(eps(Float64)))
    Hcandidate = H * candidate
    density = diag(candidate)
    exact_density = diag(exact_projector)
    density_error = density - exact_density
    bond_error = ComplexF64[
        candidate[site, neighbour] - exact_projector[site, neighbour]
        for (site, neighbour, _) in bonds
    ]
    horizontal = findall(index -> bonds[index][3] == :horizontal, eachindex(bonds))
    vertical = findall(index -> bonds[index][3] == :vertical, eachindex(bonds))
    hermitian_candidate = Hermitian((candidate + candidate') / 2)
    occupations = eigvals(hermitian_candidate)
    energy = real(tr(Hcandidate))
    exact_energy = real(tr(H * exact_projector))
    return (
        trace=real(tr(candidate)),
        relative_projector_error=norm(candidate - exact_projector) / exact_norm,
        idempotency_residual=norm(candidate * candidate - candidate) / candidate_norm,
        hermiticity_residual=norm(candidate - candidate') / candidate_norm,
        stationarity_residual=norm(Hcandidate - candidate * H) /
                              max(norm(Hcandidate), sqrt(eps(Float64))),
        density_max_abs_error=maximum(abs, density_error),
        density_rms_error=rms(density_error),
        horizontal_bond_max_abs_error=maximum(abs, bond_error[horizontal]),
        horizontal_bond_rms_error=rms(bond_error[horizontal]),
        vertical_bond_max_abs_error=maximum(abs, bond_error[vertical]),
        vertical_bond_rms_error=rms(bond_error[vertical]),
        checkerboard_order=checkerboard_order(density, L),
        energy=energy,
        energy_error=energy - exact_energy,
        occupation_min=first(occupations),
        occupation_max=last(occupations),
        density=density,
    )
end

mean_link_dimension(mpo::MPO) = length(mpo) == 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) / (length(mpo) - 1)
link_dimensions(mpo::MPO) = length(mpo) == 1 ? Int[] :
    [dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)]

operational_itensors_tol = isnothing(itensors_tol_override) ?
    spec.params.itensors_tol : itensors_tol_override
params = with_numerics(spec.params, maxdim, operational_itensors_tol)
bounds = validate_spectral_bounds(spec.spectral_bounds...)
N = 2^params.L
Ne = round(Int, params.density * N)
sp2_absolute_trace_tolerance = isnothing(sp2_relative_trace_tolerance) ?
    nothing : sp2_relative_trace_tolerance * Ne
bonds = collect(square_undirected_bonds(params.L))
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
started_at = now(UTC)
mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))

H = direct_initial_hamiltonian(params)
dense_reference = @timed begin
    eigenpairs = eigen(Symmetric(H))
    occupied = @view eigenpairs.vectors[:, 1:Ne]
    projector = Matrix(occupied * occupied')
    (; eigenpairs, projector)
end
eigenpairs = dense_reference.value.eigenpairs
exact_projector = dense_reference.value.projector

sys = System(params)
sites = sys.sites
source_mpo = nothing
purification = nothing
if mode == :sp2
    rho0 = construct_rho_0(
        sys, params, bounds...;
        method=:sp2,
        verify_spectral_bounds=false,
    )
    purification = @timed open(joinpath(output, "sp2_progress.txt"), "w") do progress
        perform_purification(
            rho0, params;
            method=:sp2,
            verbose=1,
            io=progress,
            overwrite_progress=false,
            sp2_idempotency_tolerance=sp2_idempotency_tolerance,
            sp2_trace_tolerance=sp2_absolute_trace_tolerance,
            spectral_bounds=bounds,
            spectral_bounds_validation=:supplied_analytical,
        )
    end
    source_mpo = purification.value.rho
end

summary_header = (
    "source", "cutoff", "requested_maxdim", "max_chi", "mean_chi",
    "trace", "trace_error", "relative_projector_error", "idempotency_residual",
    "hermiticity_residual", "stationarity_residual", "density_max_abs_error",
    "density_rms_error", "horizontal_bond_max_abs_error",
    "horizontal_bond_rms_error", "vertical_bond_max_abs_error",
    "vertical_bond_rms_error", "checkerboard_order", "energy", "energy_error",
    "occupation_min", "occupation_max", "compression_time_s",
    "compression_allocations_bytes", "dense_conversion_time_s",
    "dense_conversion_allocations_bytes",
)

open(joinpath(output, "summary.csv"), "w") do summary_io
    open(joinpath(output, "bond_dimensions.csv"), "w") do bonds_io
        open(joinpath(output, "density_profiles.csv"), "w") do density_io
            write_csv_row(summary_io, summary_header)
            write_csv_row(bonds_io, ("source", "cutoff", "bond", "dimension"))
            write_csv_row(density_io, (
                "source", "cutoff", "site", "x", "y", "exact_density",
                "candidate_density", "error",
            ))

            candidate_descriptors = if mode == :exact
                [(label="exact_compressed", cutoff=cutoff) for cutoff in cutoffs]
            else
                vcat(
                    [(label="sp2_source", cutoff=NaN)],
                    [(label="sp2_recompressed", cutoff=cutoff) for cutoff in cutoffs],
                )
            end

            exact_density = diag(exact_projector)
            for descriptor in candidate_descriptors
                calculation = if descriptor.label == "exact_compressed"
                    @timed matrix_to_interleaved_mpo(
                        exact_projector, sites;
                        cutoff=descriptor.cutoff,
                        maxdim=maxdim,
                    )
                elseif descriptor.label == "sp2_source"
                    (; value=copy(source_mpo), time=0.0, bytes=0)
                else
                    @timed begin
                        candidate = copy(source_mpo)
                        truncate!(candidate; cutoff=descriptor.cutoff, maxdim=maxdim)
                        candidate
                    end
                end
                candidate = calculation.value
                dense_candidate = @timed interleaved_mpo_to_matrix(candidate, sites)
                diagnostics = matrix_diagnostics(
                    dense_candidate.value, exact_projector, H, bonds, params.L,
                )
                write_csv_row(summary_io, (
                    descriptor.label,
                    descriptor.cutoff,
                    maxdim,
                    maxlinkdim(candidate),
                    mean_link_dimension(candidate),
                    diagnostics.trace,
                    abs(diagnostics.trace - Ne),
                    diagnostics.relative_projector_error,
                    diagnostics.idempotency_residual,
                    diagnostics.hermiticity_residual,
                    diagnostics.stationarity_residual,
                    diagnostics.density_max_abs_error,
                    diagnostics.density_rms_error,
                    diagnostics.horizontal_bond_max_abs_error,
                    diagnostics.horizontal_bond_rms_error,
                    diagnostics.vertical_bond_max_abs_error,
                    diagnostics.vertical_bond_rms_error,
                    diagnostics.checkerboard_order,
                    diagnostics.energy,
                    diagnostics.energy_error,
                    diagnostics.occupation_min,
                    diagnostics.occupation_max,
                    calculation.time,
                    calculation.bytes,
                    dense_candidate.time,
                    dense_candidate.bytes,
                ))
                for (bond, dimension) in enumerate(link_dimensions(candidate))
                    write_csv_row(bonds_io, (
                        descriptor.label, descriptor.cutoff, bond, dimension,
                    ))
                end
                for site in 1:N
                    x, y = square_lattice_decoder(site - 1, params.L)
                    write_csv_row(density_io, (
                        descriptor.label,
                        descriptor.cutoff,
                        site,
                        x,
                        y,
                        exact_density[site],
                        diagnostics.density[site],
                        diagnostics.density[site] - exact_density[site],
                    ))
                end
                flush(summary_io)
                flush(bonds_io)
                flush(density_io)
            end
        end
    end
end

metadata = Dict(
    "started_at" => string(started_at),
    "finished_at" => string(now(UTC)),
    "campaign" => string(campaign.name),
    "label" => string(spec.label),
    "task_index" => task_index,
    "mode" => string(mode),
    "matrix_dimension" => N,
    "target_particles" => Ne,
    "maxdim" => maxdim,
    "cutoffs" => cutoffs,
    "campaign_itensors_tol" => spec.params.itensors_tol,
    "itensors_tol" => params.itensors_tol,
    "sp2_idempotency_tolerance" => sp2_idempotency_tolerance,
    "sp2_relative_trace_tolerance" => (
        isnothing(sp2_relative_trace_tolerance) ?
        "automatic" : sp2_relative_trace_tolerance
    ),
    "sp2_absolute_trace_tolerance" => (
        isnothing(sp2_absolute_trace_tolerance) ?
        "automatic" : sp2_absolute_trace_tolerance
    ),
    "spectral_lower" => bounds[1],
    "spectral_upper" => bounds[2],
    "exact_lambda_min" => first(eigenpairs.values),
    "exact_lambda_max" => last(eigenpairs.values),
    "exact_homo" => eigenpairs.values[Ne],
    "exact_lumo" => eigenpairs.values[Ne + 1],
    "exact_fermi_gap" => eigenpairs.values[Ne + 1] - eigenpairs.values[Ne],
    "dense_reference_time_s" => dense_reference.time,
    "dense_reference_allocations_bytes" => dense_reference.bytes,
    "julia_version" => string(VERSION),
    "threads" => Threads.nthreads(),
    "git_revision" => git_revision(repo_root),
    "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
    "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
    "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
    "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", "local"),
)
if !isnothing(purification)
    result = purification.value
    merge!(metadata, Dict(
        "sp2_converged" => result.converged,
        "sp2_termination_reason" => string(result.termination_reason),
        "sp2_iterations" => result.iterations,
        "sp2_trace" => result.trace,
        "sp2_trace_error" => result.trace_error,
        "sp2_idempotency_residual" => result.idempotency_residual,
        "sp2_final_max_chi" => result.final_bond_dimension,
        "sp2_purification_time_s" => purification.time,
        "sp2_purification_allocations_bytes" => purification.bytes,
    ))
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

println(
    "Projector compression comparison: label=$(spec.label) mode=$mode " *
    "maxdim=$maxdim cutoffs=$(join(cutoffs, ','))",
)
println("Result directory: $output")
