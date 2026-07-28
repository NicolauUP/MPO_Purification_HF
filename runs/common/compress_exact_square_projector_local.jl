#!/usr/bin/env julia

"""Compress an exact occupied projector and compare local observables.

Unlike the small-system compression diagnostic, this script never converts a
candidate MPO back to a dense matrix. The exact dense diagonalization and
projector are formed once. Each requested QTT cap is then checked through its
trace, every site density, and every physical nearest-neighbor bond order.
"""

using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using TOML

length(ARGS) == 6 || error(
    "usage: compress_exact_square_projector_local.jl " *
    "CAMPAIGN_FILE TASK_INDEX MAXDIMS CUTOFF OUTPUT_DIRECTORY BLAS_THREADS",
)
campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
maxdims = parse.(Int, split(ARGS[3], ','))
cutoff = parse(Float64, ARGS[4])
output = abspath(ARGS[5])
blas_threads = parse(Int, ARGS[6])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
all(>(0), maxdims) || error("all MAXDIMS must be positive")
issorted(maxdims) || error("MAXDIMS must be in increasing order")
length(unique(maxdims)) == length(maxdims) || error("MAXDIMS must be unique")
isfinite(cutoff) && cutoff > 0 || error("CUTOFF must be finite and positive")
blas_threads > 0 || error("BLAS_THREADS must be positive")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
params = spec.params
params isa ParametersSquare || error("comparison requires ParametersSquare")

BLAS.set_num_threads(blas_threads)

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
rms(values) = sqrt(sum(abs2, values) / length(values))
mean_link_dimension(mpo::MPO) = length(mpo) == 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) /
    (length(mpo) - 1)

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

function interleaved_projector_tensor(matrix::AbstractMatrix, sites)
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
    return ITensor(values, physical_indices...)
end

function measure_local_projector(rho::MPO, sys::System, bonds)
    N = 2^sys.params.L
    density = [
        real(MatrixChecker(
            rho, sys.sites, site, site, sys.bra_states, sys.ket_states,
        )) for site in 1:N
    ]
    bond_order = ComplexF64[
        MatrixChecker(
            rho, sys.sites, site, neighbour, sys.bra_states, sys.ket_states,
        ) for (site, neighbour, _) in bonds
    ]
    return (; density, bond_order)
end

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
started_at = now(UTC)
N = 2^params.L
Ne = round(Int, params.density * N)
bonds = collect(square_undirected_bonds(params.L))
H = direct_initial_hamiltonian(params)
dense_reference = @timed begin
    eigenpairs = eigen(Symmetric(H))
    occupied = @view eigenpairs.vectors[:, 1:Ne]
    projector = Matrix(occupied * occupied')
    (; eigenpairs, projector)
end
eigenpairs = dense_reference.value.eigenpairs
exact_projector = dense_reference.value.projector
exact_density = real.(diag(exact_projector))
exact_bond_order = ComplexF64[
    exact_projector[site, neighbour] for (site, neighbour, _) in bonds
]
exact_energy = sum(@view eigenpairs.values[1:Ne])

sys = System(params)
tensor_build = @timed interleaved_projector_tensor(exact_projector, sys.sites)
projector_tensor = tensor_build.value

open(joinpath(output, "spectrum.csv"), "w") do io
    write_csv_row(io, ("state", "eigenvalue", "occupation"))
    for state in eachindex(eigenpairs.values)
        write_csv_row(io, (
            state, eigenpairs.values[state], state <= Ne ? 1 : 0,
        ))
    end
end

open(joinpath(output, "summary.csv"), "w") do summary
    write_csv_row(summary, (
        "cutoff", "requested_maxdim", "max_chi", "mean_chi", "trace",
        "trace_error", "density_max_abs_error", "density_rms_error",
        "horizontal_bond_max_abs_error", "horizontal_bond_rms_error",
        "vertical_bond_max_abs_error", "vertical_bond_rms_error",
        "one_body_energy", "one_body_energy_error", "compression_time_s",
        "compression_allocations_bytes", "measurement_time_s",
        "measurement_allocations_bytes",
    ))
    open(joinpath(output, "bond_dimensions.csv"), "w") do dimensions
        write_csv_row(dimensions, ("requested_maxdim", "bond", "dimension"))
        open(joinpath(output, "density_profiles.csv"), "w") do profiles
            write_csv_row(profiles, (
                "requested_maxdim", "site", "x", "y", "exact", "compressed",
                "error",
            ))
            for maxdim in maxdims
                compression = @timed MPO(
                    copy(projector_tensor), sys.sites;
                    cutoff=cutoff, maxdim=maxdim,
                )
                candidate = compression.value
                measurement = @timed measure_local_projector(candidate, sys, bonds)
                density = measurement.value.density
                bond_order = measurement.value.bond_order
                density_error = density - exact_density
                bond_error = bond_order - exact_bond_order
                horizontal = findall(
                    index -> bonds[index][3] == :horizontal, eachindex(bonds),
                )
                vertical = findall(
                    index -> bonds[index][3] == :vertical, eachindex(bonds),
                )
                onsite_energy = sum(H[site, site] * density[site] for site in 1:N)
                hopping_energy = sum(
                    2real(H[site, neighbour] * conj(bond_order[index]))
                    for (index, (site, neighbour, _)) in enumerate(bonds)
                )
                energy = onsite_energy + hopping_energy
                write_csv_row(summary, (
                    cutoff, maxdim, maxlinkdim(candidate),
                    mean_link_dimension(candidate), real(tr(candidate)),
                    abs(real(tr(candidate)) - Ne),
                    maximum(abs, density_error), rms(density_error),
                    maximum(abs, bond_error[horizontal]),
                    rms(bond_error[horizontal]),
                    maximum(abs, bond_error[vertical]),
                    rms(bond_error[vertical]), energy, energy - exact_energy,
                    compression.time, compression.bytes,
                    measurement.time, measurement.bytes,
                ))
                for bond in 1:(length(candidate) - 1)
                    write_csv_row(dimensions, (
                        maxdim, bond, dim(linkind(candidate, bond)),
                    ))
                end
                for site in 1:N
                    x, y = square_lattice_decoder(site - 1, params.L)
                    write_csv_row(profiles, (
                        maxdim, site, x, y, exact_density[site],
                        density[site], density_error[site],
                    ))
                end
                flush(summary)
                flush(dimensions)
                flush(profiles)
            end
        end
    end
end

open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "started_at" => string(started_at),
        "finished_at" => string(now(UTC)),
        "campaign" => string(campaign.name),
        "label" => string(spec.label),
        "task_index" => task_index,
        "matrix_dimension" => N,
        "target_particles" => Ne,
        "cutoff" => cutoff,
        "maxdims" => maxdims,
        "blas_threads" => BLAS.get_num_threads(),
        "exact_lambda_min" => first(eigenpairs.values),
        "exact_lambda_max" => last(eigenpairs.values),
        "exact_homo" => eigenpairs.values[Ne],
        "exact_lumo" => eigenpairs.values[Ne + 1],
        "exact_fermi_gap" => eigenpairs.values[Ne + 1] - eigenpairs.values[Ne],
        "dense_diagonalization_time_s" => dense_reference.time,
        "dense_diagonalization_allocations_bytes" => dense_reference.bytes,
        "tensor_build_time_s" => tensor_build.time,
        "tensor_build_allocations_bytes" => tensor_build.bytes,
    ))
end

println("Exact square projector compression completed.")
println("N=$N Ne=$Ne exact gap=$(eigenpairs.values[Ne + 1] - eigenpairs.values[Ne])")
println("Result directory: $output")
