#!/usr/bin/env julia

"""Compare fixed-Hamiltonian square SP2 with exact diagonalization.

The dense Hamiltonian is reconstructed directly from the campaign functions.
The MPO is not densified: every site density and physical nearest-neighbor
bond order is measured directly and compared with the exact occupied
projector. This keeps the production 32x32 comparison tractable.
"""

using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using SHA
using TOML

length(ARGS) in 4:8 || error(
    "usage: compare_square_sp2_dense_projector.jl CAMPAIGN_FILE TASK_INDEX MAXDIM OUTPUT_DIRECTORY [SP2_IDEMPOTENCY_TOLERANCE [BACKEND [ITENSORS_TOL [SP2_RELATIVE_TRACE_TOLERANCE]]]]",
)
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
maxdim = tryparse(Int, ARGS[3])
output = abspath(ARGS[4])
sp2_idempotency_tolerance = length(ARGS) >= 5 ? tryparse(Float64, ARGS[5]) : 1e-3
backend = length(ARGS) >= 6 ? Symbol(lowercase(ARGS[6])) : :cpu
itensors_tol_override = length(ARGS) >= 7 ? tryparse(Float64, ARGS[7]) : nothing
sp2_relative_trace_tolerance = length(ARGS) == 8 ? tryparse(Float64, ARGS[8]) : nothing
isnothing(task_index) && error("TASK_INDEX must be an integer")
isnothing(maxdim) && error("MAXDIM must be an integer")
isnothing(sp2_idempotency_tolerance) &&
    error("SP2_IDEMPOTENCY_TOLERANCE must be a number")
maxdim > 0 || error("MAXDIM must be positive")
isfinite(sp2_idempotency_tolerance) && sp2_idempotency_tolerance > 0 ||
    error("SP2_IDEMPOTENCY_TOLERANCE must be finite and positive")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
if length(ARGS) >= 7
    isnothing(itensors_tol_override) && error("ITENSORS_TOL must be a number")
    isfinite(itensors_tol_override) && itensors_tol_override > 0 ||
        error("ITENSORS_TOL must be finite and positive")
end
if length(ARGS) == 8
    isnothing(sp2_relative_trace_tolerance) &&
        error("SP2_RELATIVE_TRACE_TOLERANCE must be a number")
    isfinite(sp2_relative_trace_tolerance) && sp2_relative_trace_tolerance > 0 ||
        error("SP2_RELATIVE_TRACE_TOLERANCE must be finite and positive")
end
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("comparison requires ParametersSquare")
ispath(output) && error("refusing to overwrite existing output directory: $output")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end
rms(values) = sqrt(sum(abs2, values) / length(values))

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

to_device = backend == :cuda ?
    (value -> ITensors.adapt(CUDA.CuArray, value)) : identity
to_host = backend == :cuda ? ITensors.cpu : identity
synchronize_backend() = backend == :cuda ? CUDA.synchronize() : nothing
device_name = backend == :cuda ? CUDA.name(CUDA.device()) : "CPU"
device_total_memory = backend == :cuda ? CUDA.total_memory() : 0
device_free_memory_before = backend == :cuda ? CUDA.free_memory() : 0

function backend_timed(f)
    synchronize_backend()
    measurement = @timed begin
        value = f()
        synchronize_backend()
        value
    end
    return measurement
end

function with_numerical_controls(
    params::ParametersSquare,
    maxdim::Int,
    itensors_tol_override::Union{Nothing,Float64},
)
    return ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol,
        itensors_tol=isnothing(itensors_tol_override) ?
            params.itensors_tol : itensors_tol_override,
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

function measure_local_projector(rho::MPO, sys::System, bonds)
    N = 2^sys.params.L
    density = Vector{Float64}(undef, N)
    for site in 1:N
        density[site] = real(MatrixChecker(
            rho, sys.sites, site, site, sys.bra_states, sys.ket_states,
        ))
    end
    bond_order = Vector{ComplexF64}(undef, length(bonds))
    for (index, (site, neighbour, _)) in enumerate(bonds)
        bond_order[index] = MatrixChecker(
            rho, sys.sites, site, neighbour, sys.bra_states, sys.ket_states,
        )
    end
    return (; density, bond_order)
end

function one_body_energy(H, density, bonds, bond_order)
    onsite = sum(H[site, site] * density[site] for site in eachindex(density))
    hopping = sum(
        2real(H[site, neighbour] * conj(bond_order[index]))
        for (index, (site, neighbour, _)) in enumerate(bonds)
    )
    return (; onsite, hopping, total=onsite + hopping)
end

function checkerboard_order(density, L)
    total = sum(eachindex(density)) do site
        x, y = square_lattice_decoder(site - 1, L)
        (iseven(x + y) ? 1.0 : -1.0) * density[site]
    end
    return total / length(density)
end

params = with_numerical_controls(spec.params, maxdim, itensors_tol_override)
bounds = validate_spectral_bounds(spec.spectral_bounds...)
N = 2^params.L
Ne = round(Int, params.density * N)
bonds = collect(square_undirected_bonds(params.L))
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
started_at = now(UTC)
mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))

H_dense = direct_initial_hamiltonian(params)
dense_calculation = @timed begin
    eigenpairs = eigen(Symmetric(H_dense))
    occupied = @view eigenpairs.vectors[:, 1:Ne]
    projector = Matrix(occupied * occupied')
    (; eigenpairs, projector)
end
eigenpairs = dense_calculation.value.eigenpairs
exact_projector = dense_calculation.value.projector
exact_density = real.(diag(exact_projector))
exact_bond_order = ComplexF64[
    exact_projector[site, neighbour] for (site, neighbour, _) in bonds
]

sys = System(params)
H_mpo_cpu = +(sys.H0, sys.VH, sys.VF;
    cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
)
transfer_to_device = backend_timed() do
    sys.H0 = to_device(sys.H0)
    sys.VH = to_device(sys.VH)
    sys.VF = to_device(sys.VF)
end
rho0_calculation = backend_timed() do
    construct_rho_0(
        sys, params, bounds...;
        method=:sp2,
        verify_spectral_bounds=false,
        to_gpu=to_device,
    )
end
rho0 = rho0_calculation.value
device_scalar_type = ITensors.scalartype(rho0)
backend == :cuda && device_scalar_type != Float64 && error(
    "CUDA SP2 comparison requires Float64 MPO tensors, got $device_scalar_type",
)
purification = backend_timed() do
    open(joinpath(output, "sp2_progress.txt"), "w") do progress
        perform_purification(
            rho0, params;
            method=:sp2,
            verbose=1,
            io=progress,
            overwrite_progress=false,
            sp2_idempotency_tolerance=sp2_idempotency_tolerance,
            sp2_trace_tolerance=isnothing(sp2_relative_trace_tolerance) ?
                nothing : Ne * sp2_relative_trace_tolerance,
            spectral_bounds=bounds,
            spectral_bounds_validation=:supplied_analytical,
        )
    end
end
result = purification.value

transfer_to_host = backend_timed(() -> to_host(result.rho))
rho_host = transfer_to_host.value
measurement = @timed measure_local_projector(rho_host, sys, bonds)
sp2_density = measurement.value.density
sp2_bond_order = measurement.value.bond_order
stationarity = @timed MPO_MeanField._scf_commutator_residual(
    H_mpo_cpu, rho_host, params,
)
device_free_memory_after = backend == :cuda ? CUDA.free_memory() : 0

if backend == :cuda
    open(joinpath(output, "cuda_status.txt"), "w") do io
        CUDA.versioninfo(io)
        println(io)
        CUDA.pool_status(io)
    end
end

density_error = sp2_density - exact_density
bond_error = sp2_bond_order - exact_bond_order
horizontal_indices = findall(index -> bonds[index][3] == :horizontal, eachindex(bonds))
vertical_indices = findall(index -> bonds[index][3] == :vertical, eachindex(bonds))
exact_energy = one_body_energy(H_dense, exact_density, bonds, exact_bond_order)
sp2_energy = one_body_energy(H_dense, sp2_density, bonds, sp2_bond_order)
eigenvalue_energy = sum(@view eigenpairs.values[1:Ne])

open(joinpath(output, "site_density_comparison.csv"), "w") do io
    write_csv_row(io, ("site", "x", "y", "exact", "sp2", "error"))
    for site in 1:N
        x, y = square_lattice_decoder(site - 1, params.L)
        write_csv_row(io, (
            site, x, y, exact_density[site], sp2_density[site], density_error[site],
        ))
    end
end
open(joinpath(output, "bond_order_comparison.csv"), "w") do io
    write_csv_row(io, (
        "site_left", "site_right", "orientation", "exact_real", "exact_imag",
        "sp2_real", "sp2_imag", "error_abs",
    ))
    for (index, (site, neighbour, orientation)) in enumerate(bonds)
        write_csv_row(io, (
            site, neighbour, orientation,
            real(exact_bond_order[index]), imag(exact_bond_order[index]),
            real(sp2_bond_order[index]), imag(sp2_bond_order[index]),
            abs(bond_error[index]),
        ))
    end
end
open(joinpath(output, "spectrum.csv"), "w") do io
    write_csv_row(io, ("state", "eigenvalue", "occupation"))
    for state in eachindex(eigenpairs.values)
        write_csv_row(io, (
            state, eigenpairs.values[state], state <= Ne ? 1 : 0,
        ))
    end
end

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => string(campaign.name),
        "label" => string(spec.label),
        "task_index" => task_index,
        "matrix_dimension" => N,
        "target_particles" => Ne,
        "itensors_tol" => params.itensors_tol,
        "itensors_maxdim" => params.itensors_maxdim,
        "backend" => string(backend),
        "device_name" => device_name,
        "device_scalar_type" => string(device_scalar_type),
        "device_total_memory_bytes" => device_total_memory,
        "device_free_memory_before_bytes" => device_free_memory_before,
        "device_free_memory_after_bytes" => device_free_memory_after,
        "sp2_idempotency_tolerance" => sp2_idempotency_tolerance,
        "sp2_relative_trace_tolerance" => isnothing(sp2_relative_trace_tolerance) ?
            "default" : sp2_relative_trace_tolerance,
        "spectral_lower" => bounds[1],
        "spectral_upper" => bounds[2],
        "exact_lambda_min" => first(eigenpairs.values),
        "exact_lambda_max" => last(eigenpairs.values),
        "exact_homo" => eigenpairs.values[Ne],
        "exact_lumo" => eigenpairs.values[Ne + 1],
        "exact_fermi_gap" => eigenpairs.values[Ne + 1] - eigenpairs.values[Ne],
        "sp2_converged" => result.converged,
        "sp2_termination_reason" => string(result.termination_reason),
        "sp2_iterations" => result.iterations,
        "sp2_selected_iteration" => result.selected_iteration,
        "sp2_trace" => result.trace,
        "sp2_trace_error" => abs(result.trace - Ne),
        "sp2_idempotency_residual" => result.idempotency_residual,
        "sp2_hermiticity_residual" => result.hermiticity_residual,
        "sp2_stationarity_residual" => stationarity.value,
        "sp2_final_max_chi" => result.final_bond_dimension,
        "sp2_work_max_chi" => result.work.max_bond_dimension,
        "sp2_work_mean_chi" => result.work.mean_bond_dimension,
        "density_max_abs_error" => maximum(abs, density_error),
        "density_rms_error" => rms(density_error),
        "horizontal_bond_max_abs_error" => maximum(abs, bond_error[horizontal_indices]),
        "horizontal_bond_rms_error" => rms(bond_error[horizontal_indices]),
        "vertical_bond_max_abs_error" => maximum(abs, bond_error[vertical_indices]),
        "vertical_bond_rms_error" => rms(bond_error[vertical_indices]),
        "exact_checkerboard_order" => checkerboard_order(exact_density, params.L),
        "sp2_checkerboard_order" => checkerboard_order(sp2_density, params.L),
        "exact_onsite_energy" => exact_energy.onsite,
        "sp2_onsite_energy" => sp2_energy.onsite,
        "exact_hopping_energy" => exact_energy.hopping,
        "sp2_hopping_energy" => sp2_energy.hopping,
        "exact_one_body_energy" => exact_energy.total,
        "exact_eigenvalue_sum" => eigenvalue_energy,
        "exact_energy_reconstruction_error" => exact_energy.total - eigenvalue_energy,
        "sp2_one_body_energy" => sp2_energy.total,
        "one_body_energy_error" => sp2_energy.total - exact_energy.total,
        "dense_diagonalization_time_s" => dense_calculation.time,
        "dense_diagonalization_allocations_bytes" => dense_calculation.bytes,
        "purification_time_s" => purification.time,
        "purification_allocations_bytes" => purification.bytes,
        "initialization_time_s" => rho0_calculation.time,
        "initialization_allocations_bytes" => rho0_calculation.bytes,
        "transfer_to_device_time_s" => transfer_to_device.time,
        "transfer_to_host_time_s" => transfer_to_host.time,
        "local_measurement_time_s" => measurement.time,
        "local_measurement_allocations_bytes" => measurement.bytes,
        "stationarity_time_s" => stationarity.time,
        "stationarity_allocations_bytes" => stationarity.bytes,
    ))
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "started_at" => string(started_at),
        "finished_at" => string(now(UTC)),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "git_revision" => git_revision(repo_root),
        "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
        "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
        "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", "local"),
    ))
end

println(
    "Fixed square SP2/dense comparison: label=$(spec.label) backend=$backend " *
    "device=$(device_name) maxdim=$maxdim " *
    "idempotency_tolerance=$sp2_idempotency_tolerance",
)
println("SP2 converged=$(result.converged) iterations=$(result.iterations)")
println("density max error=$(maximum(abs, density_error))")
println("bond max error=$(maximum(abs, bond_error))")
println("one-body energy error=$(sp2_energy.total - exact_energy.total)")
println("Result directory: $output")
