#!/usr/bin/env julia

"""Validate register-derived density, bonds, and Hartree fields at fixed H.

The QTT register and sparse-vector KPM paths use identical Walsh columns,
exact spectral scaling, chemical potential, and Chebyshev polynomial. This is
the first HF-map component test: it deliberately stops before assembling a
Fock MPO or updating the Hamiltonian.
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 8 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS PROBES MAXDIM CUTOFF BACKEND",
)
campaign_file, task, output = abspath(ARGS[1]), parse(Int, ARGS[2]), abspath(ARGS[3])
moments, probes, maxdim = parse(Int, ARGS[4]), parse(Int, ARGS[5]), parse(Int, ARGS[6])
cutoff, backend = parse(Float64, ARGS[7]), Symbol(ARGS[8])
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
ispow2(probes) || error("PROBES must be a power of two")
moments >= 2 && maxdim > 0 && cutoff > 0 || error("invalid numerical settings")
backend == :cuda && CUDA.functional() || error("this diagnostic requires functional CUDA")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
case = campaign_case(Main.campaign, task)
model = case.model
model isa SquareModel || error("diagnostic requires SquareModel")
nx, ny = model.size
nx == ny || error("diagnostic currently requires a square model")
N, levels = nx * ny, qtt_levels(model)
probes <= N || error("probe count exceeds matrix dimension")
representation = QTTSettings(encoding=:interleaved, tci_tol=case.representation.tci_tol,
    cutoff=cutoff, maxdim=maxdim)
parameters = legacy_parameters(model, representation, case.solver)
system = System(parameters)

# Reference matrix is the same seeded one-body Hamiltonian used by the MPO.
data = MPO_MeanField._kpm_data(model)
Hrow = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
qtt_to_row = [begin
    x, y = square_lattice_decoder(index, levels); x + nx * y + 1
end for index in 0:(N - 1)]
Hvector = sparse(Hrow[qtt_to_row, qtt_to_row])
dense_diagonalization_time = @elapsed spectrum = eigen(Symmetric(Matrix(Hvector)))
eigenvalues = spectrum.values
Ne = round(Int, model.filling * N)
center = (first(eigenvalues) + last(eigenvalues)) / 2
scale = 1.05 * (last(eigenvalues) - first(eigenvalues)) / 2
mu = (eigenvalues[Ne] + eigenvalues[Ne + 1]) / 2
coefficients = MPO_MeanField._kpm_coefficients(moments, (mu - center) / scale)
Hvector = (Hvector - center * sparse(I, N, N)) / scale

codes = 0:(probes - 1)
Z = [isodd(count_ones(index & code)) ? -1.0 : 1.0
     for index in 0:(N - 1), code in codes]
vector_time = @elapsed PZ = MPO_MeanField._kpm_apply(Hvector, Z, coefficients)
vector_moments = MPO_MeanField._kpm_moments(Hvector, Z, moments)
vector_mu = MPO_MeanField._kpm_mu(vector_moments, Ne; tolerance=1e-8)
# Reorder physical bonds to the QTT row ordering used by `Hvector`, `PZ`, and
# `Z` before evaluating the independent vector local estimator.
row_to_qtt = invperm(qtt_to_row)
bonds_qtt = [(row_to_qtt[i], row_to_qtt[j], orientation) for (i,j,orientation) in data.bonds]
vector_local = MPO_MeanField._kpm_local(PZ, Z, bonds_qtt)
vector_hartree, _ = MPO_MeanField._kpm_fields(vector_local.density, vector_local.order,
    (; data..., bonds=bonds_qtt), Float64(model.interaction))
exact_projector = spectrum.vectors[:, 1:Ne] * spectrum.vectors[:, 1:Ne]'
exact_density = real.(diag(exact_projector))
exact_bond = [real(exact_projector[i, j]) for (i, j, _) in bonds_qtt]
exact_hartree, _ = MPO_MeanField._kpm_fields(exact_density, exact_bond,
    (; data..., bonds=bonds_qtt), Float64(model.interaction))

Heffective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)
Hmpo = +(Heffective / scale, (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=maxdim)
register, probe_sites, _ = MPO_MeanField._qtt_hadamard_probe_register_mps(
    system.sites, trailing_zeros(probes),
)
Hregister = MPO_MeanField._qtt_extend_mpo_with_probe_register(Hmpo, system.sites, probe_sites)

# One short recurrence compiles the same CUDA code path outside the timing.
warm = MPO_MeanField._qtt_mps_chebyshev_apply(ITensors.adapt(CUDA.CuArray, Hregister),
    ITensors.adapt(CUDA.CuArray, register), coefficients[1:2]; cutoff=cutoff, maxdim=maxdim)
CUDA.synchronize(); warm = nothing; CUDA.reclaim()
Hdevice, register_device = ITensors.adapt(CUDA.CuArray, Hregister), ITensors.adapt(CUDA.CuArray, register)
CUDA.synchronize()
moment_time = @elapsed register_moments = MPO_MeanField._qtt_probe_register_moments(
    Hdevice, register_device, moments; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
register_mu = MPO_MeanField._kpm_mu(register_moments, Ne; tolerance=1e-8)
propagation_time = @elapsed propagated_device = MPO_MeanField._qtt_mps_chebyshev_apply(
    Hdevice, register_device, coefficients; cutoff=cutoff, maxdim=maxdim)
CUDA.synchronize()
propagated = ITensors.cpu(propagated_device.state)

density_time = @elapsed density = MPO_MeanField._qtt_density_mps_from_probe_register(
    propagated, register, system.sites, probe_sites; cutoff=cutoff, maxdim=maxdim)
TR, TL, TU, TD = system.translations
bond_fields = nothing
bond_time = @elapsed begin
    right = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TR, system.sites, probe_sites; cutoff=cutoff, maxdim=maxdim)
    left = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TL, system.sites, probe_sites; cutoff=cutoff, maxdim=maxdim)
    up = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TU, system.sites, probe_sites; cutoff=cutoff, maxdim=maxdim)
    down = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TD, system.sites, probe_sites; cutoff=cutoff, maxdim=maxdim)
    horizontal = +(0.5 * right, 0.5 * apply(TR, left; cutoff=cutoff, maxdim=maxdim);
        cutoff=cutoff, maxdim=maxdim)
    vertical = +(0.5 * up, 0.5 * apply(TU, down; cutoff=cutoff, maxdim=maxdim);
        cutoff=cutoff, maxdim=maxdim)
    bond_fields = (; horizontal, vertical)
end
adjacency = square_neighbour_adjacency_mpo(system.sites)
hartree_time = @elapsed hartree = Float64(model.interaction) * apply(adjacency, density;
    cutoff=cutoff, maxdim=maxdim)
fock, fock_components = nothing, nothing
fock_time = @elapsed begin
    fock, fock_components = MPO_MeanField._qtt_square_fock_mpo_from_bond_mps(
        bond_fields.horizontal, bond_fields.vertical, system.sites, system.translations,
        Float64(model.interaction); cutoff=cutoff, maxdim=maxdim,
    )
end

evaluate(psi, index) = MPO_MeanField._qtt_mps_amplitude(psi, system.sites, index)
density_qtt = [evaluate(density, index) for index in 0:(N - 1)]
hartree_qtt = [evaluate(hartree, index) for index in 0:(N - 1)]
bond_qtt = Float64[]
bond_vector = Float64[]
fock_qtt = Float64[]
bra_states, ket_states = MPO_MeanField.precompute_qtt_states(system.sites)
for (ordinal, (i, _, orientation)) in enumerate(bonds_qtt)
    push!(bond_qtt, evaluate(orientation == :horizontal ? bond_fields.horizontal : bond_fields.vertical, i - 1))
    push!(bond_vector, vector_local.order[ordinal])
    _, j, _ = bonds_qtt[ordinal]
    push!(fock_qtt, real(MPO_MeanField.MatrixChecker(
        fock, system.sites, i, j, bra_states, ket_states,
    )))
end
rms(values) = norm(values) / sqrt(length(values))
mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "diagnostic" => "fixed_field_qtt_probe_register_local_hf_estimators",
        "matrix_dimension" => N, "moments" => moments, "probes" => probes,
        "maxdim" => maxdim, "cutoff" => cutoff, "backend" => "cuda",
        "chemical_potential" => mu, "spectral_lower" => first(eigenvalues),
        "spectral_upper" => last(eigenvalues), "vector_kpm_time_s" => vector_time,
        "dense_diagonalization_time_s" => dense_diagonalization_time,
        "register_propagation_time_s" => propagation_time, "density_time_s" => density_time,
        "register_moment_time_s" => moment_time,
        "vector_scaled_chemical_potential" => vector_mu.scaled_mu,
        "register_scaled_chemical_potential" => register_mu.scaled_mu,
        "scaled_chemical_potential_difference" => register_mu.scaled_mu - vector_mu.scaled_mu,
        "bond_time_s" => bond_time, "hartree_time_s" => hartree_time,
        "fock_time_s" => fock_time,
        "register_max_chi" => maxlinkdim(propagated), "density_max_chi" => maxlinkdim(density),
        "horizontal_bond_max_chi" => maxlinkdim(bond_fields.horizontal),
        "vertical_bond_max_chi" => maxlinkdim(bond_fields.vertical), "hartree_max_chi" => maxlinkdim(hartree),
        "density_max_abs_error" => maximum(abs.(density_qtt - vector_local.density)),
        "density_rms_error" => rms(density_qtt - vector_local.density),
        "bond_max_abs_error" => maximum(abs.(bond_qtt - bond_vector)),
        "bond_rms_error" => rms(bond_qtt - bond_vector),
        "hartree_max_abs_error" => maximum(abs.(hartree_qtt - vector_hartree)),
        "hartree_rms_error" => rms(hartree_qtt - vector_hartree),
        "fock_max_chi" => maxlinkdim(fock),
        "fock_max_abs_error" => maximum(abs.(fock_qtt + Float64(model.interaction) .* bond_vector)),
        "fock_rms_error" => rms(fock_qtt + Float64(model.interaction) .* bond_vector),
        "density_trace" => sum(density_qtt), "density_trace_error" => abs(sum(density_qtt) - Ne),
        "vector_density_exact_rms_error" => rms(vector_local.density - exact_density),
        "vector_density_exact_max_abs_error" => maximum(abs.(vector_local.density - exact_density)),
        "qtt_density_exact_rms_error" => rms(density_qtt - exact_density),
        "qtt_density_exact_max_abs_error" => maximum(abs.(density_qtt - exact_density)),
        "vector_bond_exact_rms_error" => rms(vector_local.order - exact_bond),
        "vector_bond_exact_max_abs_error" => maximum(abs.(vector_local.order - exact_bond)),
        "qtt_bond_exact_rms_error" => rms(bond_qtt - exact_bond),
        "qtt_bond_exact_max_abs_error" => maximum(abs.(bond_qtt - exact_bond)),
        "vector_hartree_exact_rms_error" => rms(vector_hartree - exact_hartree),
        "vector_hartree_exact_max_abs_error" => maximum(abs.(vector_hartree - exact_hartree)),
        "qtt_hartree_exact_rms_error" => rms(hartree_qtt - exact_hartree),
        "qtt_hartree_exact_max_abs_error" => maximum(abs.(hartree_qtt - exact_hartree)),
        "vector_fock_exact_rms_error" => abs(Float64(model.interaction)) * rms(vector_local.order - exact_bond),
        "vector_fock_exact_max_abs_error" => abs(Float64(model.interaction)) * maximum(abs.(vector_local.order - exact_bond)),
        "qtt_fock_exact_rms_error" => abs(Float64(model.interaction)) * rms(fock_qtt + Float64(model.interaction) .* exact_bond),
        "qtt_fock_exact_max_abs_error" => abs(Float64(model.interaction)) * maximum(abs.(fock_qtt + Float64(model.interaction) .* exact_bond)),
    ))
end
println("fixed-field QTT register local-field diagnostic complete: $output")
