#!/usr/bin/env julia

"""Run one large fixed-H QTT KPM probe-register map without dense references.

The calculation builds the density, directed nearest-neighbour bond fields,
Hartree potential, and Fock MPO from a single binary Hadamard probe register.
It intentionally performs no SCF update: this is the scale-calibration step
before using these fields in a self-consistent loop.
"""

using Dates
using LinearAlgebra
using Printf
using Serialization
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) in (8, 10) || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS PROBES TRAJECTORY_MAXDIM CUTOFF BACKEND [FIELD_MAXDIM SAVE_FIELDS]",
)
campaign_file, task, output = abspath(ARGS[1]), parse(Int, ARGS[2]), abspath(ARGS[3])
moments, probes, maxdim = parse(Int, ARGS[4]), parse(Int, ARGS[5]), parse(Int, ARGS[6])
cutoff, backend = parse(Float64, ARGS[7]), Symbol(ARGS[8])
field_maxdim = length(ARGS) == 10 ? parse(Int, ARGS[9]) : maxdim
save_fields = length(ARGS) == 10 ? parse(Bool, ARGS[10]) : false
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
ispow2(probes) || error("PROBES must be a power of two")
moments >= 2 && maxdim > 0 && field_maxdim > 0 && cutoff > 0 || error("invalid numerical settings")
backend == :cuda && CUDA.functional() || error("this diagnostic requires functional CUDA")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
spec = Main.campaign.runs[task]
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")
N, levels = 1 << params.L, params.L
probes <= N || error("probe count exceeds N=$N")
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
scale, center = (upper - lower) / 2, (upper + lower) / 2

effective_params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=cutoff, itensors_maxdim=maxdim,
    density=params.density, purification_steps=params.purification_steps,
    scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
    scf_max_iterations=params.scf_max_iterations,
)
system = System(effective_params)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)
H_scaled = +(H_effective / scale, (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=maxdim)
register, probe_sites, _ = MPO_MeanField._qtt_hadamard_probe_register_mps(
    system.sites, trailing_zeros(probes),
)
H_register = MPO_MeanField._qtt_extend_mpo_with_probe_register(
    H_scaled, system.sites, probe_sites,
)

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    line = "[$(Dates.now())] $message"
    println(line); flush(stdout)
    open(progress_path, "a") do io
        println(io, line)
    end
end

progress("fixed-H register map: N=$N M=$moments R=$probes trajectory_chi=$maxdim field_chi=$field_maxdim")
# Compile the GPU recurrence outside the measurement.
warm = MPO_MeanField._qtt_mps_chebyshev_apply(
    ITensors.adapt(CUDA.CuArray, H_register), ITensors.adapt(CUDA.CuArray, register),
    [1.0, 0.0]; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize(); warm = nothing; CUDA.reclaim()
H_device = ITensors.adapt(CUDA.CuArray, H_register)
register_device = ITensors.adapt(CUDA.CuArray, register)
CUDA.synchronize()

moment_time = @elapsed register_moments = MPO_MeanField._qtt_probe_register_moments(
    H_device, register_device, moments; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
mu_result = MPO_MeanField._kpm_mu(register_moments, round(Int, params.density * N);
    tolerance=1e-8)
coefficients = MPO_MeanField._kpm_coefficients(moments, mu_result.scaled_mu)
progress(@sprintf("moments complete: scaled mu=%.8e", mu_result.scaled_mu))

propagation_time = @elapsed propagated_device = MPO_MeanField._qtt_mps_chebyshev_apply(
    H_device, register_device, coefficients; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
propagated = ITensors.cpu(propagated_device.state)
H_device = nothing; register_device = nothing; propagated_device = nothing
CUDA.reclaim()
progress(@sprintf("projector recurrence complete: chi=%d", maxlinkdim(propagated)))

density_time = @elapsed density = MPO_MeanField._qtt_density_mps_from_probe_register(
    propagated, register, system.sites, probe_sites; cutoff=cutoff, maxdim=field_maxdim,
)
TR, TL, TU, TD = system.translations
bond_fields = nothing
bond_time = @elapsed bond_fields = begin
    right = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TR, system.sites, probe_sites; cutoff=cutoff, maxdim=field_maxdim)
    left = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TL, system.sites, probe_sites; cutoff=cutoff, maxdim=field_maxdim)
    up = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TU, system.sites, probe_sites; cutoff=cutoff, maxdim=field_maxdim)
    down = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        propagated, register, TD, system.sites, probe_sites; cutoff=cutoff, maxdim=field_maxdim)
    horizontal = +(0.5 * right, 0.5 * apply(TR, left; cutoff=cutoff, maxdim=field_maxdim);
        cutoff=cutoff, maxdim=field_maxdim)
    vertical = +(0.5 * up, 0.5 * apply(TU, down; cutoff=cutoff, maxdim=field_maxdim);
        cutoff=cutoff, maxdim=field_maxdim)
    (; horizontal, vertical)
end
hartree_time = @elapsed hartree = Float64(params.U) * apply(
    square_neighbour_adjacency_mpo(system.sites), density; cutoff=cutoff, maxdim=field_maxdim,
)
fock = nothing
fock_time = @elapsed begin
    fock, _ = MPO_MeanField._qtt_square_fock_mpo_from_bond_mps(
        bond_fields.horizontal, bond_fields.vertical, system.sites, system.translations,
        Float64(params.U); cutoff=cutoff, maxdim=field_maxdim,
    )
end
progress("local density, bond, Hartree, and Fock fields complete")

side = 1 << (levels ÷ 2)
sample_coordinates = unique([(0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1),
    (side ÷ 2, side ÷ 2), (side ÷ 4, side ÷ 4), (3side ÷ 4, 3side ÷ 4)])
bra_states, ket_states = MPO_MeanField.precompute_qtt_states(system.sites)
open(joinpath(output, "field_samples.csv"), "w") do io
    println(io, "x,y,site,density,hartree,horizontal_bond,vertical_bond")
    for (x, y) in sample_coordinates
        site = square_lattice_index(x, y, levels)
        index = site - 1
        horizontal = MPO_MeanField._qtt_mps_amplitude(bond_fields.horizontal, system.sites, index)
        vertical = MPO_MeanField._qtt_mps_amplitude(bond_fields.vertical, system.sites, index)
        @printf(io, "%d,%d,%d,%.16g,%.16g,%.16g,%.16g\n", x, y, site,
            MPO_MeanField._qtt_mps_amplitude(density, system.sites, index),
            MPO_MeanField._qtt_mps_amplitude(hartree, system.sites, index), horizontal, vertical)
    end
end
ones = MPS([ITensor([1.0, 1.0], site) for site in system.sites])
trace = real(inner(ones, density))
field_artifact = joinpath(output, "fields.jls")
if save_fields
    serialize(field_artifact, (; density, horizontal=bond_fields.horizontal,
        vertical=bond_fields.vertical, hartree))
end
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => spec.label,
        "diagnostic" => "large_fixed_field_qtt_probe_register_local_hf",
        "matrix_dimension" => N, "moments" => moments, "probes" => probes,
        "maxdim" => maxdim, "trajectory_maxdim" => maxdim,
        "field_maxdim" => field_maxdim, "cutoff" => cutoff, "backend" => "cuda",
        "spectral_lower" => lower, "spectral_upper" => upper,
        "chemical_potential" => center + scale * mu_result.scaled_mu,
        "scaled_chemical_potential" => mu_result.scaled_mu,
        "register_moment_time_s" => moment_time,
        "register_propagation_time_s" => propagation_time,
        "density_time_s" => density_time, "bond_time_s" => bond_time,
        "hartree_time_s" => hartree_time, "fock_time_s" => fock_time,
        "register_max_chi" => maxlinkdim(propagated),
        "density_max_chi" => maxlinkdim(density),
        "horizontal_bond_max_chi" => maxlinkdim(bond_fields.horizontal),
        "vertical_bond_max_chi" => maxlinkdim(bond_fields.vertical),
        "hartree_max_chi" => maxlinkdim(hartree), "fock_max_chi" => maxlinkdim(fock),
        "density_trace" => trace,
        "density_trace_error" => abs(trace - params.density * N),
        "field_artifact" => save_fields ? field_artifact : "",
        "note" => "Fixed initial field only; no ED/vector reference and no SCF update.",
    ))
end
progress("fixed-H QTT-register local-field map complete: $output")
