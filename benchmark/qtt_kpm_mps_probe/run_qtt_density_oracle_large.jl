#!/usr/bin/env julia

"""Large fixed-H QTT density-oracle experiment for legacy square campaigns.

Unlike the 32x32 validation driver, this script never diagonalizes H, never
forms a sparse/vector KPM block, and never materializes n[i]. It assumes a
half-filled bipartite hopping model with a staggered seed, for which mu=0.
The supplied campaign spectral interval must enclose the fixed H0+seed.
"""

using Dates
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS

length(ARGS) == 7 || error(
    "usage: $(PROGRAM_FILE) LEGACY_CAMPAIGN.jl TASK OUTPUT MOMENTS PROBES MAXDIM TCI_TOL",
)
campaign_file = abspath(ARGS[1]); task = parse(Int, ARGS[2]); output = abspath(ARGS[3])
moments = parse(Int, ARGS[4]); probes = parse(Int, ARGS[5]); maxdim = parse(Int, ARGS[6]); tci_tol = parse(Float64, ARGS[7])
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments > 0 && probes > 0 && maxdim > 0 && tci_tol > 0 || error("moments, probes, maxdim, and tci_tol must be positive")

Base.include(Main, campaign_file)
@isdefined(campaign) || error("campaign file did not define `campaign`")
1 <= task <= length(campaign.runs) || error("TASK is outside campaign.runs")
spec = campaign.runs[task]
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")
N = 1 << params.L
probes <= N || error("probes exceed N=$N")
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
abs(lower + upper) <= 1e-12 || error("large oracle currently requires a symmetric spectral interval")

# μ=0 is physically specific: it is valid for the open bipartite hopping
# Hamiltonian plus staggered checkerboard seed used by this campaign.
coefficients = MPO_MeanField._kpm_coefficients(moments, 0.0)
effective_params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
    itensors_maxdim=maxdim, density=params.density,
    purification_steps=params.purification_steps, scf_mixing=params.scf_mixing,
    scf_tol=params.scf_tol, scf_max_iterations=params.scf_max_iterations,
)
system = System(effective_params)
H_eff = +(system.H0, system.VH, system.VF; cutoff=effective_params.itensors_tol, maxdim=maxdim)
halfwidth = (upper - lower) / 2
H_scaled = H_eff / halfwidth

propagated = MPS[]
propagation_time = @elapsed for row in 0:(probes - 1)
    probe = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row)
    result = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_scaled, probe, coefficients;
        cutoff=effective_params.itensors_tol, maxdim=maxdim,
    )
    push!(propagated, result.state)
end
signs = [MPO_MeanField._qtt_hadamard_probe_vector(params.L, row) for row in 0:(probes - 1)]
query_count = Ref(0)
function density_oracle(index)
    qtt_index = Int(index)
    query_count[] += 1
    return sum(signs[column][qtt_index + 1] *
               MPO_MeanField._qtt_mps_amplitude(propagated[column], system.sites, qtt_index)
               for column in eachindex(propagated)) / probes
end
tci_time = @elapsed _, _, density_mps = MPO_MeanField.Quantics_TCI(
    density_oracle, Float64, system.sites, tci_tol,
)
ones = MPS([ITensor([1.0, 1.0], site) for site in system.sites])
estimated_trace = real(inner(ones, density_mps))

mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => spec.label,
        "matrix_dimension" => N, "moments" => moments, "probes" => probes,
        "maxdim" => maxdim, "mpo_cutoff" => effective_params.itensors_tol,
        "tci_tolerance" => tci_tol, "chemical_potential" => 0.0,
        "spectral_lower" => lower, "spectral_upper" => upper,
        "propagation_time_s" => propagation_time, "tci_time_s" => tci_time,
        "tci_oracle_evaluations" => query_count[], "density_tci_max_chi" => maxlinkdim(density_mps),
        "estimated_particle_number" => estimated_trace,
        "estimated_particle_number_error" => abs(estimated_trace - params.density * N),
        "validation_note" => "fixed-H large-system oracle: no dense/vector reference or full density array was constructed",
    ))
end
println("large QTT density oracle complete: $output")
