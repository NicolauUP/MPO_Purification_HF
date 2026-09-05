#!/usr/bin/env julia

"""Hybrid GPU-propagated, host-TCI QTT density-oracle ladder.

For a fixed legacy square campaign, every rank-one Hadamard MPS probe is
propagated through the KPM polynomial on CUDA. Its final compressed MPS is
then copied to host before the next probe, so device memory never retains a
batch of propagated states. The same host-resident probe states are reused for
every requested Quantics-TCI tolerance.

This is a fixed-H diagnostic: it does not run SCF and does not form a dense or
N-site density vector. The validation CSV evaluates only deterministic sample
sites after each TCI fit.
"""

using Dates
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 8 || error(
    "usage: $(PROGRAM_FILE) LEGACY_CAMPAIGN.jl TASK OUTPUT MOMENTS PROBES MAXDIM TCI_TOLS PROBE_ROW_OFFSET",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
probes = parse(Int, ARGS[5])
maxdim = parse(Int, ARGS[6])
tci_tolerances = parse.(Float64, split(ARGS[7], ','))
probe_row_offset = parse(Int, ARGS[8])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 2 || error("moments must be at least two")
probes > 0 || error("probes must be positive")
maxdim > 0 || error("maxdim must be positive")
!isempty(tci_tolerances) && all(>(0), tci_tolerances) || error("TCI tolerances must be positive")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

# Campaign closures are loaded before any numerical work. Legacy QTT MPO
# construction invokes these function-valued hoppings/seeds inside TCI.
Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
spec = Main.campaign.runs[task]
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")
N = 1 << params.L
probes + probe_row_offset <= N || error("requested probe rows exceed N=$N")
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
abs(lower + upper) <= 1e-12 || error("hybrid large oracle requires a symmetric spectral interval")

effective_params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
    itensors_maxdim=maxdim, density=params.density,
    purification_steps=params.purification_steps, scf_mixing=params.scf_mixing,
    scf_tol=params.scf_tol, scf_max_iterations=params.scf_max_iterations,
)
system = System(effective_params)
H_effective = +(system.H0, system.VH, system.VF; cutoff=effective_params.itensors_tol, maxdim=maxdim)
H_device = ITensors.adapt(CUDA.CuArray, H_effective / ((upper - lower) / 2))
coefficients = MPO_MeanField._kpm_coefficients(moments, 0.0)

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    stamped = "[$(Dates.now())] $message"
    println(stamped)
    flush(stdout)
    open(progress_path, "a") do io
        println(io, stamped)
    end
end

progress("starting GPU propagation: N=$N M=$moments R=$probes χmax=$maxdim")
propagated = MPS[]
probe_rows = collect(probe_row_offset:(probe_row_offset + probes - 1))
propagation_rows = NamedTuple[]
for (ordinal, row) in enumerate(probe_rows)
    probe_host = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row)
    probe_device = ITensors.adapt(CUDA.CuArray, probe_host)
    CUDA.synchronize()
    free_before = CUDA.free_memory()
    elapsed = @elapsed result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_device, probe_device, coefficients;
        cutoff=effective_params.itensors_tol, maxdim=maxdim,
    )
    CUDA.synchronize()
    result_host = ITensors.cpu(result_device.state)
    push!(propagated, result_host)
    final_chi = maxlinkdim(result_host)
    probe_device = nothing
    result_device = nothing
    GC.gc()
    CUDA.reclaim()
    free_after = CUDA.free_memory()
    push!(propagation_rows, (
        ordinal=ordinal, row=row, time_s=elapsed, final_max_chi=final_chi,
        free_before=free_before, free_after=free_after,
    ))
    progress(@sprintf("probe %d/%d row=%d complete: %.3f s, χ=%d, free GPU %.2f GiB",
        ordinal, probes, row, elapsed, final_chi, free_after / 2.0^30))
end

signs = [MPO_MeanField._qtt_hadamard_probe_vector(params.L, row) for row in probe_rows]
function raw_density(index::Int)
    return sum(signs[column][index + 1] *
               MPO_MeanField._qtt_mps_amplitude(propagated[column], system.sites, index)
               for column in eachindex(propagated)) / probes
end

# Fixed stratified validation points cover the QTT index space without forming
# a full density vector. They are excluded from the TCI oracle-evaluation count.
validation_count = min(1024, N)
validation_indices = unique(round.(Int, range(0, N - 1; length=validation_count)))
validation_raw = [raw_density(index) for index in validation_indices]

open(joinpath(output, "propagation.csv"), "w") do io
    println(io, "ordinal,probe_row,time_s,final_max_chi,free_gpu_before_bytes,free_gpu_after_bytes")
    for row in propagation_rows
        println(io, "$(row.ordinal),$(row.row),$(row.time_s),$(row.final_max_chi),$(row.free_before),$(row.free_after)")
    end
end

ladder_rows = NamedTuple[]
validation_columns = Dict{Float64,Vector{Float64}}()
for tolerance in tci_tolerances
    query_count = Ref(0)
    function counting_oracle(index)
        query_count[] += 1
        return raw_density(Int(index))
    end
    progress("starting host TCI fit at tolerance=$tolerance")
    elapsed = @elapsed _, _, density_mps = MPO_MeanField.Quantics_TCI(
        counting_oracle, Float64, system.sites, tolerance,
    )
    fit_queries = query_count[]
    estimated_trace = real(inner(MPS([ITensor([1.0, 1.0], site) for site in system.sites]), density_mps))
    validation_fit = [MPO_MeanField._qtt_mps_amplitude(density_mps, system.sites, index)
                      for index in validation_indices]
    validation_error = validation_fit .- validation_raw
    validation_columns[tolerance] = validation_fit
    push!(ladder_rows, (
        tolerance=tolerance, time_s=elapsed, queries=fit_queries,
        max_chi=maxlinkdim(density_mps), trace=estimated_trace,
        trace_error=abs(estimated_trace - params.density * N),
        validation_maximum=maximum(abs, validation_error),
        validation_rms=norm(validation_error) / sqrt(length(validation_error)),
    ))
    progress(@sprintf("TCI tolerance=%g complete: %.3f s, queries=%d, χ=%d, sample RMS=%.3e",
        tolerance, elapsed, fit_queries, maxlinkdim(density_mps), ladder_rows[end].validation_rms))
end

open(joinpath(output, "tci_ladder.csv"), "w") do io
    println(io, "tci_tolerance,tci_time_s,tci_oracle_evaluations,density_tci_max_chi,estimated_particle_number,estimated_particle_number_error,validation_max_abs_error,validation_rms_error")
    for row in ladder_rows
        println(io, "$(row.tolerance),$(row.time_s),$(row.queries),$(row.max_chi),$(row.trace),$(row.trace_error),$(row.validation_maximum),$(row.validation_rms)")
    end
end
open(joinpath(output, "validation_density.csv"), "w") do io
    headers = ["qtt_index", "raw_probe_density"]
    append!(headers, ["tci_" * replace(string(tolerance), "." => "p", "-" => "m") for tolerance in tci_tolerances])
    println(io, join(headers, ','))
    for (position, index) in enumerate(validation_indices)
        row = [string(index), @sprintf("%.16g", validation_raw[position])]
        append!(row, [@sprintf("%.16g", validation_columns[tolerance][position]) for tolerance in tci_tolerances])
        println(io, join(row, ','))
    end
end

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file,
        "label" => spec.label, "matrix_dimension" => N, "moments" => moments,
        "probes" => probes, "probe_row_offset" => probe_row_offset,
        "maxdim" => maxdim, "mpo_cutoff" => effective_params.itensors_tol,
        "tci_tolerances" => tci_tolerances, "chemical_potential" => 0.0,
        "spectral_lower" => lower, "spectral_upper" => upper,
        "device_name" => CUDA.name(CUDA.device()),
        "device_total_memory_bytes" => CUDA.total_memory(),
        "total_propagation_time_s" => sum(row.time_s for row in propagation_rows),
        "host_tci_note" => "final GPU MPS states were copied to host individually; each TCI fit is host-resident",
        "validation_note" => "fixed-H diagnostic: samples compare each TCI fit to the same R-probe MPS oracle, not to ED or full vector KPM",
    ))
end
progress("hybrid GPU-propagated QTT density-oracle ladder complete: $output")
