#!/usr/bin/env julia

"""Fixed-H feasibility diagnostic for QTT KPM MPS Hadamard probes.

Usage:
    julia --project=. benchmark/qtt_kpm_mps_probe/diagnose_qtt_kpm_mps_probe.jl \
      CAMPAIGN.jl TASK OUTPUT [--moments M] [--probes R] [--maxdim CHI] [--cutoff EPS]

The campaign must define `campaign::CampaignSpec` and its selected case must
be an equal square. This is intentionally *not* an SCF driver. It compares:

  1. MPO--MPS Chebyshev propagation of rank-one Walsh--Hadamard probes;
  2. the identical dense-vector recurrence in QTT site order; and
  3. the exact dense fixed-H projector.

Thus `mps_vs_vector_*` measures MPO--MPS truncation, while
`vector_vs_exact_*` combines Jackson-polynomial and finite-probe errors.
"""

using Dates
using DelimitedFiles
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS

# Campaign callables must enter the world before `main` is compiled: the
# legacy MPO construction invokes hopping/seed closures from deep inside TCI.
if !isempty(ARGS) && isfile(ARGS[1])
    Base.include(Main, abspath(ARGS[1]))
end

function usage()
    println(stderr, "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT [--moments M] [--probes R] [--maxdim CHI] [--cutoff EPS]")
end

function parse_options(arguments)
    length(arguments) >= 3 || (usage(); error("missing required arguments"))
    options = Dict{String,String}("moments" => "256", "probes" => "16", "maxdim" => "128", "cutoff" => "1e-10")
    position = 4
    while position <= length(arguments)
        key = arguments[position]
        startswith(key, "--") || error("expected option beginning with --, got $key")
        position < length(arguments) || error("missing value for $key")
        normalized = key[3:end]
        haskey(options, normalized) || error("unknown option $key")
        options[normalized] = arguments[position + 1]
        position += 2
    end
    return (
        campaign=abspath(arguments[1]), task=parse(Int, arguments[2]), output=abspath(arguments[3]),
        moments=parse(Int, options["moments"]), probes=parse(Int, options["probes"]),
        maxdim=parse(Int, options["maxdim"]), cutoff=parse(Float64, options["cutoff"]),
    )
end

relative_error(value, reference) = norm(value - reference) / max(norm(reference), sqrt(eps(Float64)))
rms_error(value, reference) = norm(value - reference) / sqrt(length(value))

function qtt_row_permutation(levels::Int)
    side = 1 << (levels ÷ 2)
    return [begin
        x, y = square_lattice_decoder(index, levels)
        x + side * y + 1
    end for index in 0:((1 << levels) - 1)]
end

function local_estimates(PZ, Z, bonds)
    R = size(Z, 2)
    density = vec(sum(PZ .* Z; dims=2)) ./ R
    order = Vector{Float64}(undef, length(bonds))
    for k in eachindex(bonds)
        i, j = bonds[k]
        order[k] = (dot(@view(PZ[i, :]), @view(Z[j, :])) + dot(@view(PZ[j, :]), @view(Z[i, :]))) / (2R)
    end
    trace = sum(PZ .* Z) / R
    trace_squared = sum(abs2, PZ) / R
    return (; density, order, trace, trace_squared,
        idempotency=abs(trace - trace_squared) / max(abs(trace), sqrt(eps(Float64))))
end

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function main(arguments)
    config = parse_options(arguments)
    isfile(config.campaign) || error("campaign file does not exist: $(config.campaign)")
    ispath(config.output) && error("refusing to overwrite existing output directory: $(config.output)")
    config.moments >= 1 || error("moments must be positive")
    config.probes >= 1 || error("probes must be positive")
    config.maxdim >= 1 || error("maxdim must be positive")
    config.cutoff > 0 || error("cutoff must be positive")
    isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
    case = campaign_case(Main.campaign, config.task)
    model = case.model
    model isa SquareModel || error("this QTT MPO diagnostic currently accepts SquareModel only")
    nx, ny = model.size
    nx == ny || error("the legacy QTT MPO construction supports equal squares only")
    N = nx * ny
    ispow2(N) || error("QTT probe diagnostic requires a power-of-two site count")
    config.probes <= N || error("probes must not exceed N=$N")
    levels = qtt_levels(model)
    iseven(levels) || error("square QTT level count must be even")

    mkpath(config.output)
    started_at = string(now())
    representation = QTTSettings(
        encoding=:interleaved, tci_tol=case.representation.tci_tol,
        cutoff=config.cutoff, maxdim=config.maxdim,
    )
    parameters = legacy_parameters(model, representation, case.solver)
    system = System(parameters)
    H_effective = +(system.H0, system.VH, system.VF; cutoff=config.cutoff, maxdim=config.maxdim)

    # The sparse KPM helper supplies a direct, independent fixed-H reference.
    # Reordering it to the QTT basis makes its Hadamard rows exactly rank one.
    data = MPO_MeanField._kpm_data(model)
    H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
    qtt_to_row = qtt_row_permutation(levels)
    row_to_qtt = invperm(qtt_to_row)
    H_qtt = Matrix(H_row[qtt_to_row, qtt_to_row])
    spectrum = eigen(Symmetric(H_qtt))
    Ne = round(Int, model.filling * N)
    0 < Ne < N || error("fixed-H dense reference requires 0 < Ne < N")
    lower, upper = spectrum.values[1], spectrum.values[end]
    half_width = (upper - lower) / 2
    half_width > 0 || error("Hamiltonian has a degenerate spectral interval")
    center = (upper + lower) / 2
    padding = 0.05
    scale = half_width * (1 + padding)
    scaled_mu = ((spectrum.values[Ne] + spectrum.values[Ne + 1]) / 2 - center) / scale
    coefficients = MPO_MeanField._kpm_coefficients(config.moments, scaled_mu)
    H_scaled_dense = (H_qtt - center * Matrix{Float64}(I, N, N)) / scale
    H_scaled_mpo = +(
        H_effective / scale,
        (-center / scale) * Identity_MPO(system.sites);
        cutoff=config.cutoff,
        maxdim=config.maxdim,
    )

    # The QTT order is an interleaving of y/x bits for an equal square, hence
    # the coordinate Hadamard color code equals the QTT basis index.
    bonds_qtt = [(row_to_qtt[i], row_to_qtt[j]) for (i, j, _) in data.bonds]
    probes = hcat([
        MPO_MeanField._qtt_hadamard_probe_vector(levels, row)
        for row in 0:(config.probes - 1)
    ]...)
    vector_output = MPO_MeanField._kpm_apply(H_scaled_dense, probes, coefficients)
    mps_output = Matrix{Float64}(undef, N, config.probes)
    trajectory_rows = Vector{Vector{String}}()
    mps_time = 0.0
    for column in 1:config.probes
        probe = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, column - 1)
        elapsed = @elapsed result = MPO_MeanField._qtt_mps_chebyshev_apply(
            H_scaled_mpo, probe, coefficients; cutoff=config.cutoff, maxdim=config.maxdim,
        )
        mps_time += elapsed
        mps_output[:, column] .= MPO_MeanField._qtt_mps_amplitudes(result.state, system.sites)
        for record in result.trajectory
            push!(trajectory_rows, [
                string(column), string(record.order), string(record.state_max_chi),
                @sprintf("%.16g", record.state_mean_chi), string(record.result_max_chi),
                @sprintf("%.16g", record.result_mean_chi), @sprintf("%.16g", elapsed),
            ])
        end
    end
    mps_local = local_estimates(mps_output, probes, bonds_qtt)
    vector_local = local_estimates(vector_output, probes, bonds_qtt)
    exact_projector = spectrum.vectors[:, 1:Ne] * spectrum.vectors[:, 1:Ne]'
    exact_density = real.(diag(exact_projector))
    exact_order = [real((exact_projector[i, j] + exact_projector[j, i]) / 2) for (i, j) in bonds_qtt]

    # QTT/TCI is applied only after local fields are estimated: it is not used
    # in the recurrence and therefore cannot hide MPO--MPS propagation error.
    _, _, density_mps = MPO_MeanField.Quantics_TCI(
        index -> mps_local.density[Int(index) + 1], Float64, system.sites, config.cutoff,
    )
    density_tci = MPO_MeanField._qtt_mps_amplitudes(density_mps, system.sites)

    write_csv(joinpath(config.output, "trajectory.csv"),
        ["probe", "order", "chebyshev_max_chi", "chebyshev_mean_chi", "accumulator_max_chi", "accumulator_mean_chi", "probe_time_s"],
        trajectory_rows,
    )
    write_csv(joinpath(config.output, "density.csv"), ["qtt_index", "x", "y", "mps", "vector_kpm", "exact"], [
        [string(index), string(square_lattice_decoder(index, levels)[1]), string(square_lattice_decoder(index, levels)[2]),
         @sprintf("%.16g", mps_local.density[index + 1]), @sprintf("%.16g", vector_local.density[index + 1]),
         @sprintf("%.16g", exact_density[index + 1])]
        for index in 0:(N - 1)
    ])
    write_csv(joinpath(config.output, "bond_order.csv"), ["bond", "qtt_i", "qtt_j", "mps", "vector_kpm", "exact"], [
        [string(k), string(i - 1), string(j - 1), @sprintf("%.16g", mps_local.order[k]),
         @sprintf("%.16g", vector_local.order[k]), @sprintf("%.16g", exact_order[k])]
        for (k, (i, j)) in enumerate(bonds_qtt)
    ])
    max_chi = maximum(parse.(Int, getindex.(trajectory_rows, 3)))
    max_accumulator_chi = maximum(parse.(Int, getindex.(trajectory_rows, 5)))
    summary = Dict(
        "created_at" => started_at, "campaign" => config.campaign, "label" => case.label,
        "matrix_dimension" => N, "levels" => levels, "moments" => config.moments,
        "probes" => config.probes, "maxdim" => config.maxdim, "cutoff" => config.cutoff,
        "spectral_lower_exact" => lower, "spectral_upper_exact" => upper,
        "scaled_chemical_potential" => scaled_mu, "target_particles" => Ne,
        "mps_total_time_s" => mps_time, "mps_mean_time_per_probe_s" => mps_time / config.probes,
        "mps_chebyshev_max_chi" => max_chi, "mps_accumulator_max_chi" => max_accumulator_chi,
        "mps_vs_vector_state_relative_error" => relative_error(mps_output, vector_output),
        "mps_vs_vector_density_max_abs_error" => maximum(abs.(mps_local.density - vector_local.density)),
        "mps_vs_vector_density_rms_error" => rms_error(mps_local.density, vector_local.density),
        "mps_vs_vector_bond_max_abs_error" => maximum(abs.(mps_local.order - vector_local.order)),
        "mps_vs_vector_bond_rms_error" => rms_error(mps_local.order, vector_local.order),
        "mps_trace" => mps_local.trace, "mps_trace_idempotency_defect" => mps_local.idempotency,
        "vector_trace" => vector_local.trace, "vector_trace_idempotency_defect" => vector_local.idempotency,
        "vector_vs_exact_density_max_abs_error" => maximum(abs.(vector_local.density - exact_density)),
        "vector_vs_exact_density_rms_error" => rms_error(vector_local.density, exact_density),
        "vector_vs_exact_bond_max_abs_error" => maximum(abs.(vector_local.order - exact_order)),
        "vector_vs_exact_bond_rms_error" => rms_error(vector_local.order, exact_order),
        "mps_vs_exact_density_max_abs_error" => maximum(abs.(mps_local.density - exact_density)),
        "mps_vs_exact_density_rms_error" => rms_error(mps_local.density, exact_density),
        "density_tci_max_chi" => maxlinkdim(density_mps),
        "density_tci_relative_fit_error" => relative_error(density_tci, mps_local.density),
    )
    open(joinpath(config.output, "summary.toml"), "w") do io
        TOML.print(io, summary)
    end
    println("QTT KPM MPS-probe diagnostic complete: $(config.output)")
    println("  MPO-MPS vs vector state relative error = $(summary["mps_vs_vector_state_relative_error"])")
    println("  MPO-MPS max χ: Chebyshev=$(max_chi), accumulator=$(max_accumulator_chi)")
    println("  vector KPM vs exact density RMS (polynomial + R probes) = $(summary["vector_vs_exact_density_rms_error"])")
end

main(ARGS)
