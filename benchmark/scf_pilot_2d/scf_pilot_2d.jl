#!/usr/bin/env julia

"""Bounded single-case 2D square Hartree--Fock SCF pilot.

This is intentionally a *single-run* diagnostic, not a parameter sweep. It
writes sufficient numerical history and provenance to assess whether a chosen
MPO cap, truncation tolerance, and initial seed support a trustworthy later
campaign. To avoid accidental production-scale runs, `side_level` is capped at
six (`side = 64`, `N = 4096`).
"""

using Dates
using ITensors, ITensorMPS
using MPO_MeanField
using Printf
using SHA

_csv(value) = "\"" * replace(string(value), '\"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))
_toml_string(value) = "\"" * replace(string(value), '\"' => "\\\"") * "\""
_optional(value) = isnothing(value) ? "" : string(value)

function parse_arguments(arguments)
    config = (
        output=nothing, side_level=6, tx=-0.6, ty=-0.35, U=0.3,
        density=0.5, potential=:checkerboard, potential_amplitude=0.6,
        seed=:uniform, seed_amplitude=0.05, tci_tol=1e-10,
        itensors_tol=1e-14, maxdim=256, steps=50, padding=0.5,
        scf_mixing=0.5, scf_tol=0.1, scf_max_iterations=30,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (config = merge(config, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (config = merge(config, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--tx" && (config = merge(config, (tx=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--ty" && (config = merge(config, (ty=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--U" && (config = merge(config, (U=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--density" && (config = merge(config, (density=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--potential" && (config = merge(config, (potential=Symbol(arguments[index + 1]),)); index += 2; continue)
        argument == "--potential-amplitude" && (config = merge(config, (potential_amplitude=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--seed" && (config = merge(config, (seed=Symbol(arguments[index + 1]),)); index += 2; continue)
        argument == "--seed-amplitude" && (config = merge(config, (seed_amplitude=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--tci-tol" && (config = merge(config, (tci_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (config = merge(config, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdim" && (config = merge(config, (maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (config = merge(config, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--padding" && (config = merge(config, (padding=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-mixing" && (config = merge(config, (scf_mixing=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-tol" && (config = merge(config, (scf_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-max-iterations" && (config = merge(config, (scf_max_iterations=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(config.output) && throw(ArgumentError("--output DIRECTORY is required"))
    1 <= config.side_level <= 6 || throw(ArgumentError("side level must lie in 1:6 for this bounded SCF pilot"))
    config.potential in (:none, :checkerboard) || throw(ArgumentError("--potential must be none or checkerboard"))
    config.seed in (:uniform, :checkerboard_plus, :checkerboard_minus) || throw(ArgumentError("--seed must be uniform, checkerboard_plus, or checkerboard_minus"))
    all(isfinite, (config.tx, config.ty, config.U, config.density, config.potential_amplitude, config.seed_amplitude, config.tci_tol, config.itensors_tol, config.padding, config.scf_mixing, config.scf_tol)) || throw(ArgumentError("all floating-point inputs must be finite"))
    config.potential_amplitude >= 0 || throw(ArgumentError("potential amplitude must be nonnegative"))
    config.seed_amplitude >= 0 || throw(ArgumentError("seed amplitude must be nonnegative"))
    config.tci_tol > 0 && config.itensors_tol > 0 || throw(ArgumentError("TCI and ITensor tolerances must be positive"))
    config.maxdim > 0 && config.steps > 0 && config.scf_max_iterations > 0 || throw(ArgumentError("maxdim, steps, and SCF iterations must be positive"))
    0 < config.density < 1 || throw(ArgumentError("density must lie strictly between zero and one"))
    0 <= config.scf_mixing <= 1 && config.scf_tol > 0 && config.padding >= 0 || throw(ArgumentError("invalid SCF mixing, tolerance, or padding"))
    return config
end

function _checkerboard(amplitude, sign=1.0)
    (x, y) -> sign * (iseven(Int(x) + Int(y)) ? amplitude : -amplitude)
end

function _parameters(config)
    potential = config.potential == :none ? nothing : _checkerboard(config.potential_amplitude)
    seed = if config.seed == :uniform
        nothing
    elseif config.seed == :checkerboard_plus
        _checkerboard(config.seed_amplitude)
    else
        _checkerboard(config.seed_amplitude, -1.0)
    end
    return ParametersSquare(
        L=2config.side_level, t=(config.tx, config.ty), U=config.U,
        W=potential, S=seed, tci_tol=config.tci_tol,
        itensors_tol=config.itensors_tol, itensors_maxdim=config.maxdim,
        density=config.density, purification_steps=config.steps,
        scf_mixing=config.scf_mixing, scf_tol=config.scf_tol,
        scf_max_iterations=config.scf_max_iterations,
    )
end

function _write_metadata(path, config, params, bounds)
    open(path, "w") do io
        println(io, "format_version = 1")
        println(io, "created_at_utc = ", _toml_string(Dates.format(now(UTC), DateFormat("yyyy-mm-ddTHH:MM:SS"))))
        println(io, "runner = ", _toml_string("benchmark/scf_pilot_2d/scf_pilot_2d.jl"))
        println(io, "project_toml_sha256 = ", _toml_string(bytes2hex(sha256(read("Project.toml")))))
        println(io, "side_level = $(config.side_level)")
        println(io, "L_total = $(params.L)")
        println(io, "N = $(2^params.L)")
        println(io, "target_particles = $(round(Int, config.density * 2^params.L))")
        for name in (:tx, :ty, :U, :density, :potential_amplitude, :seed_amplitude, :tci_tol, :itensors_tol, :padding, :scf_mixing, :scf_tol)
            println(io, "$(name) = $(getfield(config, name))")
        end
        println(io, "potential = ", _toml_string(config.potential))
        println(io, "seed = ", _toml_string(config.seed))
        println(io, "itensors_maxdim = $(config.maxdim)")
        println(io, "purification_steps = $(config.steps)")
        println(io, "scf_max_iterations = $(config.scf_max_iterations)")
        println(io, "spectral_lower = $(bounds[1])")
        println(io, "spectral_upper = $(bounds[2])")
    end
end

function _write_history(path, diagnostics)
    open(path, "w") do io
        _row(io, ("iteration", "trace", "vh_residual", "vf_residual", "rho_residual", "commutator_residual", "two_cycle_residual", "purification_converged", "purification_termination_reason", "purification_iterations", "rho_bond_dimension", "hartree_bond_dimension", "fock_bond_dimension", "effective_hamiltonian_bond_dimension", "energy_total"))
        for record in diagnostics.history
            _row(io, (record.iteration, record.trace, record.vh_residual, record.vf_residual, record.rho_residual, record.commutator_residual, record.two_cycle_residual, record.purification_converged, record.purification_termination_reason, record.purification_iterations, _optional(record.rho_bond_dimension), _optional(record.hartree_bond_dimension), _optional(record.fock_bond_dimension), _optional(record.effective_hamiltonian_bond_dimension), _optional(record.energy_total)))
        end
    end
end

function _probe_sites(side, L)
    coordinates = unique(((0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1), (div(side, 2), div(side, 2))))
    [(x, y, MPO_MeanField.square_lattice_index(x, y, L)) for (x, y) in coordinates]
end

function _write_observables(path, sys, converged, bounds)
    params = sys.params
    N = 2^params.L
    side = 2^div(params.L, 2)
    density_sum = 0.0
    density_min, density_max = Inf, -Inf
    checkerboard_sum = 0.0
    for site in 1:N
        value = real(MatrixChecker(sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states))
        x, y = MPO_MeanField.square_lattice_decoder(site - 1, params.L)
        density_sum += value
        density_min = min(density_min, value)
        density_max = max(density_max, value)
        checkerboard_sum += iseven(x + y) ? value : -value
    end
    energy = nearest_neighbor_hf_energy_square(sys)
    effective_hamiltonian = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    idempotency = MPO_MeanField._mpo_relative_change(
        apply(sys.ρ, sys.ρ; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim), sys.ρ, params,
    )
    stationarity = MPO_MeanField._scf_commutator_residual(effective_hamiltonian, sys.ρ, params)
    open(path, "w") do io
        println(io, "scf_converged = $(converged)")
        println(io, "scf_termination_reason = ", _toml_string(scf_diagnostics(sys).termination_reason))
        println(io, "spectral_lower = $(bounds[1])")
        println(io, "spectral_upper = $(bounds[2])")
        println(io, "particle_number = $(density_sum)")
        println(io, "particle_number_error = $(density_sum - round(Int, params.density * N))")
        println(io, "density_min = $(density_min)")
        println(io, "density_max = $(density_max)")
        println(io, "checkerboard_density_order = $(checkerboard_sum / N)")
        println(io, "idempotency_residual = $(idempotency)")
        println(io, "stationarity_residual = $(stationarity)")
        println(io, "rho_bond_dimension = $(maxlinkdim(sys.ρ))")
        println(io, "hartree_bond_dimension = $(maxlinkdim(sys.VH))")
        println(io, "fock_bond_dimension = $(maxlinkdim(sys.VF))")
        println(io, "effective_hamiltonian_bond_dimension = $(maxlinkdim(effective_hamiltonian))")
        for name in (:kinetic, :hartree, :fock, :interaction, :total)
            println(io, "energy_$(name) = $(getfield(energy, name))")
        end
        for (x, y, site) in _probe_sites(side, params.L)
            println(io, "\n[[density_probes]]")
            println(io, "x = $x\ny = $y\nsite = $site")
            println(io, "density = ", real(MatrixChecker(sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states)))
            println(io, "stored_hartree = ", real(MatrixChecker(sys.VH, sys.sites, site, site, sys.bra_states, sys.ket_states)))
        end
    end
end

function run_pilot(config)
    output = abspath(config.output)
    mkpath(output)
    params = _parameters(config)
    potential_bounds = config.potential == :none ? (0.0, 0.0) : (-config.potential_amplitude, config.potential_amplitude)
    bounds = square_scf_spectral_bounds(params; potential_bounds=potential_bounds, margin=config.padding)
    _write_metadata(joinpath(output, "metadata.toml"), config, params, bounds)
    sys = System(params)
    converged = open(joinpath(output, "progress.txt"), "w") do progress
        run_scf!(sys, bounds...; purification_method=:sp2, verbose=:all, io=progress,
            overwrite_progress=false, record_energy=true)
    end
    _write_history(joinpath(output, "scf_history.csv"), scf_diagnostics(sys))
    _write_observables(joinpath(output, "observables.toml"), sys, converged, bounds)
    println("SCF pilot completed: converged=$converged termination=$(scf_diagnostics(sys).termination_reason)")
    println("Output directory: $output")
    return converged
end

function main(arguments)
    config = parse_arguments(arguments)
    if isnothing(config)
        println("Usage: scf_pilot_2d.jl --output DIRECTORY [options]")
        println("Options: --side-level 1..6 --potential none|checkerboard --seed uniform|checkerboard_plus|checkerboard_minus")
        return
    end
    run_pilot(config)
end

main(ARGS)
