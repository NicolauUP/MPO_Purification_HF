#!/usr/bin/env julia

"""
Diagnose the square binary-carry Hartree construction at selected physical sites.

The smooth density is known analytically, while the direct MPO probes use the
actual density MPO. This distinguishes source-density approximation, directional
binary-carry component errors, and error introduced by field assembly.

Example (large cluster case):

    julia --startup-file=no --project=. \
      benchmark/extraction_scaling_2d/diagnose_tensorial_hartree.jl \
      --side-level 10 --output binary_carry_hartree_Lside10.txt

Repeat with a different `--cutoff` or `--maxdim` rather than interpreting a
single global MPO norm as a local accuracy guarantee.
"""

include(joinpath(@__DIR__, "extraction_scaling_2d.jl"))

function _diagnostic_arguments(arguments)
    configuration = (side_level=10, cutoff=1e-12, maxdim=128, output=nothing)
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--side-level" && (configuration = merge(configuration, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--cutoff" && (configuration = merge(configuration, (cutoff=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdim" && (configuration = merge(configuration, (maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    2 <= configuration.side_level <= 10 || throw(ArgumentError("--side-level must lie in 2:10"))
    isfinite(configuration.cutoff) && configuration.cutoff > 0 || throw(ArgumentError("--cutoff must be positive"))
    configuration.maxdim > 0 || throw(ArgumentError("--maxdim must be positive"))
    configuration
end

function _smooth_system_for_diagnostic(side_level::Int, cutoff::Float64, maxdim::Int)
    total_bits = 2side_level
    side = 2^side_level
    params = ParametersSquare(
        L=total_bits, t=(-0.6, -0.35), U=0.3, W=nothing, S=nothing,
        tci_tol=1e-10, itensors_tol=cutoff, itensors_maxdim=maxdim,
        density=0.5, purification_steps=50, scf_mixing=0.5, scf_tol=0.1,
        scf_max_iterations=5,
    )
    sys = System(params)
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            0.45 + 0.08 * cospi(x / (side - 1)) + 0.05 * sinpi(y / (side - 1))
        end, Float64, sys.sites, params.tci_tol,
    )
    horizontal = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            x < side - 1 ? 0.03 * sinpi((x + 2y) / (2side - 2)) : 0.0
        end, Float64, sys.sites, params.tci_tol,
    )
    vertical = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            y < side - 1 ? 0.025 * cospi((2x + y) / (2side - 2)) : 0.0
        end, Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L, T_U, T_D = sys.translations
    rho_right = apply(horizontal, T_R; cutoff=cutoff, maxdim=maxdim)
    rho_left = apply(T_L, ITensors.dag(horizontal); cutoff=cutoff, maxdim=maxdim)
    rho_up = apply(vertical, T_U; cutoff=cutoff, maxdim=maxdim)
    rho_down = apply(T_D, ITensors.dag(vertical); cutoff=cutoff, maxdim=maxdim)
    sys.ρ = +(diagonal, rho_right, rho_left, rho_up, rho_down; cutoff=cutoff, maxdim=maxdim)
    return sys
end

function _rho_diagonal(sys, site::Int)
    real(MatrixChecker(sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states))
end

function _field_diagonal(field::MPO, sys, site::Int)
    real(MatrixChecker(field, sys.sites, site, site, sys.bra_states, sys.ket_states))
end

function _analytic_density(x::Int, y::Int, side::Int)
    0.45 + 0.08 * cospi(x / (side - 1)) + 0.05 * sinpi(y / (side - 1))
end

function _probes(side::Int)
    unique([
        (0, 0), (1, 0), (0, 1),
        (side - 1, 0), (side - 2, 0),
        (0, side - 1), (0, side - 2),
        (side - 1, side - 1), (side - 2, side - 1), (side - 1, side - 2),
        (div(side, 2), div(side, 2)),
        (max(1, div(side, 3)), max(1, div(2side, 3))),
    ])
end

function _direct_hartree(sys, site::Int)
    sum(_rho_diagonal(sys, neighbour) for neighbour in values(square_neighbours(site, sys.params.L)) if !isnothing(neighbour)) * sys.params.U
end

function _write_report(io::IO, configuration)
    sys = _smooth_system_for_diagnostic(configuration.side_level, configuration.cutoff, configuration.maxdim)
    params = sys.params
    side = 2^configuration.side_level
    L = params.L

    density_tensors = MPO_MeanField._density_diagonal_qtt_tensors(sys.ρ, sys.sites)
    components = (
        right=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :right),
            sys.sites, params; symmetrize=false,
        ),
        left=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :left),
            sys.sites, params; symmetrize=false,
        ),
        up=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :up),
            sys.sites, params; symmetrize=false,
        ),
        down=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :down),
            sys.sites, params; symmetrize=false,
        ),
    )
    raw_hartree = +(params.U * components.right, params.U * components.left,
        params.U * components.up, params.U * components.down;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    final_hartree = +(0.5 * raw_hartree, 0.5 * ITensors.dag(raw_hartree);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    public_hartree = extract_hartree_mpo_binary_carry_square(sys)

    println(io, "Square binary-carry Hartree diagnostic")
    println(io, "=======================================")
    println(io, "side_level=$(configuration.side_level), L_total=$L, side=$side, N=$(2^L)")
    println(io, "density source=smooth synthetic; tci_tol=$(params.tci_tol), cutoff=$(params.itensors_tol), maxdim=$(params.itensors_maxdim)")
    println(io, "rho chi: max=$(maxlinkdim(sys.ρ)), mean=$(_mean_chi(sys.ρ))")
    println(io, "field chi: raw=$(maxlinkdim(raw_hartree)), final=$(maxlinkdim(final_hartree)), public=$(maxlinkdim(public_hartree))")
    println(io, "global relative raw-vs-final=$(_relative_difference(raw_hartree, final_hartree, params))")
    println(io, "global relative final-vs-public=$(_relative_difference(final_hartree, public_hartree, params))")
    println(io)
    println(io, "Each component is expected to equal the actual rho diagonal at the named neighbour; a missing neighbour contributes zero.")
    println(io, @sprintf("%-13s %13s %13s %13s %13s", "probe (x,y)", "direct VH", "raw error", "final error", "rho diag error"))

    largest_component_error = (error=-Inf, probe=(0, 0), component=:none)
    for (x, y) in _probes(side)
        site = square_lattice_index(x, y, L)
        direct = _direct_hartree(sys, site)
        raw_value = _field_diagonal(raw_hartree, sys, site)
        final_value = _field_diagonal(final_hartree, sys, site)
        rho_value = _rho_diagonal(sys, site)
        analytic = _analytic_density(x, y, side)
        println(io, @sprintf("(%4d,%4d) %13.6e %13.6e %13.6e %13.6e", x, y, direct,
            raw_value - direct, final_value - direct, rho_value - analytic))

        neighbours = square_neighbours(site, L)
        for component in (:right, :left, :up, :down)
            neighbour = getproperty(neighbours, component)
            expected = isnothing(neighbour) ? 0.0 : _rho_diagonal(sys, neighbour)
            observed = _field_diagonal(getproperty(components, component), sys, site)
            error = observed - expected
            abs(error) > largest_component_error.error &&
                (largest_component_error = (error=abs(error), probe=(x, y), component=component))
            println(io, @sprintf("  %-5s expected=% .9e observed=% .9e error=%+.3e%s", string(component), expected, observed, error,
                isnothing(neighbour) ? " (open boundary)" : ""))
        end
    end
    println(io)
    println(io, "Largest binary-carry component error: $(largest_component_error.error) at probe $(largest_component_error.probe), component $(largest_component_error.component)")
    println(io, "Interpretation:")
    println(io, "  * rho diag error isolates approximation in the synthetic density MPO.")
    println(io, "  * component errors isolate a directional carry or boundary error before Hartree assembly.")
    println(io, "  * raw error isolates the four-component MPO addition; final error additionally includes Hermitian symmetrization.")
    println(io, "  * Repeat with --cutoff 1e-14 and --maxdim 256 (one change at a time) to test truncation sensitivity.")
end

function print_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/diagnose_tensorial_hartree.jl [options]")
    println("Diagnoses binary-carry components despite the legacy filename.")
    println("  --side-level N   2:10, default 10")
    println("  --cutoff X       MPO cutoff, default 1e-12")
    println("  --maxdim N       MPO maximum bond dimension, default 128")
    println("  --output FILE    also write the report to FILE")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = _diagnostic_arguments(ARGS)
    if isnothing(configuration)
        print_usage()
    else
        buffer = IOBuffer()
        _write_report(buffer, configuration)
        report = String(take!(buffer))
        print(stdout, report)
        if !isnothing(configuration.output)
            open(configuration.output, "w") do io
                print(io, report)
            end
            println("Report written to $(abspath(configuration.output))")
        end
    end
end
