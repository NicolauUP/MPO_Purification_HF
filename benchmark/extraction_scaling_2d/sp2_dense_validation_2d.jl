#!/usr/bin/env julia

"""Dense-reference validation of bounded square-MPO SP2 cases.

This is a diagnosis tool, not an extraction benchmark or production SCF run.
It deliberately limits `side_level <= 4`, so the largest dense calculation is
the 256-by-256 (`16 x 16`) checkerboard Hamiltonian.
"""

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using SHA

_csv(value) = "\"" * replace(string(value), '\"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))

function dense_matrix(mpo::MPO, sys::System)
    N = 2^sys.params.L
    [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states) for i in 1:N, j in 1:N]
end

function parse_arguments(arguments)
    configuration = (output=nothing, side_level=4, itensors_tol=1e-14,
        itensors_maxdim=128, steps=50, padding=0.5)
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (configuration = merge(configuration, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (configuration = merge(configuration, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-maxdim" && (configuration = merge(configuration, (itensors_maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (configuration = merge(configuration, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--padding" && (configuration = merge(configuration, (padding=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(configuration.output) && throw(ArgumentError("--output DIRECTORY is required"))
    1 <= configuration.side_level <= 4 || throw(ArgumentError(
        "dense validation is intentionally limited to side levels 1:4; got $(configuration.side_level)",
    ))
    isfinite(configuration.itensors_tol) && configuration.itensors_tol > 0 || throw(ArgumentError("itensors tolerance must be positive and finite"))
    configuration.itensors_maxdim > 0 || throw(ArgumentError("itensors maxdim must be positive"))
    configuration.steps > 0 || throw(ArgumentError("steps must be positive"))
    isfinite(configuration.padding) && configuration.padding > 0 || throw(ArgumentError("padding must be positive and finite"))
    configuration
end

function _parameters(configuration)
    total_bits = 2configuration.side_level
    ParametersSquare(
        L=total_bits, t=(-0.6, -0.35), U=0.0,
        W=(x, y) -> iseven(Int(x) + Int(y)) ? 0.6 : -0.6,
        S=nothing, tci_tol=1e-10,
        itensors_tol=configuration.itensors_tol,
        itensors_maxdim=configuration.itensors_maxdim,
        density=0.5, purification_steps=configuration.steps,
        scf_mixing=0.5, scf_tol=0.1, scf_max_iterations=5,
    )
end

function _write_metadata(output, configuration)
    project = Base.active_project()
    open(joinpath(output, "metadata.txt"), "w") do io
        println(io, "started_at_utc=", Dates.now(Dates.UTC))
        println(io, "active_project=", project)
        println(io, "project_sha256=", bytes2hex(sha256(read(project))))
        println(io, "julia_version=", VERSION)
        println(io, "cpu_name=", Sys.CPU_NAME)
        println(io, "side_level=", configuration.side_level)
        println(io, "itensors_tol=", configuration.itensors_tol)
        println(io, "itensors_maxdim=", configuration.itensors_maxdim)
        println(io, "purification_steps=", configuration.steps)
        println(io, "spectral_padding=", configuration.padding)
        println(io, "potential=checkerboard_plus_minus_0.6")
    end
end

function run_validation(; output, side_level, itensors_tol, itensors_maxdim, steps, padding)
    output = abspath(output)
    mkpath(output)
    configuration = (; output, side_level, itensors_tol, itensors_maxdim, steps, padding)
    _write_metadata(output, configuration)
    try
        params = _parameters(configuration)
        sys = System(params)
        H = dense_matrix(sys.H0, sys)
        hermitian_H = Hermitian((H + adjoint(H)) / 2)
        eigen_H = eigen(hermitian_H)
        N = size(H, 1)
        Ne = round(Int, N * params.density)
        fermi_gap = eigen_H.values[Ne + 1] - eigen_H.values[Ne]
        bounds = (minimum(eigen_H.values) - padding, maximum(eigen_H.values) + padding)
        exact = eigen_H.vectors[:, 1:Ne] * adjoint(eigen_H.vectors[:, 1:Ne])

        # The public exact-bound guard is deliberately limited to N <= 16.
        # Here the dense spectrum above is the independent validation, so do
        # not request that guard again during MPO construction.
        rho0 = construct_rho_0(sys, params, bounds...; method=:sp2, verify_spectral_bounds=false)
        scaled_dense = (bounds[2] * I - H) / (bounds[2] - bounds[1])
        initial_scaling_error = opnorm(dense_matrix(rho0, sys) - scaled_dense)

        history_path = joinpath(output, "sp2_history.txt")
        result = open(history_path, "w") do history
            perform_purification(
                rho0, params;
                method=:sp2,
                verbose=1,
                io=history,
                overwrite_progress=false,
                spectral_bounds=bounds,
                spectral_bounds_validation=:dense_reference_external_to_guard,
            )
        end
        rho = dense_matrix(result.rho, sys)
        rho_eigenvalues = eigvals(Hermitian((rho + adjoint(rho)) / 2))
        values = (
            "checkerboard_cdw", side_level, params.L, N, Ne,
            itensors_tol, itensors_maxdim, steps, padding,
            result.converged, result.termination_reason, result.iterations,
            result.final_bond_dimension, result.work.max_bond_dimension,
            minimum(eigen_H.values), maximum(eigen_H.values), fermi_gap,
            initial_scaling_error,
            real(tr(rho)), abs(real(tr(rho)) - Ne),
            opnorm(rho - rho * rho), opnorm(rho - exact),
            opnorm(H * rho - rho * H), opnorm(rho - adjoint(rho)),
            minimum(rho_eigenvalues), maximum(rho_eigenvalues),
        )
        open(joinpath(output, "summary.csv"), "w") do summary
            _row(summary, (
                "source", "L_side", "L_total", "N", "Ne", "itensors_tol", "itensors_maxdim", "purification_steps", "spectral_padding",
                "converged", "termination_reason", "iterations", "final_chi", "max_chi",
                "lambda_min", "lambda_max", "fermi_gap", "initial_scaling_operator_error",
                "trace", "trace_error", "idempotency_operator_error", "projector_operator_error",
                "commutator_operator_error", "hermiticity_operator_error", "rho_eigenvalue_min", "rho_eigenvalue_max",
            ))
            _row(summary, values)
        end
        println("Dense-reference SP2 validation written to $output")
    catch error
        open(joinpath(output, "error.txt"), "w") do io
            showerror(io, error, catch_backtrace())
            println(io)
        end
        rethrow()
    end
end

function print_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/sp2_dense_validation_2d.jl --output DIRECTORY [options]")
    println("  --side-level 4       (restricted to 1:4; level 4 is 16 x 16, N=256)")
    println("  --itensors-tol 1e-14 --itensors-maxdim 128 --steps 50 --padding 0.5")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_arguments(ARGS)
    isnothing(configuration) ? print_usage() : run_validation(; configuration...)
end
