#!/usr/bin/env julia

"""Measure extra SP2 polynomials after the current stopping rule returns.

This diagnostic does not alter production purification. It answers whether the
current `idempotency_threshold` stops before the dense-projector error reaches
the desired scientific accuracy on a bounded checkerboard square case.
"""

include(joinpath(@__DIR__, "sp2_dense_validation_2d.jl"))

function parse_refinement_arguments(arguments)
    configuration = (output=nothing, side_level=4, itensors_tol=1e-14,
        itensors_maxdim=256, steps=50, padding=0.5, extra_iterations=8)
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (configuration = merge(configuration, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (configuration = merge(configuration, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-maxdim" && (configuration = merge(configuration, (itensors_maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (configuration = merge(configuration, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--padding" && (configuration = merge(configuration, (padding=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--extra-iterations" && (configuration = merge(configuration, (extra_iterations=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(configuration.output) && throw(ArgumentError("--output DIRECTORY is required"))
    1 <= configuration.side_level <= 4 || throw(ArgumentError("dense validation is intentionally limited to side levels 1:4"))
    isfinite(configuration.itensors_tol) && configuration.itensors_tol > 0 || throw(ArgumentError("itensors tolerance must be positive and finite"))
    configuration.itensors_maxdim > 0 || throw(ArgumentError("itensors maxdim must be positive"))
    configuration.steps > 0 || throw(ArgumentError("steps must be positive"))
    0 <= configuration.extra_iterations <= 16 || throw(ArgumentError("extra iterations must lie in 0:16"))
    configuration
end

function _metrics(rho::MPO, H, exact, sys)
    dense_rho = dense_matrix(rho, sys)
    eigenvalues = eigvals(Hermitian((dense_rho + adjoint(dense_rho)) / 2))
    return (
        trace=real(tr(dense_rho)),
        idempotency=opnorm(dense_rho - dense_rho * dense_rho),
        projector=opnorm(dense_rho - exact),
        commutator=opnorm(H * dense_rho - dense_rho * H),
        hermiticity=opnorm(dense_rho - adjoint(dense_rho)),
        eigenvalue_min=minimum(eigenvalues),
        eigenvalue_max=maximum(eigenvalues),
        chi=maxlinkdim(rho),
    )
end

function _write_row(io, stage, extra_step, branch, metrics, Ne)
    _row(io, (
        stage, extra_step, branch, metrics.chi,
        metrics.trace, abs(metrics.trace - Ne), metrics.idempotency,
        metrics.projector, metrics.commutator, metrics.hermiticity,
        metrics.eigenvalue_min, metrics.eigenvalue_max,
    ))
end

function run_refinement(; output, side_level, itensors_tol, itensors_maxdim, steps, padding, extra_iterations)
    output = abspath(output)
    mkpath(output)
    configuration = (; output, side_level, itensors_tol, itensors_maxdim, steps, padding)
    _write_metadata(output, configuration)
    open(joinpath(output, "metadata.txt"), "a") do io
        println(io, "extra_iterations=", extra_iterations)
        println(io, "note=extra steps are diagnostic only; production SP2 is unchanged")
    end
    try
        params = _parameters(configuration)
        sys = System(params)
        H = dense_matrix(sys.H0, sys)
        eigen_H = eigen(Hermitian((H + adjoint(H)) / 2))
        N = size(H, 1)
        Ne = round(Int, N * params.density)
        bounds = (minimum(eigen_H.values) - padding, maximum(eigen_H.values) + padding)
        exact = eigen_H.vectors[:, 1:Ne] * adjoint(eigen_H.vectors[:, 1:Ne])
        rho0 = construct_rho_0(sys, params, bounds...; method=:sp2, verify_spectral_bounds=false)

        result = open(joinpath(output, "sp2_history_until_stop.txt"), "w") do history
            perform_purification(
                rho0, params;
                method=:sp2, verbose=1, io=history, overwrite_progress=false,
                spectral_bounds=bounds,
                spectral_bounds_validation=:dense_reference_external_to_guard,
            )
        end
        rho = result.rho
        trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)

        open(joinpath(output, "refinement.csv"), "w") do summary
            _row(summary, (
                "stage", "extra_step", "branch", "field_max_chi", "trace", "trace_error",
                "idempotency_operator_error", "projector_operator_error", "commutator_operator_error",
                "hermiticity_operator_error", "rho_eigenvalue_min", "rho_eigenvalue_max",
            ))
            _write_row(summary, "normal_sp2_stop", 0, "not_applied", _metrics(rho, H, exact, sys), Ne)

            for extra_step in 1:extra_iterations
                rho_squared = apply(rho, rho;
                    cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
                )
                branch = MPO_MeanField._sp2_branch(
                    real(tr(rho)), real(tr(rho_squared)), Ne, trace_tolerance,
                )
                rho = branch == :square ? rho_squared : +(
                    2.0 * rho, -rho_squared;
                    cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
                )
                _write_row(summary, "forced_refinement", extra_step, branch, _metrics(rho, H, exact, sys), Ne)
                flush(summary)
            end
        end
        println("SP2 refinement diagnostics written to $output")
    catch error
        open(joinpath(output, "error.txt"), "w") do io
            showerror(io, error, catch_backtrace())
            println(io)
        end
        rethrow()
    end
end

function print_refinement_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/sp2_refinement_validation_2d.jl --output DIRECTORY [options]")
    println("  --side-level 4 --itensors-tol 1e-14 --itensors-maxdim 256")
    println("  --steps 50 --padding 0.5 --extra-iterations 8")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_refinement_arguments(ARGS)
    isnothing(configuration) ? print_refinement_usage() : run_refinement(; configuration...)
end
