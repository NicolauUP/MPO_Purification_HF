@testset "Public MPO solve compatibility wrapper" begin
    model = ChainModel(
        size=8,
        hopping=-0.7,
        interaction=0.0,
        potential=x -> 0.08 * Int(x),
        filling=0.5,
    )
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64)
    solver = SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1,
        maxiter=5, purification_maxiter=35,
        sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6,
        record_energy=false,
        stable_iterations=2,
        detect_two_cycles=true,
    )
    reference_params = legacy_parameters(model, representation, solver)
    reference_system = System(reference_params)
    reference_hamiltonian = [
        MatrixChecker(
            reference_system.H0, reference_system.sites, i, j,
            reference_system.bra_states, reference_system.ket_states,
        ) for i in 1:8, j in 1:8
    ]
    reference_spectrum = eigvals(Hermitian(reference_hamiltonian))
    bounds = (minimum(reference_spectrum) - 0.5, maximum(reference_spectrum) + 0.5)

    result = solve(
        model, representation, solver;
        spectral_bounds=bounds, verbose=:nothing, verify_spectral_bounds=true,
    )
    @test result isa SolveResult
    @test result.method == :mpo
    @test result.runtime.backend == :cpu
    @test result.runtime_report.execution_path == :cpu_end_to_end
    @test result.spectral_bounds == bounds
    @test result.target_particles == 4
    @test result.converged
    @test result.termination_reason == :converged
    @test result.diagnostics.converged == result.converged
    @test result.observables !== nothing
    @test isapprox(result.observables.particle_number, 4.0; atol=3e-3)
    @test occursin("SolveResult(method=mpo", sprint(show, result))

    sys = System(reference_params)
    @test run_scf!(
        sys, bounds...;
        verbose=:nothing, purification_method=:sp2, verify_spectral_bounds=true,
    )
    legacy_observables = observables(sys)
    @test result.converged == scf_diagnostics(sys).converged
    @test result.termination_reason == scf_diagnostics(sys).termination_reason
    @test length(result.diagnostics.history) == length(scf_diagnostics(sys).history)
    @test isapprox(
        result.observables.site_density, legacy_observables.site_density;
        atol=1e-12, rtol=1e-12,
    )
    @test isapprox(
        result.observables.particle_number, legacy_observables.particle_number;
        atol=1e-12, rtol=1e-12,
    )
    @test isapprox(
        result.observables.energy.total, legacy_observables.energy.total;
        atol=1e-12, rtol=1e-12,
    )

    pulay_system = System(reference_params)
    applied_mixing_methods = Symbol[]
    @test run_scf!(
        pulay_system, bounds...;
        verbose=:nothing, purification_method=:sp2,
        sp2_idempotency_tolerance=2e-4, sp2_trace_tolerance=4e-6,
        mixing_method=:pulay, pulay_history=4, pulay_warmup=3,
        verify_spectral_bounds=true,
        phase_callback=timing -> push!(applied_mixing_methods, timing.mixing_method),
    )
    @test applied_mixing_methods[1] == :direct
    @test :linear_warmup in applied_mixing_methods
    @test :pulay in applied_mixing_methods
    @test isapprox(observables(pulay_system).particle_number, 4.0; atol=3e-3)

    dense_result = solve(model, representation, solver; method=:dense)
    @test dense_result.method == :dense
    @test dense_result.spectral_bounds === nothing
    @test dense_result.converged
    @test dense_result.termination_reason == :converged
    @test isapprox(dense_result.observables.site_density, result.observables.site_density; atol=3e-3)
    @test isapprox(dense_result.observables.energy.total, result.observables.energy.total; atol=3e-3)

    mktempdir() do temporary_directory
        output = joinpath(temporary_directory, "public_result")
        @test write_result(result, output) == abspath(output)
        for filename in (
            "input.toml", "metadata.toml", "observables.toml", "scf_history.csv",
            "site_density.csv", "bond_order.csv",
        )
            @test isfile(joinpath(output, filename))
        end
        stored = read_result(output)
        @test stored isa StoredResult
        @test stored.metadata["format_version"] == 1
        @test stored.metadata["scf_converged"] == result.converged
        @test stored.metadata["runtime_active_backend"] == "cpu"
        @test stored.metadata["runtime_execution_path"] == "cpu_end_to_end"
        @test stored.input["function_serialization"] isa String
        @test stored.input["square_fock_method"] == "binary_carry"
        @test stored.input["sp2_idempotency_tolerance"] == 2e-4
        @test stored.input["sp2_relative_trace_tolerance"] == 1e-6
        @test stored.input["sp2_absolute_trace_tolerance"] == 4e-6
        @test !stored.input["record_energy"]
        @test stored.input["stable_iterations"] == 2
        @test stored.input["detect_two_cycles"]
        @test stored.input["mixing_method"] == "linear"
        @test stored.input["pulay_history"] == 4
        @test stored.input["pulay_warmup"] == 4
        @test stored.input["pulay_regularization"] == 1e-12
        @test stored.input["pulay_coefficient_limit"] == 8.0
        @test stored.input["pulay_step_limit"] == 20.0
        @test stored.observables["available"]
        @test isapprox(stored.site_density, result.observables.site_density; atol=1e-12)
        @test length(stored.bond_order) == length(result.observables.bond_order)
        @test isapprox(stored.bond_order[1].value, result.observables.bond_order[1]; atol=1e-12)
        @test length(stored.history) == length(result.diagnostics.history)
        @test_throws ArgumentError write_result(result, output)
    end

    mktempdir() do temporary_directory
        output = joinpath(temporary_directory, "public_dense_result")
        write_result(dense_result, output)
        stored = read_result(output)
        @test stored.metadata["method"] == "dense"
        @test stored.metadata["solver"] == "public_dense_hf_compatibility"
        @test stored.metadata["spectral_bounds_status"] == "not_used"
        @test !haskey(stored.metadata, "spectral_lower")
        @test stored.input["spectral_bounds_status"] == "not_used"
        @test !stored.input["purification_applied"]
    end

    @test_throws ArgumentError solve(model, representation, solver; spectral_bounds=nothing)
    @test_throws ArgumentError solve(model, representation, solver; method=:dense, spectral_bounds=bounds)
    @test_throws ArgumentError solve(
        model, representation, solver;
        spectral_bounds=bounds, purification_method=:palser_manolopoulos,
    )
end

@testset "Public dense square compatibility wrapper" begin
    model = SquareModel(
        size=(4, 4),
        hopping=(-0.6, -0.35),
        interaction=0.15,
        potential=(x, y) -> 0.11 * Int(x) + 0.07 * Int(y) + 0.013 * Int(x) * Int(y),
        filling=0.5,
    )
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64)
    solver = SCFSettings(
        purification=:sp2, mixing=0.4, tolerance=0.1,
        maxiter=40, purification_maxiter=35,
    )
    dense = solve(model, representation, solver; method=:dense)
    mpo = solve(
        model, representation, solver;
        method=:mpo, spectral_bounds=(-5.0, 5.0), verbose=:nothing,
    )
    @test dense.converged
    @test mpo.converged
    @test dense.spectral_bounds === nothing
    @test isapprox(dense.observables.site_density, mpo.observables.site_density; atol=4e-3)
    @test isapprox(dense.observables.horizontal_bond_order, mpo.observables.horizontal_bond_order; atol=4e-3)
    @test isapprox(dense.observables.vertical_bond_order, mpo.observables.vertical_bond_order; atol=4e-3)
    @test isapprox(dense.observables.energy.total, mpo.observables.energy.total; atol=4e-3)
    @test dense.observables.idempotency_residual < 1e-12

    dense_preflight = preflight(model, representation, solver; method=:dense)
    @test dense_preflight.runnable
    @test dense_preflight.spectral_bounds === nothing
    @test dense_preflight.spectral_bounds_status == :not_required_by_selected_purification
    dense_bounds_preflight = preflight(
        model, representation, solver; method=:dense, spectral_bounds=(-5.0, 5.0),
    )
    @test !dense_bounds_preflight.runnable
    @test any(issue -> issue.code == :spectral_bounds_not_used_by_dense, dense_bounds_preflight.issues)
end
