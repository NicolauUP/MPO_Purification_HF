@testset "P4.1 SCF scalar diagnostics history" begin
    params = parameters_1d(
        L=2,
        t=-0.7,
        W=x -> (0.0, 0.13, -0.08, 0.21)[Int(x) + 1],
        U=0.2,
        S=nothing,
        purification_steps=35,
        scf_mixing=0.4,
        scf_tol=0.1,
        scf_max_iterations=20,
    )
    sys = System(params)
    @test scf_diagnostics(sys).termination_reason == :not_run
    @test run_scf!(sys, -5.0, 5.0;
        purification_method=:sp2,
        record_energy=true,
        verbose=:nothing,
    )

    diagnostics = scf_diagnostics(sys)
    @test diagnostics.converged
    @test diagnostics.termination_reason == :converged
    @test length(diagnostics.history) >= 2
    @test all(record -> record.purification_converged, diagnostics.history)
    @test all(record -> record.purification_termination_reason == :converged, diagnostics.history)
    @test all(record -> isapprox(record.trace, 2.0; atol=3e-3), diagnostics.history)
    @test all(record -> !isnothing(record.energy_total), diagnostics.history)
    @test isfinite(diagnostics.history[end].vh_residual)
    @test isfinite(diagnostics.history[end].vf_residual)
    @test isfinite(diagnostics.history[end].rho_residual)
    @test isfinite(diagnostics.history[end].commutator_residual)
    @test isinf(diagnostics.history[1].two_cycle_residual)

    limited = System(parameters_1d(
        L=2,
        t=-0.7,
        W=x -> (0.0, 0.13, -0.08, 0.21)[Int(x) + 1],
        U=0.2,
        S=nothing,
        purification_steps=35,
        scf_mixing=0.4,
        scf_tol=0.1,
        scf_max_iterations=1,
    ))
    @test !run_scf!(limited, -5.0, 5.0; purification_method=:sp2, verbose=:nothing)
    limited_diagnostics = scf_diagnostics(limited)
    @test !limited_diagnostics.converged
    @test limited_diagnostics.termination_reason == :max_iterations
    @test length(limited_diagnostics.history) == 1
end

@testset "P4.2 SCF stability and two-cycle criteria" begin
    function record(; rho_residual, two_cycle_residual=Inf)
        return SCFIterationRecord(
            1,
            2.0,
            5e-4,
            6e-4,
            rho_residual,
            7e-4,
            two_cycle_residual,
            true,
            :converged,
            10,
            nothing,
        )
    end

    tolerance = 1e-3
    stable = record(rho_residual=8e-4)
    unstable = record(rho_residual=2e-3)
    @test MPO_MeanField._scf_record_within_tolerance(stable, tolerance)
    @test !MPO_MeanField._scf_record_within_tolerance(unstable, tolerance)
    @test !MPO_MeanField._scf_has_stable_history([stable], tolerance; required_iterations=2)
    @test MPO_MeanField._scf_has_stable_history([stable, stable], tolerance; required_iterations=2)
    @test !MPO_MeanField._scf_has_stable_history([stable, unstable], tolerance; required_iterations=2)
    @test_throws ArgumentError MPO_MeanField._scf_has_stable_history([stable], tolerance; required_iterations=0)

    two_cycle = record(rho_residual=2e-3, two_cycle_residual=5e-4)
    converged_step = record(rho_residual=5e-4, two_cycle_residual=5e-4)
    nonrepeating_step = record(rho_residual=2e-3, two_cycle_residual=2e-3)
    @test MPO_MeanField._scf_is_two_cycle(two_cycle, tolerance)
    @test !MPO_MeanField._scf_is_two_cycle(converged_step, tolerance)
    @test !MPO_MeanField._scf_is_two_cycle(nonrepeating_step, tolerance)
end
