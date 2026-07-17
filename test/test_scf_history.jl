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
