@testset "M7.1 host GC policy" begin
    @test !MPO_MeanField._should_collect_garbage(1, :automatic, 2, 1)
    @test MPO_MeanField._should_collect_garbage(1, :forced, 2, typemax(Int))
    @test !MPO_MeanField._should_collect_garbage(1, :periodic, 2, 1)
    @test MPO_MeanField._should_collect_garbage(2, :periodic, 2, 1)
    @test MPO_MeanField._should_collect_garbage(1, :threshold, 2, 1)
    @test_throws ArgumentError maybe_collect_garbage!(1; gc_policy=:unknown)
    @test_throws ArgumentError maybe_collect_garbage!(1; gc_policy=:periodic, gc_period=0)

    params = parameters_1d(
        t=-0.7,
        W=x -> (-0.2, 0.1, -0.05, 0.25)[Int(x) + 1],
        U=0.3,
        purification_steps=35,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=20,
    )
    reference = nothing
    for (policy, period, threshold) in (
        (:automatic, 10, 1 << 30),
        (:forced, 10, 1 << 30),
        (:periodic, 2, 1 << 30),
        (:threshold, 10, 1),
    )
        sys = System(params)
        @test run_scf!(sys, -5.0, 5.0;
            purification_method=:sp2,
            gc_policy=policy,
            gc_period=period,
            gc_threshold_bytes=threshold,
        )
        rho = dense_matrix(sys.ρ, sys)
        if isnothing(reference)
            reference = rho
        else
            @test opnorm(rho - reference) < 1e-10
        end
    end

    # PM has its own inner iteration boundary; accepting the non-default
    # policy here verifies propagation from run_scf! through the dispatcher.
    sys_pm = System(params)
    @test run_scf!(sys_pm, -5.0, 5.0;
        purification_method=:palser_manolopoulos,
        gc_policy=:periodic,
        gc_period=2,
    )
end
