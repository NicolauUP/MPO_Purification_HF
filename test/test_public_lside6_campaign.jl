@testset "Public 64x64 quasiperiodic MPO campaign" begin
    campaign_source = joinpath(
        @__DIR__, "..", "runs", "2d", "separable_aubry_andre_lside6", "campaign_public.jl",
    )
    campaign_module = Module(:PublicLside6Campaign)
    Core.eval(campaign_module, :(using MPO_MeanField))
    Base.include(campaign_module, campaign_source)
    campaign = getfield(campaign_module, :campaign)

    @test campaign isa CampaignSpec
    @test campaign.name == "separable_aubry_andre_lside6_public_mpo"
    @test length(campaign.cases) == 1
    case = campaign_case(campaign, 1)
    @test case.label == "v2_0p5_u_1_chi_512"
    @test case.model.size == (64, 64)
    @test case.model.interaction == 1.0
    @test case.model.filling == 0.5
    @test case.representation.encoding == :interleaved
    @test case.representation.tci_tol == 1e-10
    @test case.representation.cutoff == 1e-10
    @test case.representation.maxdim == 512
    @test case.spectral_bounds == (-8.5, 12.5)
    @test case.solver.purification == :sp2
    @test case.solver.square_fock_method == :binary_carry
    @test case.solver.sp2_idempotency_tolerance == 2e-4
    @test case.solver.sp2_relative_trace_tolerance == 1e-6
    @test !case.solver.record_energy
    @test case.solver.stable_iterations == 2
    @test case.solver.detect_two_cycles

    report = preflight(campaign; task=1, method=:mpo, runtime=RuntimeSettings(backend=:cpu))
    @test report.runnable
    @test report.physical_size == (64, 64)
    @test report.physical_sites == 4096
    @test report.target_particles == 2048
    @test report.spectral_bounds == (-8.5, 12.5)

    pulay_source = joinpath(
        @__DIR__, "..", "runs", "2d", "separable_aubry_andre_lside6",
        "campaign_public_pulay_chi256.jl",
    )
    pulay_module = Module(:PublicLside6PulayCampaign)
    Core.eval(pulay_module, :(using MPO_MeanField))
    Base.include(pulay_module, pulay_source)
    pulay_campaign = getfield(pulay_module, :campaign)
    pulay_case = campaign_case(pulay_campaign, 1)
    @test pulay_campaign.name == "separable_aubry_andre_lside6_public_mpo_pulay_chi256"
    @test pulay_case.representation.maxdim == 256
    @test pulay_case.solver.mixing_method == :pulay
    @test pulay_case.solver.pulay_history == 4
    @test pulay_case.solver.pulay_warmup == 3
end
