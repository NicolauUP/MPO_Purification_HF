@testset "Public generic campaign runner" begin
    model = ChainModel(size=8, hopping=-0.7, interaction=0.0, filling=0.5)
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64)
    solver = SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1,
        maxiter=5, purification_maxiter=35,
    )
    campaign = CampaignSpec(
        name="public campaign smoke",
        cases=[CaseSpec(
            label="clean chain", model=model, representation=representation, solver=solver,
            spectral_bounds=(-3.0, 3.0),
        )],
    )

    @test campaign_case(campaign, 1).label == "clean chain"
    @test_throws ArgumentError campaign_case(campaign, 2)
    @test preflight(campaign; task=1, method=:mpo).runnable
    @test preflight(campaign; task=1, method=:dense).runnable

    mktempdir() do temporary_directory
        source = joinpath(temporary_directory, "campaign.jl")
        write(source, """
            using MPO_MeanField
            campaign = CampaignSpec(
                name=\"public campaign smoke\",
                cases=[CaseSpec(
                    label=\"clean chain\",
                    model=ChainModel(size=8, hopping=-0.7, interaction=0.0, filling=0.5),
                    representation=QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64),
                    solver=SCFSettings(purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=5, purification_maxiter=35),
                    spectral_bounds=(-3.0, 3.0),
                )],
            )
        """)
        output_root = joinpath(temporary_directory, "results")
        script = joinpath(@__DIR__, "..", "bin", "mpohf.jl")
        command = `$(Base.julia_cmd()) --startup-file=no --project=$(joinpath(@__DIR__, "..")) $script run $source --task 1 --method dense --backend cpu --output-root $output_root`
        output = read(command, String)
        @test occursin("SCF: converged=true", output)
        directory = joinpath(output_root, "public_campaign_smoke", "task_0001_clean_chain")
        @test isfile(joinpath(directory, "campaign.jl"))
        stored = read_result(directory)
        @test stored.metadata["method"] == "dense"
        @test stored.metadata["scf_converged"]
        @test isapprox(sum(stored.site_density), 4.0; atol=1e-12)
    end
end
