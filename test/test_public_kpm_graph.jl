@testset "Public KPM explicit sparse graph adapter" begin
    hopping = sparse(
        [1, 2, 2, 3, 3, 4, 4, 1],
        [2, 1, 3, 2, 4, 3, 1, 4],
        [-1.0, -1.0, -0.7, -0.7, -0.9, -0.9, -0.5, -0.5],
        4, 4,
    )
    model = GraphModel(
        hopping=hopping, interaction=0.2,
        potential=[0.1, -0.05, 0.08, -0.03],
        seed=[0.01, -0.01, 0.01, -0.01],
        probe_codes=[0, 2, 1, 3],
    )
    representation = QTTSettings()
    solver = SCFSettings(mixing=0.5, tolerance=0.1, maxiter=30)
    settings = KPMSettings(moments=1024, probes=4, audit_moments=1024, audit_probes=4)

    data = MPO_MeanField._kpm_data(model)
    @test data.N == 4
    @test length(data.bonds) == 4
    @test data.codes == [0, 2, 1, 3]
    @test Matrix(MPO_MeanField._kpm_hamiltonian(
        data, zeros(4), zeros(4),
    )) == Matrix(hopping) + Diagonal(model.potential)

    report = preflight(model, representation, solver; method=:kpm, kpm=settings)
    @test report.runnable
    @test report.model_kind == :graph
    @test report.geometry == :graph
    @test !preflight(model, representation, solver; method=:dense).runnable
    @test_throws ArgumentError legacy_parameters(model, representation, solver)

    result = solve(model, representation, solver; method=:kpm, kpm=settings)
    @test result.converged
    @test isapprox(result.observables.particle_number, 2.0; atol=2e-3)
    @test length(result.observables.bonds) == 4
    @test all(bond -> bond[3] == :graph, result.observables.bonds)

    mktempdir() do temporary_directory
        output = joinpath(temporary_directory, "graph_kpm")
        write_result(result, output)
        stored = read_result(output)
        @test stored.metadata["geometry"] == "graph"
        @test stored.metadata["kpm_probe_ordering"] == "explicit_probe_codes"
        @test all(bond -> bond.orientation == :graph, stored.bond_order)
    end

    @test_throws ArgumentError GraphModel(hopping=sparse([1, 2], [2, 1], [1.0, 0.5], 4, 4), interaction=0.0)
    @test_throws ArgumentError GraphModel(hopping=hopping, interaction=0.0, probe_codes=[0, 1, 2, 2])
end
