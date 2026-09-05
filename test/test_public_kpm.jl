@testset "Public KPM equal-square compatibility wrapper" begin
    model = SquareModel(
        size=(4, 4), hopping=(-1.0, -1.0), interaction=0.2,
        potential=(x, y) -> 0.03 * (Int(x) + Int(y)), filling=0.5,
    )
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-10, maxdim=32)
    solver = SCFSettings(mixing=0.5, tolerance=0.1, maxiter=30)
    settings = KPMSettings(
        moments=1024, probes=16, audit_moments=1024, audit_probes=16,
    )

    report = preflight(model, representation, solver; method=:kpm, kpm=settings)
    @test report.runnable
    @test report.spectral_bounds === nothing
    @test report.runtime_report.execution_path == :cpu_end_to_end

    result = solve(model, representation, solver; method=:kpm, kpm=settings)
    dense = solve(model, representation, solver; method=:dense)
    @test result.method == :kpm
    @test result.kpm === settings
    @test result.diagnostics isa KPMDiagnostics
    @test result.converged
    @test isapprox(result.observables.particle_number, 8.0; atol=2e-3)
    @test abs(result.observables.energy.total - dense.observables.energy.total) < 2e-3
    @test result.observables.idempotency_residual < 2e-3
    @test result.diagnostics.audit.trace_error < 2e-3
    @test result.diagnostics.audit.density_rms_difference < 2e-3

    @test_throws ArgumentError solve(
        ChainModel(size=8, hopping=-1.0, interaction=0.0),
        representation, solver; method=:kpm, kpm=KPMSettings(probes=8, audit_probes=8),
    )
    @test_throws ArgumentError KPMSettings(probes=0)

    mktempdir() do temporary_directory
        output = joinpath(temporary_directory, "kpm_result")
        write_result(result, output)
        stored = read_result(output)
        @test stored.metadata["method"] == "kpm"
        @test stored.metadata["solver"] == "public_sparse_vector_kpm_open_graph_hf"
        @test stored.input["kpm_moments"] == settings.moments
        @test length(stored.history) == length(result.diagnostics.history)
        @test length(stored.site_density) == 16
    end
end
