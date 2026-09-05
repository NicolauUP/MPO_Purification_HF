function _rectangular_dense_reference(model::SquareModel)
    data = MPO_MeanField._kpm_data(model)
    H = Matrix(MPO_MeanField._kpm_hamiltonian(data, zeros(data.N), zeros(length(data.bonds))))
    occupation = zeros(Float64, data.N)
    occupation[1:round(Int, model.filling * data.N)] .= 1.0
    vectors = eigen(Symmetric(H)).vectors
    rho = vectors * Diagonal(occupation) * vectors'
    density = diag(rho)
    order = [real((rho[i, j] + rho[j, i]) / 2) for (i, j, _) in data.bonds]
    return data, density, order, MPO_MeanField._kpm_energy(
        data, density, order, Float64(model.interaction),
    )
end

function _rectangular_dense_hf_reference(model::SquareModel, solver::SCFSettings)
    data = MPO_MeanField._kpm_data(model)
    Ne = round(Int, model.filling * data.N)
    VH, VF = copy(data.seed), zeros(Float64, length(data.bonds))
    density = zeros(Float64, data.N)
    order = zeros(Float64, length(data.bonds))
    for iteration in 1:solver.maxiter
        H = Matrix(MPO_MeanField._kpm_hamiltonian(data, VH, VF))
        vectors = eigen(Symmetric(H)).vectors[:, 1:Ne]
        rho = vectors * vectors'
        density = diag(rho)
        order = [real((rho[i, j] + rho[j, i]) / 2) for (i, j, _) in data.bonds]
        next_VH = zeros(Float64, data.N)
        next_VF = zeros(Float64, length(data.bonds))
        for (k, (i, j, _)) in enumerate(data.bonds)
            next_VH[i] += model.interaction * density[j]
            next_VH[j] += model.interaction * density[i]
            next_VF[k] = -model.interaction * order[k]
        end
        if iteration > 1 && max(
            norm(next_VH - VH) / max(norm(VH), sqrt(eps())),
            norm(next_VF - VF) / max(norm(VF), sqrt(eps())),
        ) <= solver.tolerance / 100
            break
        end
        if iteration == 1
            VH, VF = next_VH, next_VF
        else
            VH = solver.mixing .* next_VH .+ (1 - solver.mixing) .* VH
            VF = solver.mixing .* next_VF .+ (1 - solver.mixing) .* VF
        end
    end
    return data, density, order, MPO_MeanField._kpm_energy(
        data, density, order, Float64(model.interaction),
    )
end

@testset "Public KPM rectangular sparse-graph adapter" begin
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-10, maxdim=32)
    solver = SCFSettings(mixing=0.5, tolerance=0.1, maxiter=30)

    for size in ((4, 2), (2, 4), (8, 4))
        model = SquareModel(
            size=size, hopping=(-1.0, -0.7), interaction=0.0,
            potential=(x, y) -> 0.07 * x - 0.03 * y,
        )
        data = MPO_MeanField._kpm_data(model)
        nx, ny = size
        @test data.N == nx * ny
        @test length(data.bonds) == (nx - 1) * ny + nx * (ny - 1)
        @test sort(data.codes) == collect(0:(data.N - 1))
        @test data.codes[1] == 0

        reference_hamiltonian = zeros(Float64, data.N, data.N)
        for i in 1:data.N
            reference_hamiltonian[i, i] = data.onsite[i]
        end
        for (bond, hopping) in zip(data.bonds, data.hopping)
            i, j, _ = bond
            reference_hamiltonian[i, j] = hopping
            reference_hamiltonian[j, i] = hopping
        end
        @test Matrix(MPO_MeanField._kpm_hamiltonian(
            data, zeros(data.N), zeros(length(data.bonds)),
        )) == reference_hamiltonian
    end

    model = SquareModel(
        size=(4, 2), hopping=(-1.0, -0.7), interaction=0.0,
        potential=(x, y) -> 0.11 * x - 0.07 * y,
    )
    settings = KPMSettings(moments=2048, probes=8, audit_moments=2048, audit_probes=8)
    report = preflight(model, representation, solver; method=:kpm, kpm=settings)
    @test report.runnable
    @test report.geometry == :rectangle
    @test !any(issue -> issue.code == :rectangular_backend_unavailable, report.issues)

    result = solve(model, representation, solver; method=:kpm, kpm=settings)
    data, density, order, energy = _rectangular_dense_reference(model)
    @test result.converged
    @test isapprox(result.observables.particle_number, 4.0; atol=2e-3)
    @test norm(result.observables.site_density - density) / sqrt(data.N) < 2e-3
    @test norm(real.(result.observables.horizontal_bond_order) |> collect) >= 0
    @test abs(result.observables.energy.total - energy.total) < 2e-3

    interacting = SquareModel(
        size=(4, 2), hopping=(-1.0, -0.7), interaction=0.2,
        potential=(x, y) -> 0.11 * x - 0.07 * y,
        seed=(x, y) -> iseven(x + y) ? 0.03 : -0.03,
    )
    interacting_result = solve(interacting, representation, solver; method=:kpm, kpm=settings)
    data, density, order, energy = _rectangular_dense_hf_reference(interacting, solver)
    @test interacting_result.converged
    @test norm(interacting_result.observables.site_density - density) / sqrt(data.N) < 3e-3
    @test abs(interacting_result.observables.energy.total - energy.total) < 3e-3

    mktempdir() do temporary_directory
        output = joinpath(temporary_directory, "rectangular_kpm")
        write_result(result, output)
        stored = read_result(output)
        @test stored.metadata["geometry"] == "rectangle"
        @test stored.metadata["kpm_probe_ordering"] == "interleaved_coordinate_bits"
        @test stored.input["kpm_site_storage"] == "row_major_x_fastest"
        @test length(stored.site_density) == 8
    end
end
