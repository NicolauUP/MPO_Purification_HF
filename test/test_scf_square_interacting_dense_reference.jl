function dense_hf_square(
    H0,
    U,
    Ne,
    L;
    mixing=0.5,
    tolerance=1e-3,
    max_iterations=40,
    initial_vh=nothing,
    initial_vf=nothing,
)
    N = size(H0, 1)
    VH = isnothing(initial_vh) ? zeros(Float64, N, N) : copy(initial_vh)
    VF = isnothing(initial_vf) ? zeros(Float64, N, N) : copy(initial_vf)
    rho_previous = zeros(Float64, N, N)
    for iteration in 1:max_iterations
        rho = exact_occupied_projector(H0 + VH + VF, Ne)
        new_VH = zeros(Float64, N, N)
        new_VF = zeros(Float64, N, N)
        for site in 1:N
            for neighbour in values(square_neighbours(site, L))
                !isnothing(neighbour) && (new_VH[site, site] += U * rho[neighbour, neighbour])
            end
        end
        for (site, neighbour, _) in square_undirected_bonds(L)
            new_VF[site, neighbour] = -U * real(rho[site, neighbour])
            new_VF[neighbour, site] = new_VF[site, neighbour]
        end
        if iteration > 1
            residuals = (
                norm(new_VH - VH) / max(norm(VH), sqrt(eps())),
                norm(new_VF - VF) / max(norm(VF), sqrt(eps())),
                norm(rho - rho_previous) / max(norm(rho_previous), sqrt(eps())),
                norm((H0 + VH + VF) * rho - rho * (H0 + VH + VF)) /
                    max(norm((H0 + VH + VF) * rho), sqrt(eps())),
            )
            maximum(residuals) < tolerance && return rho, VH, VF
        end
        rho_previous = rho
        VH = mixing * new_VH + (1 - mixing) * VH
        VF = mixing * new_VF + (1 - mixing) * VF
    end
    error("dense square HF reference did not converge")
end

@testset "P2.2 interacting square SCF dense reference" begin
    params = parameters_square(
        L=4,
        t=(-0.6, -0.35),
        W=(x, y) -> 0.11x + 0.07y + 0.013x * y,
        U=0.15,
        S=nothing,
        density=0.5,
        purification_steps=35,
        itensors_maxdim=64,
        scf_mixing=0.4,
        scf_tol=0.1,
        scf_max_iterations=40,
    )
    sys = System(params)
    H0 = dense_matrix(sys.H0, sys)
    rho_dense, vh_dense, vf_dense = dense_hf_square(
        H0,
        params.U,
        8,
        params.L;
        mixing=params.scf_mixing,
        max_iterations=params.scf_max_iterations,
        initial_vh=dense_matrix(sys.VH, sys),
        initial_vf=dense_matrix(sys.VF, sys),
    )
    @test run_scf!(sys, -5.0, 5.0; purification_method=:sp2, verbose=:nothing)
    rho = dense_matrix(sys.ρ, sys)
    vh = dense_matrix(sys.VH, sys)
    vf = dense_matrix(sys.VF, sys)
    @test opnorm(rho - rho_dense) < 4e-3
    @test opnorm(vh - vh_dense) < 4e-3
    @test opnorm(vf - vf_dense) < 4e-3
    @test isapprox(real(tr(rho)), 8.0; atol=3e-3)
    @test opnorm((H0 + vh + vf) * rho - rho * (H0 + vh + vf)) < 1e-5
end
