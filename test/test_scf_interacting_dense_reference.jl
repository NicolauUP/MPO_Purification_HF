function dense_hf_1d(H0, U, Ne; mixing=0.5, tolerance=1e-3, max_iterations=20, initial_vh=nothing, initial_vf=nothing)
    N = size(H0, 1)
    VH = isnothing(initial_vh) ? zeros(Float64, N, N) : copy(initial_vh)
    VF = isnothing(initial_vf) ? zeros(Float64, N, N) : copy(initial_vf)
    rho_previous = zeros(Float64, N, N)
    for iteration in 1:max_iterations
        rho = exact_occupied_projector(H0 + VH + VF, Ne)
        new_VH = zeros(Float64, N, N)
        new_VF = zeros(Float64, N, N)
        for i in 1:N
            for j in (i - 1, i + 1)
                1 <= j <= N && (new_VH[i, i] += U * rho[j, j])
            end
            if i < N
                new_VF[i, i + 1] = -U * real(rho[i, i + 1])
                new_VF[i + 1, i] = new_VF[i, i + 1]
            end
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
    error("dense HF reference did not converge")
end

function verify_interacting_scf_against_dense(params)
    sys = System(params)
    H0 = dense_matrix(sys.H0, sys)
    rho_dense, vh_dense, vf_dense = dense_hf_1d(
        H0, params.U, 2;
        mixing=params.scf_mixing,
        max_iterations=params.scf_max_iterations,
        initial_vh=dense_matrix(sys.VH, sys),
        initial_vf=dense_matrix(sys.VF, sys),
    )
    @test run_scf!(sys, -5.0, 5.0; verbose=:nothing, purification_method=:sp2)
    @test opnorm(dense_matrix(sys.ρ, sys) - rho_dense) < 4e-3
    @test opnorm(dense_matrix(sys.VH, sys) - vh_dense) < 4e-3
    @test opnorm(dense_matrix(sys.VF, sys) - vf_dense) < 4e-3
end

@testset "M5.4 weak-coupling 1D SCF dense reference" begin
    params = parameters_1d(
        t=-0.7,
        W=x -> (-0.2, 0.1, -0.05, 0.25)[Int(x) + 1],
        U=0.3,
        S=nothing,
        purification_steps=35,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=20,
    )
    verify_interacting_scf_against_dense(params)

    staggered_seed = parameters_1d(
        t=-0.7,
        W=x -> (-0.2, 0.1, -0.05, 0.25)[Int(x) + 1],
        U=0.6,
        S=x -> iseven(Int(x)) ? 0.15 : -0.15,
        purification_steps=35,
        scf_mixing=0.3,
        scf_tol=0.1,
        scf_max_iterations=40,
    )
    verify_interacting_scf_against_dense(staggered_seed)
end
