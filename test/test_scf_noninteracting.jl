function verify_noninteracting_scf(params, method)
    sys = System(params)
    initial_H = dense_matrix(+(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    ), sys)
    bounds = exact_spectral_bounds(initial_H; padding=0.5)
    exact = exact_occupied_projector(dense_matrix(sys.H0, sys), 2)
    @test run_scf!(sys, bounds...; verbose=:nothing, purification_method=method)
    @test opnorm(dense_matrix(sys.ρ, sys) - exact) < 3e-3
    @test norm(dense_matrix(sys.VH, sys)) == 0.0
    @test norm(dense_matrix(sys.VF, sys)) == 0.0
end

@testset "M5.3 non-interacting SCF limit" begin
    configurations = (
        parameters_1d(U=0.0, S=nothing, purification_steps=30, scf_max_iterations=5),
        parameters_1d(U=0.0, S=x -> iseven(Int(x)) ? 0.2 : -0.2, purification_steps=30, scf_max_iterations=5),
        parameters_1d(U=0.0, W=x -> 0.25, purification_steps=30, scf_max_iterations=5),
        parameters_1d(U=0.0, t=x -> iseven(Int(x)) ? -0.7 : -0.3, purification_steps=30, scf_max_iterations=5),
    )
    for params in configurations, method in (:sp2, :palser_manolopoulos)
        verify_noninteracting_scf(params, method)
    end

    for mixing in (0.2, 0.8), method in (:sp2, :palser_manolopoulos)
        verify_noninteracting_scf(
            parameters_1d(U=0.0, S=nothing, scf_mixing=mixing, purification_steps=30, scf_max_iterations=5),
            method,
        )
    end
end
