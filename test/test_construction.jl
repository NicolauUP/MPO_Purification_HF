function parameters_1d(; kwargs...)
    defaults = (
        L=2, t=-1.0, U=0.3, W=nothing, S=nothing,
        tci_tol=1e-10, itensors_tol=1e-12, itensors_maxdim=32,
        density=0.5, purification_steps=20, scf_mixing=0.5,
        scf_tol=0.1, scf_max_iterations=5,
    )
    return Parameters1D(; merge(defaults, (; kwargs...))...)
end

function parameters_square(; kwargs...)
    defaults = (
        L=4, t=(-1.0, -1.0), U=0.3, W=nothing, S=nothing,
        tci_tol=1e-10, itensors_tol=1e-12, itensors_maxdim=32,
        density=0.5, purification_steps=20, scf_mixing=0.5,
        scf_tol=0.1, scf_max_iterations=5,
    )
    return ParametersSquare(; merge(defaults, (; kwargs...))...)
end

@testset "M0 parameter validation" begin
    @test System(parameters_1d()) isa System
    @test System(parameters_square()) isa System
    mixed_square = System(parameters_square(t=((x, y) -> -1.0, -0.5)))
    functional_square = System(parameters_square(t=((x, y) -> -1.0, (x, y) -> -0.5)))
    @test mixed_square isa System
    @test isapprox(
        dense_matrix(mixed_square.H0, mixed_square),
        dense_matrix(functional_square.H0, functional_square);
        atol=1e-10,
        rtol=1e-10,
    )

    zero_hopping = System(parameters_1d(t=x -> 0.0))
    @test norm(dense_matrix(zero_hopping.H0, zero_hopping)) == 0.0

    zero_potential = System(parameters_1d(W=x -> 0.0))
    no_potential = System(parameters_1d(W=nothing))
    @test isapprox(
        dense_matrix(zero_potential.H0, zero_potential),
        dense_matrix(no_potential.H0, no_potential);
        atol=1e-12,
        rtol=1e-12,
    )

    for invalid in (
        parameters_1d(L=0),
        parameters_1d(L=8*sizeof(Int)),
        parameters_1d(U=x -> x),
        parameters_1d(U=Inf),
        parameters_1d(t=NaN),
        parameters_1d(tci_tol=0.0),
        parameters_1d(tci_tol=Inf),
        parameters_1d(itensors_tol=0.0),
        parameters_1d(itensors_maxdim=0),
        parameters_1d(density=0.0),
        parameters_1d(density=1.0),
        parameters_1d(density=NaN),
        parameters_1d(purification_steps=0),
        parameters_1d(scf_mixing=-0.1),
        parameters_1d(scf_mixing=1.1),
        parameters_1d(scf_mixing=NaN),
        parameters_1d(scf_tol=0.0),
        parameters_1d(scf_tol=Inf),
        parameters_1d(scf_max_iterations=0),
        parameters_square(L=3),
        parameters_square(t=(-1.0,)),
        parameters_square(t=(-1.0, nothing)),
        parameters_square(t=(-1.0, Inf)),
    )
        @test_throws ArgumentError System(invalid)
    end
end
