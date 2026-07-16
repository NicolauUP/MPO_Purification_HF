@testset "M4.2 square SCF is explicitly blocked" begin
    sys = System(parameters_square(U=0.3))
    @test_throws ArgumentError run_scf!(
        sys, -5.0, 5.0;
        purification_method=:sp2,
        verbose=:nothing,
    )
end
