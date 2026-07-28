@testset "square TCI/binary-carry full SCF equivalence" begin
    checkerboard(amplitude) =
        (x, y) -> iseven(Int(x) + Int(y)) ? amplitude : -amplitude
    params = parameters_square(
        L=4,
        t=(-0.6, -0.35),
        U=0.3,
        W=checkerboard(0.6),
        S=checkerboard(0.05),
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=64,
        density=0.5,
        purification_steps=50,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=30,
    )
    tci = System(params)
    carry = System(params)
    bounds = square_scf_spectral_bounds(
        params; potential_bounds=(-0.6, 0.6), margin=0.5,
    )
    common_options = (
        purification_method=:sp2,
        sp2_idempotency_tolerance=2e-4,
        sp2_trace_tolerance=8e-6,
        verbose=:nothing,
    )
    @test run_scf!(
        tci, bounds...; common_options..., square_fock_method=:tci,
    )
    @test run_scf!(
        carry, bounds...; common_options..., square_fock_method=:binary_carry,
    )

    tci_rho = dense_matrix(tci.ρ, tci)
    carry_rho = dense_matrix(carry.ρ, carry)
    tci_vh = dense_matrix(tci.VH, tci)
    carry_vh = dense_matrix(carry.VH, carry)
    tci_vf = dense_matrix(tci.VF, tci)
    carry_vf = dense_matrix(carry.VF, carry)
    tci_observables = observables_square(tci)
    carry_observables = observables_square(carry)
    staggered_order(observables) = sum(begin
        x, y = square_lattice_decoder(site - 1, params.L)
        (-1)^(x + y) * (
            observables.site_density[site] -
            sum(observables.site_density) / length(observables.site_density)
        )
    end for site in eachindex(observables.site_density)) /
        length(observables.site_density)

    # Both paths solve the same physical fixed point. These tolerances are
    # much tighter than the configured 1e-3 SCF residual and 2e-4 SP2
    # idempotency thresholds; they detect a convention or boundary mismatch.
    @test length(scf_diagnostics(tci).history) ==
          length(scf_diagnostics(carry).history)
    @test opnorm(tci_rho - carry_rho) < 1e-8
    @test opnorm(tci_vh - carry_vh) < 1e-8
    @test opnorm(tci_vf - carry_vf) < 1e-8
    @test maximum(abs.(
        tci_observables.site_density .- carry_observables.site_density,
    )) < 1e-8
    @test abs(
        tci_observables.energy.total - carry_observables.energy.total,
    ) < 1e-8
    @test abs(
        staggered_order(tci_observables) -
        staggered_order(carry_observables),
    ) < 1e-8
    @test tci_observables.idempotency_residual < 2e-4
    @test carry_observables.idempotency_residual < 2e-4
    @test tci_observables.stationarity_residual < 1e-6
    @test carry_observables.stationarity_residual < 1e-6

    # The unqualified square selector now chooses the validated fast path.
    default_hartree, default_fock = extract_mean_fields(carry)
    explicit_hartree, explicit_fock = extract_mean_fields(
        carry; square_fock_method=:binary_carry,
    )
    @test opnorm(
        dense_matrix(default_hartree, carry) -
        dense_matrix(explicit_hartree, carry),
    ) < 1e-12
    @test opnorm(
        dense_matrix(default_fock, carry) -
        dense_matrix(explicit_fock, carry),
    ) < 1e-12
end
