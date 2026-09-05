using MPO_MeanField

# Matched 1024x1024 control for the V2=0.1 quasiperiodic-hopping calculation.
# Only V2 is changed: U=1, filling, +2 checkerboard seed, QTT tolerance,
# spectral interval, SP2 thresholds, Pulay policy, and maxdim are identical.
checkerboard_seed_l10_u1_v2zero(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside10_u1_v2zero_chi768",
    cases=[CaseSpec(
        label="v2_0_u_1_seed_plus2_chi_768",
        model=SquareModel(
            size=(1024, 1024),
            hopping=((x, y) -> -1.0, (x, y) -> -1.0),
            interaction=1.0, seed=checkerboard_seed_l10_u1_v2zero(+2.0), filling=0.5,
        ),
        representation=QTTSettings(
            encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=768,
        ),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=20,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            record_energy=false, stable_iterations=3, require_stationarity=false,
            measure_stationarity=false, detect_two_cycles=true,
            mixing_method=:pulay, pulay_history=4, pulay_warmup=3,
            pulay_regularization=1e-12, pulay_coefficient_limit=8.0, pulay_step_limit=20.0,
        ),
        spectral_bounds=(-23.0, 23.0),
    )],
)
