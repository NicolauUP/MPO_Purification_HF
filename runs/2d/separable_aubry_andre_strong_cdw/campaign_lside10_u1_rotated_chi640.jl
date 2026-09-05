using MPO_MeanField

# Memory-safe rank retry of the theta=pi/4 non-separable hopping model.
# This keeps all physical and numerical controls of the chi=768 run fixed;
# only the representation cap is reduced to test whether its late-SP2
# temporary contractions fit on one H100.
const TAU_ROTATED_CDW_L10_CHI640 = sqrt(2.0) - 5.0 / 6.0
const ROTATION_ROTATED_CDW_L10_CHI640 = pi / 4

function rotated_aubry_andre_hopping_l10_chi640(V2::Real)
    amplitude = Float64(V2)
    c, s = cos(ROTATION_ROTATED_CDW_L10_CHI640), sin(ROTATION_ROTATED_CDW_L10_CHI640)
    return (
        (x, y) -> -1.0 - amplitude * cos(2pi * TAU_ROTATED_CDW_L10_CHI640 *
            ((Float64(x) + 0.5) * c + Float64(y) * s)),
        (x, y) -> -1.0 - amplitude * cos(2pi * TAU_ROTATED_CDW_L10_CHI640 *
            (-Float64(x) * s + (Float64(y) + 0.5) * c)),
    )
end

checkerboard_seed_rotated_l10_chi640(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

campaign = CampaignSpec(
    name="rotated_aubry_andre_strong_cdw_lside10_u1_chi640",
    cases=[
        CaseSpec(
            label="v2_0p1_u_1_theta_pi4_seed_plus2_chi_640",
            model=SquareModel(
                size=(1024, 1024), hopping=rotated_aubry_andre_hopping_l10_chi640(0.1),
                interaction=1.0, seed=checkerboard_seed_rotated_l10_chi640(+2.0), filling=0.5,
            ),
            representation=QTTSettings(
                encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=640,
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
        ),
    ],
)
