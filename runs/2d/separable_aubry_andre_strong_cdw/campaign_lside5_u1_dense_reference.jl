using MPO_MeanField

# Independent 32 x 32 reference for the lower-U cap study.  This matches the
# physical hopping, interaction, filling, boundary, and +2 checkerboard seed
# of the large U=1 calculation; dense HF is deliberately CPU-only.
const TAU_STRONG_CDW_L5_U1 = sqrt(2.0) - 5.0 / 6.0

separable_aubry_andre_hopping_l5_u1(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L5_U1 * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L5_U1 * (Float64(y) + 0.5)),
)
checkerboard_seed_l5_u1(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

model = SquareModel(
    size=(32, 32),
    hopping=separable_aubry_andre_hopping_l5_u1(0.1),
    interaction=1.0,
    seed=checkerboard_seed_l5_u1(+2.0),
    filling=0.5,
)

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside5_u1_dense_reference",
    cases=[CaseSpec(
        label="v2_0p1_u_1_seed_plus2_dense",
        model=model,
        # Dense HF does not use this representation, but an explicit setting
        # keeps the result input self-describing and physically matched.
        representation=QTTSettings(encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=1024),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=80,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            record_energy=true, stable_iterations=3, require_stationarity=true,
            measure_stationarity=true, detect_two_cycles=true,
            mixing_method=:linear,
        ),
    )],
)
