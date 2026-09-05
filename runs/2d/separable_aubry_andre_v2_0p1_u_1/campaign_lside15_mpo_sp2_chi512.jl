using MPO_MeanField

# Billion-site extension of the converged 2048x2048 V2=0.1, U=1 campaign.
# Only the physical lattice size changes from the validated parent campaign.
const TAU = sqrt(2.0) - 5.0 / 6.0

separable_hopping(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU * (Float64(y) + 0.5)),
)

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

model = SquareModel(
    size=(32768, 32768),
    hopping=separable_hopping(0.1),
    interaction=1.0,
    seed=checkerboard_seed(+2.0),
    filling=0.5,
)

representation = QTTSettings(
    encoding=:interleaved,
    tci_tol=1e-10,
    cutoff=1e-10,
    maxdim=512,
)

solver = SCFSettings(
    purification=:sp2,
    mixing=0.5,
    tolerance=0.1,
    maxiter=20,
    purification_maxiter=60,
    square_fock_method=:binary_carry,
    sp2_idempotency_tolerance=2e-4,
    sp2_relative_trace_tolerance=1e-6,
    record_energy=false,
    stable_iterations=3,
    require_stationarity=false,
    measure_stationarity=false,
    detect_two_cycles=true,
    mixing_method=:pulay,
    pulay_history=4,
    pulay_warmup=3,
    pulay_regularization=1e-12,
    pulay_coefficient_limit=8.0,
    pulay_step_limit=20.0,
)

campaign = CampaignSpec(
    name="separable_aubry_andre_v2_0p1_u_1_lside15_mpo_sp2_chi512",
    cases=[CaseSpec(
        label="v2_0p1_u_1_seed_plus2_chi512",
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-23.0, 23.0),
    )],
)
