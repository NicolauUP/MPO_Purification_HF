using MPO_MeanField

# Full 4096 x 4096 SCF scaling gate for the strong-CDW branch. This keeps the
# validated 1024² physical model and solver conventions unchanged; only the
# QTT lattice depth grows from L=20 to L=24.
const TAU_STRONG_CDW_L12_SCF = sqrt(2.0) - 5.0 / 6.0

separable_aubry_andre_hopping_l12_scf(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L12_SCF * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L12_SCF * (Float64(y) + 0.5)),
)
checkerboard_seed_l12_scf(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

model = SquareModel(
    size=(4096, 4096),
    hopping=separable_aubry_andre_hopping_l12_scf(0.1),
    interaction=2.0,
    seed=checkerboard_seed_l12_scf(+2.0),
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
    name="separable_aubry_andre_strong_cdw_lside12_mpo_field_converged_chi512",
    cases=[CaseSpec(
        label="v2_0p1_u_2_seed_plus2",
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-23.0, 23.0),
    )],
)
