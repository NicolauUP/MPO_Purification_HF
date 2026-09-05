using MPO_MeanField

const SIDE = parse(Int, get(ENV, "MPO_SCALE_SIDE", "128"))
const SEED_SIGN = parse(Float64, get(ENV, "MPO_SCALE_SEED", "1.0"))
const MAXDIM = parse(Int, get(ENV, "MPO_SCALE_MAXDIM", "256"))
const SEED_LABEL = SEED_SIGN > 0 ? "plus" : "minus"
const V2 = 0.0
const U = 1.0

separable_hopping(V::Real) = (
    (x, y) -> -1.0 - Float64(V) * cos(2π * (sqrt(2.0) - 5.0 / 6.0) * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V) * cos(2π * (sqrt(2.0) - 5.0 / 6.0) * (Float64(y) + 0.5)),
)

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

case = CaseSpec(
    label="scaling_control_v2_0_u_1_seed$(SEED_LABEL)_lside$(SIDE)_chi$(MAXDIM)",
    model=SquareModel(
        size=(SIDE, SIDE), hopping=separable_hopping(V2), interaction=U,
        seed=checkerboard_seed(SEED_SIGN), filling=0.5,
    ),
    representation=QTTSettings(
        encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=MAXDIM,
    ),
    solver=SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=30,
        purification_maxiter=50, square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6,
        record_energy=false, stable_iterations=2, detect_two_cycles=true,
        mixing_method=:linear,
    ),
    spectral_bounds=(-8.5, 12.5),
)

campaign = CampaignSpec(
    name="separable_aubry_andre_scaling_control_mpo_sp2",
    cases=[case],
)
