using MPO_MeanField

# Two nearby irrational C8 reciprocal-space stars. Their difference wave
# number creates a supermoiré beat length of approximately 226 lattice sites,
# while each individual star retains the short C8 hopping modulation.
const DOUBLE_C8_TAU_1 = sqrt(2.0) / 8.0
const DOUBLE_C8_TAU_2 = (41.0 / 40.0) * DOUBLE_C8_TAU_1
const DOUBLE_C8_ANGLES = (0.0, pi / 4, pi / 2, 3pi / 4)
const DOUBLE_C8_V_1 = 0.05
const DOUBLE_C8_V_2 = 0.05

function c8_star(x::Real, y::Real, tau::Real)
    return sum(cos(2pi * Float64(tau) *
        (Float64(x) * cos(angle) + Float64(y) * sin(angle)))
        for angle in DOUBLE_C8_ANGLES) / length(DOUBLE_C8_ANGLES)
end

function double_c8_modulation(x::Real, y::Real)
    return DOUBLE_C8_V_1 * c8_star(x, y, DOUBLE_C8_TAU_1) +
           DOUBLE_C8_V_2 * c8_star(x, y, DOUBLE_C8_TAU_2)
end

double_c8_hopping = (
    (x, y) -> -1.0 - double_c8_modulation(Float64(x) + 0.5, Float64(y)),
    (x, y) -> -1.0 - double_c8_modulation(Float64(x), Float64(y) + 0.5),
)

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function dense_solver()
    return SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=1e-8, maxiter=100,
        purification_maxiter=60, square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
        record_energy=true, stable_iterations=3, require_stationarity=true,
        measure_stationarity=true, detect_two_cycles=true, mixing_method=:linear,
    )
end

function mpo_solver()
    return SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=60,
        purification_maxiter=60, square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
        record_energy=false, stable_iterations=3, require_stationarity=false,
        measure_stationarity=false, detect_two_cycles=true, mixing_method=:linear,
    )
end

case_for(side::Int, maxdim::Int, label::String, solver::SCFSettings) = CaseSpec(
    label=label,
    model=SquareModel(
        size=(side, side), hopping=double_c8_hopping, interaction=2.0,
        seed=checkerboard_seed(+2.0), filling=0.5,
    ),
    representation=QTTSettings(
        encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim,
    ),
    solver=solver,
    spectral_bounds=(-8.5, 12.5),
)

campaign = CampaignSpec(
    name="double_c8_supermoire_validation",
    cases=[
        # Task 1: independent dense reference.
        case_for(64, 512, "double_c8_l64_u2_v005x2_dense", dense_solver()),
        # Task 2: same model in MPO form, with a representation cap ladder.
        case_for(128, 256, "double_c8_l128_u2_v005x2_chi256", mpo_solver()),
        case_for(128, 512, "double_c8_l128_u2_v005x2_chi512", mpo_solver()),
    ],
)
