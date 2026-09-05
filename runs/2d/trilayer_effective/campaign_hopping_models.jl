using MPO_MeanField

# Four effective hopping models for a single active square-lattice band. The
# substrate is common to all cases and has independent x/y wavelengths.
const A_HOP = 0.10
const TAU_FAST = sqrt(2.0) - 5.0 / 6.0
const TAU_SECOND = sqrt(5.0) - 2.0
const TAU_PRODUCT_Y = sqrt(3.0) - 1.0
const TAU_SUBSTRATE = 1.0 / 32.0
const THETA_ROTATED = pi / 9
const THETA_SECOND = 2pi / 7

wave(x, y, tau, theta) = cos(2pi * tau *
    (Float64(x) * cos(theta) + Float64(y) * sin(theta)))

# Substrate: unrotated and isotropic in the square coordinates. The substrate
# wavelength differs from the active-layer hopping wavelength, but tau_x=tau_y
# as requested for the aligned graphene/hBN analogue.
substrate_potential(x, y) = 0.20 * cos(2pi * TAU_SUBSTRATE * Float64(x)) +
                            0.20 * cos(2pi * TAU_SUBSTRATE * Float64(y))

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

# 1. One-dimensional hopping control (x-bond modulation only).
function hopping_1d_x()
    tx = (x, y) -> -1.0 - A_HOP * cos(2pi * TAU_FAST * (Float64(x) + 0.5))
    ty = (x, y) -> -1.0
    return (tx, ty)
end

# 2. A single generic-angle rotated wave, evaluated at bond centers.
function hopping_rotated()
    tx = (x, y) -> -1.0 - A_HOP * wave(Float64(x) + 0.5, y, TAU_FAST, THETA_ROTATED)
    ty = (x, y) -> -1.0 - A_HOP * wave(x, Float64(y) + 0.5, TAU_FAST, THETA_ROTATED)
    return (tx, ty)
end

# 3. A genuinely mixed, non-separable product modulation.
function hopping_product()
    product_wave(x, y) = cos(2pi * TAU_FAST * Float64(x)) *
                         cos(2pi * TAU_PRODUCT_Y * Float64(y))
    tx = (x, y) -> -1.0 - A_HOP * product_wave(Float64(x) + 0.5, y)
    ty = (x, y) -> -1.0 - A_HOP * product_wave(x, Float64(y) + 0.5)
    return (tx, ty)
end

# 4. Two independent layer-like waves: a supermoiré hopping control.
function hopping_two_layer()
    two_layer(x, y) = 0.5 * wave(x, y, TAU_FAST, 0.0) +
                      0.5 * wave(x, y, TAU_SECOND, THETA_SECOND)
    tx = (x, y) -> -1.0 - A_HOP * two_layer(Float64(x) + 0.5, y)
    ty = (x, y) -> -1.0 - A_HOP * two_layer(x, Float64(y) + 0.5)
    return (tx, ty)
end

function mpo_solver()
    return SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=30,
        purification_maxiter=60, square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
        record_energy=false, stable_iterations=3, require_stationarity=false,
        measure_stationarity=false, detect_two_cycles=true, mixing_method=:pulay,
        pulay_history=4, pulay_warmup=3, pulay_regularization=1e-12,
        pulay_coefficient_limit=8.0, pulay_step_limit=20.0,
    )
end

function model_case(label, hopping)
    return CaseSpec(
        label=label,
        model=SquareModel(
            size=(128, 128), hopping=hopping, interaction=1.0,
            potential=substrate_potential, seed=checkerboard_seed(+2.0), filling=0.5,
        ),
        representation=QTTSettings(
            encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=512,
        ),
        # For U=1 and |t|<=1.1, the conservative HF interval with 0.5 padding
        # is [-13.3,13.3]; the earlier [-12,12] interval was too narrow.
        solver=mpo_solver(), spectral_bounds=(-13.3, 13.3),
    )
end

campaign = CampaignSpec(
    name="trilayer_effective_hopping_model_comparison",
    cases=[
        model_case("hopping_1d_x", hopping_1d_x()),
        model_case("hopping_rotated_generic_angle", hopping_rotated()),
        model_case("hopping_product_nonseparable", hopping_product()),
        model_case("hopping_two_layer_supermoire", hopping_two_layer()),
    ],
)
