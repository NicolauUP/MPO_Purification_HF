using MPO_MeanField

const A_HOP = 0.10
const TAU_FAST = sqrt(2.0) - 5.0 / 6.0
const TAU_SECOND = sqrt(5.0) - 2.0
const TAU_PRODUCT_Y = sqrt(3.0) - 1.0
const TAU_SUBSTRATE = 1.0 / 32.0
const THETA_ROTATED = pi / 9
const THETA_SECOND = 2pi / 7

wave(x, y, tau, theta) = cos(2pi * tau * (Float64(x) * cos(theta) + Float64(y) * sin(theta)))
substrate_potential(x, y) = 0.20 * cos(2pi * TAU_SUBSTRATE * Float64(x)) + 0.20 * cos(2pi * TAU_SUBSTRATE * Float64(y))
checkerboard_seed(amplitude::Real) = (x, y) -> iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)
hopping_1d_x() = ((x, y) -> -1.0 - A_HOP * cos(2pi * TAU_FAST * (Float64(x) + 0.5)), (x, y) -> -1.0)
hopping_rotated() = ((x, y) -> -1.0 - A_HOP * wave(Float64(x) + 0.5, y, TAU_FAST, THETA_ROTATED), (x, y) -> -1.0 - A_HOP * wave(x, Float64(y) + 0.5, TAU_FAST, THETA_ROTATED))
function hopping_product()
    product_wave(x, y) = cos(2pi * TAU_FAST * Float64(x)) * cos(2pi * TAU_PRODUCT_Y * Float64(y))
    return ((x, y) -> -1.0 - A_HOP * product_wave(Float64(x) + 0.5, y), (x, y) -> -1.0 - A_HOP * product_wave(x, Float64(y) + 0.5))
end
function hopping_two_layer()
    two_layer(x, y) = 0.5 * wave(x, y, TAU_FAST, 0.0) + 0.5 * wave(x, y, TAU_SECOND, THETA_SECOND)
    return ((x, y) -> -1.0 - A_HOP * two_layer(Float64(x) + 0.5, y), (x, y) -> -1.0 - A_HOP * two_layer(x, Float64(y) + 0.5))
end

function solver()
    SCFSettings(purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=30, purification_maxiter=60,
        square_fock_method=:binary_carry, sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6, record_energy=false, stable_iterations=3,
        require_stationarity=false, measure_stationarity=false, detect_two_cycles=true,
        mixing_method=:pulay, pulay_history=4, pulay_warmup=3, pulay_regularization=1e-12,
        pulay_coefficient_limit=8.0, pulay_step_limit=20.0)
end

function cap_case(label, hopping, maxdim)
    model = SquareModel(size=(128, 128), hopping=hopping, interaction=1.0,
        potential=substrate_potential, seed=checkerboard_seed(+2.0), filling=0.5)
    representation = QTTSettings(encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim)
    params = legacy_parameters(model, representation, solver())
    bounds = square_scf_spectral_bounds(params; potential_bounds=(-0.4, 0.4),
        hopping_abs_bounds=(1.1, 1.1), margin=0.5)
    return CaseSpec(label=label, model=model, representation=representation,
        solver=solver(), spectral_bounds=bounds)
end

model_specs = [("hopping_1d_x", hopping_1d_x()),
    ("hopping_rotated_generic_angle", hopping_rotated()),
    ("hopping_product_nonseparable", hopping_product()),
    ("hopping_two_layer_supermoire", hopping_two_layer())]
cases = CaseSpec[]
for (model_label, hopping) in model_specs, cap in (512, 768, 1024)
    push!(cases, cap_case("$(model_label)_chi$(cap)", hopping, cap))
end
campaign = CampaignSpec(name="trilayer_effective_hopping_cap_ladder", cases=cases)
