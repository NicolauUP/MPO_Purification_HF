using MPO_MeanField

# Physical two-scale C8 supermoiré hopping field. Each star has a microscopic
# wavelength ~1/tau_center, while their small radial mismatch produces a
# 512-site beat envelope: 1/(tau_fast - tau_slow) = 512.
const C8_TAU_CENTER = sqrt(2.0) / 8.0
const C8_TAU_FAST = C8_TAU_CENTER + 1.0 / 1024.0
const C8_TAU_SLOW = C8_TAU_CENTER - 1.0 / 1024.0
const C8_SUPERMOIRE_PERIOD = 1.0 / (C8_TAU_FAST - C8_TAU_SLOW)
const C8_ANGLES = (0.0, pi / 4, pi / 2, 3pi / 4)
const C8_V_FAST = 0.05
const C8_V_SLOW = 0.05

function c8_star(x::Real, y::Real, tau::Real)
    return sum(cos(2pi * Float64(tau) *
        (Float64(x) * cos(angle) + Float64(y) * sin(angle)))
        for angle in C8_ANGLES) / length(C8_ANGLES)
end

function double_c8_supermoire_modulation(x::Real, y::Real)
    return C8_V_FAST * c8_star(x, y, C8_TAU_FAST) +
           C8_V_SLOW * c8_star(x, y, C8_TAU_SLOW)
end

double_c8_supermoire_hopping = (
    (x, y) -> -1.0 - double_c8_supermoire_modulation(Float64(x) + 0.5, Float64(y)),
    (x, y) -> -1.0 - double_c8_supermoire_modulation(Float64(x), Float64(y) + 0.5),
)

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

campaign = CampaignSpec(
    name="double_c8_lside11_supermoire_mpo_sp2",
    cases=[CaseSpec(
        label="double_c8_tau_fast_slow_l2048_u2_v005x2_chi256",
        model=SquareModel(
            size=(2048, 2048), hopping=double_c8_supermoire_hopping,
            interaction=2.0, seed=checkerboard_seed(+2.0), filling=0.5,
        ),
        representation=QTTSettings(
            encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=256,
        ),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=60,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            record_energy=false, stable_iterations=3, require_stationarity=false,
            measure_stationarity=false, detect_two_cycles=true, mixing_method=:linear,
        ),
        spectral_bounds=(-8.5, 12.5),
    )],
)
