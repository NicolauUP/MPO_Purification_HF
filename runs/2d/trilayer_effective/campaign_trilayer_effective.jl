using MPO_MeanField

# Effective single-band model for two active layers on an aligned substrate.
# There is no explicit layer index: the coupling of the active layers is
# represented by one incommensurate hopping field, while the substrate enters
# as an onsite modulation. The substrate is not rotated, but has a different
# lattice spacing (a different modulation wave number).
const TAU_HOPPING = sqrt(2.0) - 5.0 / 6.0
const THETA_HOPPING = pi / 4
const HOPPING_AMPLITUDE = 0.10
const TAU_SUBSTRATE = 1.0 / 32.0
const THETA_SUBSTRATE = 0.0
const SUBSTRATE_AMPLITUDE = 0.20

function incommensurate_wave(x::Real, y::Real, theta::Real, tau::Real)
    c, s = cos(theta), sin(theta)
    return cos(2pi * Float64(tau) * (Float64(x) * c + Float64(y) * s))
end

function effective_hopping()
    return (
        (x, y) -> -1.0 - HOPPING_AMPLITUDE * incommensurate_wave(
            Float64(x) + 0.5, Float64(y), THETA_HOPPING, TAU_HOPPING),
        (x, y) -> -1.0 - HOPPING_AMPLITUDE * incommensurate_wave(
            Float64(x), Float64(y) + 0.5, THETA_HOPPING, TAU_HOPPING),
    )
end

function substrate_potential(x::Real, y::Real)
    return SUBSTRATE_AMPLITUDE * incommensurate_wave(
        Float64(x), Float64(y), THETA_SUBSTRATE, TAU_SUBSTRATE)
end

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function dense_solver()
    return SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=1e-8, maxiter=60,
        purification_maxiter=60, square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
        record_energy=true, stable_iterations=3, require_stationarity=true,
        measure_stationarity=true, detect_two_cycles=true, mixing_method=:pulay,
        pulay_history=4, pulay_warmup=3, pulay_regularization=1e-12,
        pulay_coefficient_limit=8.0, pulay_step_limit=20.0,
    )
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

function trilayer_case(side::Int, maxdim::Int, label::String, solver::SCFSettings)
    return CaseSpec(
        label=label,
        model=SquareModel(
            size=(side, side), hopping=effective_hopping(), interaction=1.0,
            potential=substrate_potential, seed=checkerboard_seed(+2.0), filling=0.5,
        ),
        representation=QTTSettings(
            encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim,
        ),
        solver=solver,
        # Conservative bound for |t_x|+|t_y|+|v|+U at this first test.
        spectral_bounds=(-12.0, 12.0),
    )
end

campaign = CampaignSpec(
    name="trilayer_effective_rotated_hopping_substrate",
    cases=[
        trilayer_case(64, 512, "trilayer_effective_l64_dense", dense_solver()),
        trilayer_case(128, 256, "trilayer_effective_l128_mpo_chi256", mpo_solver()),
        trilayer_case(128, 512, "trilayer_effective_l128_mpo_chi512", mpo_solver()),
    ],
)
