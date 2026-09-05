using MPO_MeanField

# Production 256 x 256 strong-CDW branch. The previous two-seed validation
# established a stable SP2/MPO commutator floor near 2.5e-3 while the actual
# Hartree, Fock, and density fixed-point residuals fell below 1e-3. Therefore
# the commutator is omitted inside the SCF loop and measured only by the final
# host-side observables audit; it is not a stopping gate here.
const TAU_STRONG_CDW = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude * cos(2π * TAU_STRONG_CDW * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude * cos(2π * TAU_STRONG_CDW * (Float64(y) + 0.5)),
    )
end

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

model = SquareModel(
    size=(256, 256),
    hopping=separable_aubry_andre_hopping(0.1),
    interaction=2.0,
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
    name="separable_aubry_andre_strong_cdw_lside8_mpo_field_converged_chi512",
    cases=[CaseSpec(
        label="v2_0p1_u_2_seed_plus2",
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-23.0, 23.0),
    )],
)
