using MPO_MeanField

# Seed-stability study for the 256x256 open square. The hopping modulation is
# separable and quasiperiodic along each coordinate. The three seeds contain a
# positive, negative, or weak explicit checkerboard component; an exactly
# uniform seed is deliberately excluded because it remains in the
# symmetry-preserving SCF sector by construction.
const TAU_SEED_STABILITY = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_seed_stability(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_SEED_STABILITY * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_SEED_STABILITY * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_seed_stability(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

const SEED_STABILITY_CASES = (
    # The moment orders come from the completed representative-point ladder.
    (label="v2_0_u_0p5", V2=0.0, U=0.5, moments=1200, audit_moments=1600),
    (label="v2_0p5_u_0p5", V2=0.5, U=0.5, moments=1600, audit_moments=2000),
    (label="v2_2_u_0p3", V2=2.0, U=0.3, moments=2400, audit_moments=3200),
)
const SEED_STABILITY_SEEDS = (
    (label="checkerboard_plus", amplitude=0.1),
    (label="checkerboard_minus", amplitude=-0.1),
    (label="checkerboard_weak", amplitude=1e-3),
)

case_index = parse(Int, get(ENV, "KPM_SEED_STABILITY_CASE", "1"))
seed_index = parse(Int, get(ENV, "KPM_SEED_STABILITY_SEED", "1"))
1 <= case_index <= length(SEED_STABILITY_CASES) ||
    error("KPM_SEED_STABILITY_CASE must lie in 1:$(length(SEED_STABILITY_CASES))")
1 <= seed_index <= length(SEED_STABILITY_SEEDS) ||
    error("KPM_SEED_STABILITY_SEED must lie in 1:$(length(SEED_STABILITY_SEEDS))")
case = SEED_STABILITY_CASES[case_index]
seed_spec = SEED_STABILITY_SEEDS[seed_index]

params = ParametersSquare(
    L=16,
    t=separable_aubry_andre_hopping_seed_stability(case.V2),
    U=case.U,
    W=nothing,
    S=checkerboard_seed_seed_stability(seed_spec.amplitude),
    tci_tol=1e-10,
    itensors_tol=1e-10,
    itensors_maxdim=256,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=60,
)

campaign = (
    name="separable_aubry_andre_lside8_seed_stability",
    runs=[(
        label="$(case.label)_$(seed_spec.label)",
        params=params,
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
