using MPO_MeanField

# Conference production candidate: a 256 x 256 weakly quasiperiodic hopping
# model in the strong-U checkerboard regime. The two cases select the two
# symmetry-related broken-sublattice branches. The seed is initialization only.
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

function strong_cdw_case(label::String, seed_amplitude::Real)
    model = SquareModel(
        size=(256, 256),
        hopping=separable_aubry_andre_hopping(0.1),
        interaction=2.0,
        seed=checkerboard_seed(seed_amplitude),
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
        maxiter=60,
        purification_maxiter=60,
        square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6,
        record_energy=false,
        stable_iterations=2,
        detect_two_cycles=true,
        mixing_method=:pulay,
        pulay_history=4,
        pulay_warmup=3,
        pulay_regularization=1e-12,
        pulay_coefficient_limit=8.0,
        pulay_step_limit=20.0,
    )
    # Hopping radius <= 4.4. The nearest-neighbour Hartree and Fock fields
    # contribute at most 8|U|=16 in this implementation. Include |S|=2 for
    # iteration one and a 0.5 numerical margin: 4.4 + 16 + 2 + 0.5 < 23.
    return CaseSpec(
        label=label,
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-23.0, 23.0),
    )
end

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside8_mpo_pulay_chi512",
    cases=[
        strong_cdw_case("v2_0p1_u_2_seed_plus2", +2.0),
        strong_cdw_case("v2_0p1_u_2_seed_minus2", -2.0),
    ],
)
