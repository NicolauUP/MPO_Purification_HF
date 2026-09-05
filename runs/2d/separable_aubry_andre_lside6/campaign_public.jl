using MPO_MeanField

# Production public-API campaign for the 64×64 separable quasiperiodic MPO
# reference. The explicit enclosure covers the initial seed and conservative
# Hartree/Fock bounds; it is recorded verbatim in the normalized result.
const TAU_AUBRY_ANDRE_LSIDE6_PUBLIC = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_lside6_public(V2::Real)
    amplitude = Float64(V2)
    frequency = TAU_AUBRY_ANDRE_LSIDE6_PUBLIC
    return (
        (x, y) -> -1.0 - amplitude * cos(2π * frequency * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude * cos(2π * frequency * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_lside6_public(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function lside6_public_case(
    ; maxdim::Integer=512,
    cutoff::Real=1e-10,
    mixing_method::Symbol=:linear,
    pulay_history::Integer=4,
    pulay_warmup::Integer=4,
    pulay_regularization::Real=1e-12,
    pulay_coefficient_limit::Real=8.0,
    pulay_step_limit::Real=20.0,
)
    model = SquareModel(
        size=(64, 64),
        hopping=separable_aubry_andre_hopping_lside6_public(0.5),
        interaction=1.0,
        seed=checkerboard_seed_lside6_public(0.1),
        filling=0.5,
    )
    representation = QTTSettings(
        encoding=:interleaved, tci_tol=1e-10, cutoff=cutoff, maxdim=maxdim,
    )
    solver = SCFSettings(
        purification=:sp2,
        mixing=0.5,
        tolerance=0.1,
        maxiter=30,
        purification_maxiter=50,
        square_fock_method=:binary_carry,
        sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6,
        # Per-iteration full MPO energies require extra contractions and are
        # deliberately disabled for this GPU production case. Final
        # observables still include the converged energy.
        record_energy=false,
        stable_iterations=2,
        detect_two_cycles=true,
        mixing_method=mixing_method,
        pulay_history=pulay_history,
        pulay_warmup=pulay_warmup,
        pulay_regularization=pulay_regularization,
        pulay_coefficient_limit=pulay_coefficient_limit,
        pulay_step_limit=pulay_step_limit,
    )
    return CaseSpec(
        label="v2_0p5_u_1_chi_$(maxdim)",
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-8.5, 12.5),
    )
end

campaign = CampaignSpec(
    name="separable_aubry_andre_lside6_public_mpo",
    cases=[lside6_public_case()],
)
