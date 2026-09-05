using MPO_MeanField

# Small-system MPO cap convergence against campaign_lside5_u1_dense_reference.
const TAU_STRONG_CDW_L5_U1_CAPS = sqrt(2.0) - 5.0 / 6.0

separable_aubry_andre_hopping_l5_u1_caps(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L5_U1_CAPS * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L5_U1_CAPS * (Float64(y) + 0.5)),
)
checkerboard_seed_l5_u1_caps(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function u1_small_mpo_case(maxdim::Int)
    CaseSpec(
        label="v2_0p1_u_1_seed_plus2_chi_$(maxdim)",
        model=SquareModel(
            size=(32, 32), hopping=separable_aubry_andre_hopping_l5_u1_caps(0.1),
            interaction=1.0, seed=checkerboard_seed_l5_u1_caps(+2.0), filling=0.5,
        ),
        representation=QTTSettings(encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=40,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            record_energy=true, stable_iterations=3, require_stationarity=false,
            measure_stationarity=false, detect_two_cycles=true,
            mixing_method=:pulay, pulay_history=4, pulay_warmup=3,
            pulay_regularization=1e-12, pulay_coefficient_limit=8.0, pulay_step_limit=20.0,
        ),
        spectral_bounds=(-23.0, 23.0),
    )
end

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside5_u1_mpo_cap_ladder",
    cases=[u1_small_mpo_case(256), u1_small_mpo_case(512), u1_small_mpo_case(768), u1_small_mpo_case(1024)],
)
