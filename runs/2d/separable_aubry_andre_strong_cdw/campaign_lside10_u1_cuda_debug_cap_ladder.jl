using MPO_MeanField

# CUDA-debug cap ladder.  Physics and solver settings match the prior L=10 U=1
# test; only the caps and campaign path differ so failed attempts never share
# a result directory with the production 512 run.
const TAU_STRONG_CDW_L10_U1_DEBUG = sqrt(2.0) - 5.0 / 6.0

separable_aubry_andre_hopping_l10_u1_debug(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L10_U1_DEBUG * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L10_U1_DEBUG * (Float64(y) + 0.5)),
)
checkerboard_seed_l10_u1_debug(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function u1_debug_case(maxdim::Int)
    CaseSpec(
        label="v2_0p1_u_1_seed_plus2_chi_$(maxdim)",
        model=SquareModel(
            size=(1024, 1024), hopping=separable_aubry_andre_hopping_l10_u1_debug(0.1),
            interaction=1.0, seed=checkerboard_seed_l10_u1_debug(+2.0), filling=0.5,
        ),
        representation=QTTSettings(encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=20,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            record_energy=false, stable_iterations=3, require_stationarity=false,
            measure_stationarity=false, detect_two_cycles=true,
            mixing_method=:pulay, pulay_history=4, pulay_warmup=3,
            pulay_regularization=1e-12, pulay_coefficient_limit=8.0, pulay_step_limit=20.0,
        ),
        spectral_bounds=(-23.0, 23.0),
    )
end

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside10_u1_cuda_debug_cap_ladder",
    cases=[u1_debug_case(640), u1_debug_case(768), u1_debug_case(1024)],
)
