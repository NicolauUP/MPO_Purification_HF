using MPO_MeanField

# This campaign deliberately duplicates the fixed-H specification instead of
# including another campaign file: public campaign loading evaluates each file
# in an isolated module where include is unavailable.
#
# The solver is otherwise identical to the L=20 V2=2, U=0.2 control. The
# active source implementation selects direct binary-carry 1D Hartree/Fock
# extraction, avoiding scalar-query TCI after each purification.
const TAU = sqrt(2.0) - 5.0 / 6.0

hopping_modulation(V2::Real) = CosineHopping(
    -1.0, -Float64(V2), 2π * TAU, π * TAU,
)

checkerboard_seed(amplitude::Real) = x ->
    iseven(Int(x)) ? Float64(amplitude) : -Float64(amplitude)

model = ChainModel(
    size=1 << 20,
    hopping=hopping_modulation(2.0),
    interaction=0.2,
    seed=checkerboard_seed(+0.1),
    filling=0.5,
)

representation = QTTSettings(
    tci_tol=1e-10,
    cutoff=1e-10,
    maxdim=1024,
)

solver = SCFSettings(
    purification=:sp2,
    mixing=0.5,
    tolerance=0.1,
    maxiter=20,
    purification_maxiter=60,
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
    name="aubry_andre_hopping_lside20_u0p2_v2_2_mpo_sp2_chi1024_binarycarry",
    cases=[CaseSpec(
        label="v2_2_u_0p2_seed_plus0p1_chi1024",
        model=model,
        representation=representation,
        solver=solver,
        spectral_bounds=(-8.0, 8.0),
    )],
)
