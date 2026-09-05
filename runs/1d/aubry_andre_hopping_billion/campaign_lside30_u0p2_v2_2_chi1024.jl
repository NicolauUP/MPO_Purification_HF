using MPO_MeanField

# Billion-site off-diagonal Aubry--Andre Hartree--Fock gate. The physical
# chain has 2^30 sites but only 30 QTT levels. No scalar onsite potential is
# present: quasiperiodicity enters solely through nearest-neighbour hopping.
const TAU = sqrt(2.0) - 5.0 / 6.0

hopping_modulation(V2::Real) = CosineHopping(
    -1.0, -Float64(V2), 2π * TAU, π * TAU,
)

checkerboard_seed(amplitude::Real) = x ->
    iseven(Int(x)) ? Float64(amplitude) : -Float64(amplitude)

model = ChainModel(
    size=1 << 30,
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
    name="aubry_andre_hopping_lside30_u0p2_v2_2_mpo_sp2_chi1024",
    cases=[CaseSpec(
        label="v2_2_u_0p2_seed_plus0p1_chi1024",
        model=model,
        representation=representation,
        solver=solver,
        # |t_i| <= 3 gives a hopping row bound of 6. The added interval margin
        # covers the bounded nearest-neighbour Hartree and Fock fields.
        spectral_bounds=(-8.0, 8.0),
    )],
)
