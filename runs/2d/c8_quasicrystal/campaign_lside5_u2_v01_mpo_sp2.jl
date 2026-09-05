using MPO_MeanField

# MPO--SP2 counterpart of the dense C8 reference. The four cosine waves
# generate the complete ±q octagonal star. The modulation is evaluated on
# bond midpoints, so the x and y hoppings are both sampled from one genuinely
# non-separable scalar field.
const C8_TAU = sqrt(2.0) - 5.0 / 6.0
const C8_ANGLES = (0.0, pi / 4, pi / 2, 3pi / 4)

function c8_modulation(x::Real, y::Real)
    return sum(cos(2pi * C8_TAU *
        (Float64(x) * cos(angle) + Float64(y) * sin(angle)))
        for angle in C8_ANGLES) / length(C8_ANGLES)
end

c8_hopping(V8::Real) = (
    (x, y) -> -1.0 - Float64(V8) * c8_modulation(Float64(x) + 0.5, Float64(y)),
    (x, y) -> -1.0 - Float64(V8) * c8_modulation(Float64(x), Float64(y) + 0.5),
)

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function c8_sp2_case(maxdim::Integer)
    return CaseSpec(
        label="c8_v8_0p1_u_2_seed_plus2_chi_$(maxdim)",
        model=SquareModel(
            size=(32, 32), hopping=c8_hopping(0.1), interaction=2.0,
            seed=checkerboard_seed(+2.0), filling=0.5,
        ),
        representation=QTTSettings(
            encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=maxdim,
        ),
        solver=SCFSettings(
            purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=60,
            purification_maxiter=60, square_fock_method=:binary_carry,
            sp2_idempotency_tolerance=2e-4, sp2_relative_trace_tolerance=1e-6,
            # Avoid an additional MPO product every SCF iteration. Energies
            # are computed only in the final host-side observable audit.
            record_energy=false, stable_iterations=3, require_stationarity=false,
            measure_stationarity=false, detect_two_cycles=true, mixing_method=:linear,
        ),
        # Conservative explicit enclosure for the ±2 checkerboard initial
        # field. This remains fixed throughout SP2, as required by its map.
        spectral_bounds=(-8.5, 12.5),
    )
end

campaign = CampaignSpec(
    name="c8_quasicrystal_lside5_u2_v01_mpo_sp2_caps",
    cases=[c8_sp2_case(maxdim) for maxdim in (256, 512, 1024)],
)
