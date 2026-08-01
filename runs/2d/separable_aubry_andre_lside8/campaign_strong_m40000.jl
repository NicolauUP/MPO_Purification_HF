using MPO_MeanField

# High-degree production SCF for the strong quasiperiodic CDW point. The
# +checkerboard seed selects the symmetry-broken branch already established
# by the seed-stability study; this run tests KPM-SCF convergence at M=40000.
const TAU_STRONG_M40000 = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_strong_m40000(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_STRONG_M40000 * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_STRONG_M40000 * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_strong_m40000(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=16,
    t=separable_aubry_andre_hopping_strong_m40000(2.0),
    U=0.3,
    W=nothing,
    S=checkerboard_seed_strong_m40000(0.1),
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
    name="separable_aubry_andre_lside8_strong_m40000",
    runs=[(
        label="v2_2_u_0p3_checkerboard_plus",
        params=params,
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
