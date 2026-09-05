using MPO_MeanField

# A 256x256 open square in interleaved QTT ordering has 16 quantics bits.
const TAU_AUBRY_ANDRE_LSIDE8 = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_lside8(V2::Real)
    amplitude = Float64(V2)
    frequency = TAU_AUBRY_ANDRE_LSIDE8
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * frequency * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * frequency * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_lside8(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=16,
    t=separable_aubry_andre_hopping_lside8(0.5),
    U=1.0,
    W=nothing,
    S=checkerboard_seed_lside8(0.1),
    tci_tol=1e-10,
    itensors_tol=1e-10,
    itensors_maxdim=256,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=30,
)

campaign = (
    name="separable_aubry_andre_lside8_seed0p1",
    runs=[(
        label="v2_0p5_u_1",
        params=params,
        # Analytically safe fixed-H enclosure: hopping row sum <= 6,
        # checkerboard seed <= 0.1, plus 0.25 margin.
        spectral_bounds=(-6.35, 6.35),
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)

