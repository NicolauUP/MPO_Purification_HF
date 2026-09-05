using MPO_MeanField

# 1024×1024 open square. The irrational is close to the previous 1D/2D
# Aubry--André choice and remains separable only for this scale calibration.
const TAU_AUBRY_ANDRE_LSIDE10 = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_lside10(V2::Real)
    amplitude = Float64(V2)
    frequency = TAU_AUBRY_ANDRE_LSIDE10
    return (
        (x, y) -> -1.0 - amplitude * cos(2π * frequency * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude * cos(2π * frequency * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_lside10(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=20,
    t=separable_aubry_andre_hopping_lside10(0.5),
    U=1.0,
    W=nothing,
    S=checkerboard_seed_lside10(0.1),
    tci_tol=1e-10,
    itensors_tol=1e-10,
    itensors_maxdim=1024,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=30,
)

campaign = (
    name="separable_aubry_andre_lside10_seed0p1",
    runs=[(
        label="v2_0p5_u_1",
        params=params,
        # Initial seeded one-body field: hopping row sum <= 6, |S| <= 0.1.
        spectral_bounds=(-6.35, 6.35),
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
