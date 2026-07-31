using MPO_MeanField

const TAU_AA_LSIDE5_STRICT = sqrt(2.0) - 5.0 / 6.0

function hopping_lside5_strict(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_AA_LSIDE5_STRICT * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_AA_LSIDE5_STRICT * (Float64(y) + 0.5)),
    )
end

seed_lside5_strict(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=10,
    t=hopping_lside5_strict(0.5),
    U=1.0,
    W=nothing,
    S=seed_lside5_strict(0.1),
    tci_tol=1e-10,
    itensors_tol=1e-10,
    itensors_maxdim=256,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.02,
    scf_max_iterations=50,
)

campaign = (
    name="separable_aubry_andre_lside5_v2_0p5_strict",
    runs=[(
        label="v2_0p5_u_1",
        params=params,
        spectral_bounds=(-6.35, 6.35),
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
