using MPO_MeanField

# A 64x64 open square in interleaved QTT ordering has 12 quantics bits.
const TAU_AUBRY_ANDRE_LSIDE6 = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_lside6(V2::Real)
    amplitude = Float64(V2)
    frequency = TAU_AUBRY_ANDRE_LSIDE6
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * frequency * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * frequency * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_lside6(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function reference_case_lside6()
    V2 = 0.5
    params = ParametersSquare(
        L=12,
        t=separable_aubry_andre_hopping_lside6(V2),
        U=1.0,
        W=nothing,
        S=checkerboard_seed_lside6(0.1),
        tci_tol=1e-10,
        itensors_tol=1e-10,
        itensors_maxdim=256,
        density=0.5,
        purification_steps=50,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=30,
    )
    # This campaign is used by the fixed-H diagnostics. H0 + S has hopping
    # row sum at most 4(1 + V2) = 6 and |S| = 0.1. The extra 0.25 protects
    # against construction/truncation error without including absent HF fields.
    bounds = (-6.35, 6.35)
    return (
        label="v2_0p5_u_1",
        params=params,
        spectral_bounds=bounds,
        purification_method=:sp2,
        verify_spectral_bounds=false,
        verbose=:all,
    )
end

campaign = (
    name="separable_aubry_andre_lside6_seed0p1",
    runs=[reference_case_lside6()],
)
