using MPO_MeanField

# `side_level=5` means a 2^5 by 2^5 open square. ParametersSquare stores the
# total number of interleaved quantics bits, so L=2*side_level=10.
const TAU_AUBRY_ANDRE_2D = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping(V2::Real)
    amplitude = Float64(V2)
    tx = (x, y) -> -1.0 - amplitude * cos(2π * TAU_AUBRY_ANDRE_2D * (Float64(x) + 0.5))
    ty = (x, y) -> -1.0 - amplitude * cos(2π * TAU_AUBRY_ANDRE_2D * (Float64(y) + 0.5))
    return (tx, ty)
end

# This selects a checkerboard HF branch but is not a physical onsite term.
checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function reference_case(label::String, V2::Real, U::Real)
    params = ParametersSquare(
        L=10,
        t=separable_aubry_andre_hopping(V2),
        U=Float64(U),
        W=nothing,
        S=checkerboard_seed(0.1),
        tci_tol=1e-10,
        itensors_tol=1e-14,
        itensors_maxdim=256,
        density=0.5,
        purification_steps=50,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=30,
    )
    bounds = square_scf_spectral_bounds(
        params;
        hopping_abs_bounds=(1.0 + abs(Float64(V2)), 1.0 + abs(Float64(V2))),
        margin=0.5,
    )
    return (
        label=label,
        params=params,
        spectral_bounds=bounds,
        purification_method=:sp2,
        verify_spectral_bounds=false,
        verbose=:all,
    )
end

campaign = (
    name="separable_aubry_andre_lside5_seed0p1",
    runs=[
        reference_case("v2_0_u_1", 0.0, 1.0),
        reference_case("v2_0p5_u_1", 0.5, 1.0),
        reference_case("v2_2_u_0p2", 2.0, 0.2),
    ],
)
