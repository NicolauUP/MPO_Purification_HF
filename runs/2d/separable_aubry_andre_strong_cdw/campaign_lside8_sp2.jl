using MPO_MeanField

# 256 x 256 fixed-H SP2 scaling case. The seed opens the initial one-particle
# gap; the weak hopping modulation is the physical quasiperiodic structure.
const TAU_STRONG_CDW = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude * cos(2π * TAU_STRONG_CDW * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude * cos(2π * TAU_STRONG_CDW * (Float64(y) + 0.5)),
    )
end

checkerboard_seed(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function strong_cdw_case(label::String, seed_amplitude::Real)
    params = ParametersSquare(
        L=16, # 256 x 256 sites
        t=separable_aubry_andre_hopping(0.1), U=2.0,
        W=nothing, S=checkerboard_seed(seed_amplitude),
        tci_tol=1e-10, itensors_tol=1e-10, itensors_maxdim=512,
        density=0.5, purification_steps=50,
        scf_mixing=0.5, scf_tol=0.01, scf_max_iterations=80,
    )
    return (
        label=label, seed_amplitude=Float64(seed_amplitude), params=params,
        spectral_bounds=(-6.9, 6.9), purification_method=:sp2,
        verify_spectral_bounds=false, verbose=:all,
    )
end

campaign = (
    name="separable_aubry_andre_strong_cdw_lside8_sp2",
    runs=[
        strong_cdw_case("v2_0p1_u_2_seed_plus2", +2.0),
        strong_cdw_case("v2_0p1_u_2_seed_minus2", -2.0),
    ],
)
