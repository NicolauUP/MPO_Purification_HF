using MPO_MeanField

# 1024 x 1024 fixed-initial-field scaling test. This deliberately preserves
# the physical model of the successful 256 x 256 strong-CDW branch: only the
# QTT lattice depth changes. It is not a self-consistent HF calculation.
const TAU_STRONG_CDW_L10 = sqrt(2.0) - 5.0 / 6.0

separable_aubry_andre_hopping_l10(V2::Real) = (
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L10 * (Float64(x) + 0.5)),
    (x, y) -> -1.0 - Float64(V2) * cos(2π * TAU_STRONG_CDW_L10 * (Float64(y) + 0.5)),
)
checkerboard_seed_l10(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=20, # 1024 x 1024 sites
    t=separable_aubry_andre_hopping_l10(0.1), U=2.0,
    W=nothing, S=checkerboard_seed_l10(+2.0),
    tci_tol=1e-10, itensors_tol=1e-10, itensors_maxdim=512,
    density=0.5, purification_steps=50,
    scf_mixing=0.5, scf_tol=0.01, scf_max_iterations=80,
)

campaign = (
    name="separable_aubry_andre_strong_cdw_lside10_fixed_sp2",
    runs=[(
        label="v2_0p1_u_2_seed_plus2",
        params=params,
        spectral_bounds=(-6.9, 6.9),
        purification_method=:sp2,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
