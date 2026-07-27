using MPO_MeanField

# Fixed 32x32 clean square used only to test canonical SP2 at a large,
# analytically controlled one-particle gap. For nearest-neighbour hopping,
# H = K + mΓ has E = ±sqrt(ε_K^2 + m^2); m=2 therefore gives gap Δ=4.
constant_square_hopping() = (
    (x, y) -> -1.0,
    (x, y) -> -1.0,
)

checkerboard_mass(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

params = ParametersSquare(
    L=10,
    t=constant_square_hopping(),
    U=1.0,
    W=nothing,
    S=checkerboard_mass(2.0),
    tci_tol=1e-10,
    itensors_tol=1e-14,
    itensors_maxdim=256,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=1,
)

campaign = (
    name="square_lside5_sp2_gap4",
    runs=[(
        label="checkerboard_mass_2_gap_4",
        target_gap=4.0,
        checkerboard_mass=2.0,
        params=params,
        spectral_bounds=(-4.75, 4.75),
        purification_method=:sp2,
    )],
)
