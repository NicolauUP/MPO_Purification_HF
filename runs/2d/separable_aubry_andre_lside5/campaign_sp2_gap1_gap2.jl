using MPO_MeanField

# Fixed 32x32 clean-square SP2 controls between the failed gap-0.2 case and
# the converged gap-4 case. For H=K+mΓ, the total gap is Δ=2|m|.
constant_square_hopping() = (
    (x, y) -> -1.0,
    (x, y) -> -1.0,
)

checkerboard_mass(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function sp2_gap_case(gap::Real)
    gap_value = Float64(gap)
    gap_value > 0 || throw(ArgumentError("gap must be positive"))
    mass = gap_value / 2
    label = "gap_" * replace(string(gap_value), "." => "p")
    params = ParametersSquare(
        L=10,
        t=constant_square_hopping(),
        U=1.0,
        W=nothing,
        S=checkerboard_mass(mass),
        tci_tol=1e-10,
        itensors_tol=1e-14,
        itensors_maxdim=256,
        density=0.5,
        purification_steps=50,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=1,
    )
    return (
        label=label,
        target_gap=gap_value,
        checkerboard_mass=mass,
        params=params,
        spectral_bounds=(-4.25, 4.25),
        purification_method=:sp2,
    )
end

campaign = (
    name="square_lside5_sp2_gap1_gap2",
    runs=[
        sp2_gap_case(1.0),
        sp2_gap_case(2.0),
    ],
)
