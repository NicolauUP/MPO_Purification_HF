using MPO_MeanField

# Open 1D nearest-neighbour Hartree--Fock reference cases for comparison with
# external ED results. The coordinate provided to `t` is zero-based, so this
# implements t(x) = -1 - V2*cos(2π*τ*(x + 1/2)) exactly as specified.
const TAU_AUBRY_ANDRE = sqrt(2.0) - 5.0 / 6.0

function aubry_andre_hopping(V2::Real)
    x -> -1.0 - Float64(V2) * cos(2π * TAU_AUBRY_ANDRE * (Float64(x) + 0.5))
end

# A temporary positive checkerboard Hartree seed. `S` initializes only the
# first SCF Hamiltonian; it is replaced by the extracted mean field after the
# first iteration and is not an external physical onsite potential.
checkerboard_seed(amplitude::Real) = x -> iseven(Int(x)) ? Float64(amplitude) : -Float64(amplitude)

# For an open chain, ||H0||∞ <= 2*max|t|. With a physical density matrix,
# each nearest-neighbour Hartree and real-exchange row contribution is bounded
# by 2|U|, giving the conservative SCF interval ±(2*(1+V2)+4|U|+0.5).
function conservative_scf_bounds(V2::Real, U::Real)
    radius = 2.0 * (1.0 + abs(Float64(V2))) + 4.0 * abs(Float64(U)) + 0.5
    return (-radius, radius)
end

function reference_case(label::String, V2::Real, U::Real)
    params = Parameters1D(
        L=10,
        t=aubry_andre_hopping(V2),
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
    return (
        label=label,
        params=params,
        spectral_bounds=conservative_scf_bounds(V2, U),
        purification_method=:sp2,
        verify_spectral_bounds=false,
        verbose=:all,
    )
end

campaign = (
    name="aubry_andre_nn_l10_seed0p1",
    runs=[
        reference_case("v2_0_u_1", 0.0, 1.0),
        reference_case("v2_0p5_u_1", 0.5, 1.0),
        reference_case("v2_2_u_0p2", 2.0, 0.2),
    ],
)
