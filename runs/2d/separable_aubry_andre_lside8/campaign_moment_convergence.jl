using MPO_MeanField

# Fixed 256x256 representative-point campaign for selecting a common KPM
# resolution before a broad interacting (U, V2) phase map.
const TAU_MOMENT_CONVERGENCE = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_moment_convergence(V2::Real)
    amplitude = Float64(V2)
    return (
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_MOMENT_CONVERGENCE * (Float64(x) + 0.5)),
        (x, y) -> -1.0 - amplitude *
            cos(2π * TAU_MOMENT_CONVERGENCE * (Float64(y) + 0.5)),
    )
end

checkerboard_seed_moment_convergence(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

const MOMENT_CASES = (
    (label="v2_0_u_0p5", V2=0.0, U=0.5),
    (label="v2_0p5_u_0p5", V2=0.5, U=0.5),
    (label="v2_2_u_0p3", V2=2.0, U=0.3),
)

case_index = parse(Int, get(ENV, "KPM_MOMENT_CASE", "1"))
1 <= case_index <= length(MOMENT_CASES) ||
    error("KPM_MOMENT_CASE must lie in 1:$(length(MOMENT_CASES))")
seed_mode = Symbol(lowercase(get(ENV, "KPM_MOMENT_SEED", "checkerboard")))
seed_mode in (:checkerboard, :uniform) ||
    error("KPM_MOMENT_SEED must be checkerboard or uniform")
case = MOMENT_CASES[case_index]

# The uniform seed preserves the unbroken-symmetry basin. The checkerboard
# seed tests whether a symmetry-broken Hartree--Fock solution is selected.
seed = seed_mode == :checkerboard ?
    checkerboard_seed_moment_convergence(0.1) : nothing

params = ParametersSquare(
    L=16,
    t=separable_aubry_andre_hopping_moment_convergence(case.V2),
    U=case.U,
    W=nothing,
    S=seed,
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
    name="separable_aubry_andre_lside8_moment_convergence",
    runs=[(
        label="$(case.label)_$(seed_mode)",
        params=params,
        purification_method=:kpm,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
