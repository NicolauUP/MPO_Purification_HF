using MPO_MeanField

# One controlled refinement of task 1 in campaign.jl. This keeps the same
# physical model, seed, spectral interval, and SCF policy while raising only
# the MPO bond-dimension cap from 256 to 512.
const TAU_AUBRY_ANDRE_2D = sqrt(2.0) - 5.0 / 6.0

tx = (x, y) -> -1.0
ty = (x, y) -> -1.0
checkerboard_seed = (x, y) -> iseven(Int(x) + Int(y)) ? 0.1 : -0.1

params = ParametersSquare(
    L=10,
    t=(tx, ty),
    U=1.0,
    W=nothing,
    S=checkerboard_seed,
    tci_tol=1e-10,
    itensors_tol=1e-14,
    itensors_maxdim=512,
    density=0.5,
    purification_steps=50,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=30,
)

campaign = (
    name="separable_aubry_andre_lside5_seed0p1_chi512",
    runs=[(
        label="v2_0_u_1",
        params=params,
        spectral_bounds=square_scf_spectral_bounds(
            params; hopping_abs_bounds=(1.0, 1.0), margin=0.5,
        ),
        purification_method=:sp2,
        verify_spectral_bounds=false,
        verbose=:all,
    )],
)
