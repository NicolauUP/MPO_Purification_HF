using MPO_MeanField

# Public-API equivalent of campaign.jl. Keep the original source untouched
# while this campaign exercises the generic CPU MPO/dense entry point.
const TAU_AUBRY_ANDRE_2D_PUBLIC = sqrt(2.0) - 5.0 / 6.0

function separable_aubry_andre_hopping_public(V2::Real)
    amplitude = Float64(V2)
    tx = (x, y) -> -1.0 - amplitude * cos(2π * TAU_AUBRY_ANDRE_2D_PUBLIC * (Float64(x) + 0.5))
    ty = (x, y) -> -1.0 - amplitude * cos(2π * TAU_AUBRY_ANDRE_2D_PUBLIC * (Float64(y) + 0.5))
    return (tx, ty)
end

checkerboard_seed_public(amplitude::Real) = (x, y) ->
    iseven(Int(x) + Int(y)) ? Float64(amplitude) : -Float64(amplitude)

function public_case(label::String, V2::Real, U::Real)
    model = SquareModel(
        size=(32, 32),
        hopping=separable_aubry_andre_hopping_public(V2),
        interaction=Float64(U),
        seed=checkerboard_seed_public(0.1),
        filling=0.5,
    )
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-14, maxdim=256)
    solver = SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1,
        maxiter=30, purification_maxiter=50,
    )
    params = legacy_parameters(model, representation, solver)
    bounds = square_scf_spectral_bounds(
        params;
        hopping_abs_bounds=(1.0 + abs(Float64(V2)), 1.0 + abs(Float64(V2))),
        margin=0.5,
    )
    return CaseSpec(
        label=label, model=model, representation=representation, solver=solver,
        spectral_bounds=bounds,
    )
end

campaign = CampaignSpec(
    name="separable_aubry_andre_lside5_seed0p1",
    cases=[
        public_case("v2_0_u_1", 0.0, 1.0),
        public_case("v2_0p5_u_1", 0.5, 1.0),
        public_case("v2_2_u_0p2", 2.0, 0.2),
    ],
)
