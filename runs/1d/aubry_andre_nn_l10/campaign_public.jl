using MPO_MeanField

# Public-API equivalent of campaign.jl. The legacy source remains unchanged as
# an expert/campaign regression reference during the migration period.
const TAU_AUBRY_ANDRE_PUBLIC = sqrt(2.0) - 5.0 / 6.0

aubry_andre_hopping_public(V2::Real) = x ->
    -1.0 - Float64(V2) * cos(2π * TAU_AUBRY_ANDRE_PUBLIC * (Float64(x) + 0.5))

checkerboard_seed_public(amplitude::Real) = x ->
    iseven(Int(x)) ? Float64(amplitude) : -Float64(amplitude)

function public_case(label::String, V2::Real, U::Real)
    model = ChainModel(
        size=1024,
        hopping=aubry_andre_hopping_public(V2),
        interaction=Float64(U),
        seed=checkerboard_seed_public(0.1),
        filling=0.5,
    )
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-14, maxdim=256)
    solver = SCFSettings(
        purification=:sp2, mixing=0.5, tolerance=0.1,
        maxiter=30, purification_maxiter=50,
    )
    radius = 2.0 * (1.0 + abs(Float64(V2))) + 4.0 * abs(Float64(U)) + 0.5
    return CaseSpec(
        label=label, model=model, representation=representation, solver=solver,
        spectral_bounds=(-radius, radius),
    )
end

campaign = CampaignSpec(
    name="aubry_andre_nn_l10_seed0p1",
    cases=[
        public_case("v2_0_u_1", 0.0, 1.0),
        public_case("v2_0p5_u_1", 0.5, 1.0),
        public_case("v2_2_u_0p2", 2.0, 0.2),
    ],
)
