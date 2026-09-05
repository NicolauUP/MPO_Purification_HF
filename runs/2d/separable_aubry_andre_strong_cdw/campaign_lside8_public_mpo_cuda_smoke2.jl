using MPO_MeanField

# Exact production model and numerical controls, limited to two SCF iterations
# to validate the public CUDA entry point and timing before a long allocation.
Base.include(
    @__MODULE__,
    joinpath(@__DIR__, "campaign_lside8_public_mpo_field_converged_chi512.jl"),
)

production_case = only(campaign.cases)
production_solver = production_case.solver
solver_names = fieldnames(typeof(production_solver))
solver_values = NamedTuple{solver_names}(
    Tuple(getfield(production_solver, name) for name in solver_names),
)
smoke_solver = SCFSettings(; merge(solver_values, (
    maxiter=2,
    stable_iterations=3,
))...)

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside8_public_cuda_smoke2",
    cases=[CaseSpec(
        label=production_case.label,
        model=production_case.model,
        representation=production_case.representation,
        solver=smoke_solver,
        spectral_bounds=production_case.spectral_bounds,
    )],
)
