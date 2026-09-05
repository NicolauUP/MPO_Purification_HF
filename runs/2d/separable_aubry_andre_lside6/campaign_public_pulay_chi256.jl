using MPO_MeanField

# Reuse the physical model and fixed numerical conventions of the public L=6
# campaign. This is a performance/stability pilot, not a second physical case.
Base.include(@__MODULE__, joinpath(@__DIR__, "campaign_public.jl"))

campaign = CampaignSpec(
    name="separable_aubry_andre_lside6_public_mpo_pulay_chi256",
    cases=[lside6_public_case(
        maxdim=256,
        mixing_method=:pulay,
        pulay_history=4,
        pulay_warmup=3,
        pulay_regularization=1e-12,
        pulay_coefficient_limit=8.0,
        pulay_step_limit=20.0,
    )],
)
