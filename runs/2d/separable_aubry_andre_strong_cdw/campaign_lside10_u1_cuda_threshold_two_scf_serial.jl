# Reuse the exact two-SCF physics and numerical settings from the threshold
# campaign, while using a fresh output namespace and only the 768/896 checks.
Base.include(@__MODULE__, joinpath(@__DIR__, "campaign_lside10_u1_cuda_threshold_two_scf.jl"))

campaign = CampaignSpec(
    name="separable_aubry_andre_strong_cdw_lside10_u1_cuda_threshold_two_scf_serial",
    cases=campaign.cases[1:2],
)
