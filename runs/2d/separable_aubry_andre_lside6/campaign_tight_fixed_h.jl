include(joinpath(@__DIR__, "campaign.jl"))

base_spec = campaign.runs[1]
campaign = (
    name="separable_aubry_andre_lside6_tight_fixed_h",
    runs=[merge(base_spec, (spectral_bounds=(-4.1, 4.1),))],
)
