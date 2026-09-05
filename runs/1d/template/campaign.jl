# Compatibility pointer for the old template location. New work should copy
# runs/template/campaign.jl and use bin/mpohf.jl. Keeping this include allows
# existing documentation links to fail safely into the canonical schema.
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "template", "campaign.jl"))
