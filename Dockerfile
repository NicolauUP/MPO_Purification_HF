FROM --platform=linux/amd64 julia:1.12.6-bookworm

WORKDIR /opt/MPO_Purification_HF

# Install the locked dependency environment before copying source files so
# Docker can reuse this layer when only Julia source or tests change.
COPY Project.toml Manifest.toml ./

ENV JULIA_DEPOT_PATH=/opt/julia_depot

RUN julia --project=. -e \
    'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY . .

# Default container action: run the maintained package test entry point.
CMD ["julia", "--project=.", "test/runtests.jl"]
