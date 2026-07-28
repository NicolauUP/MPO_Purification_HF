# CPU square mean-field profiling

These diagnostics prepare the same fixed gapped checkerboard density with SP2
and then measure mean-field extraction on CPU. They do not run SCF.

`component_profile_2d.jl` times the current untruncated adjacency Hartree path,
horizontal and vertical TCI Fock paths, their assembly, and the existing
binary-carry Fock prototypes. The binary-carry rows report their relative MPO
error against the current TCI convention. `fock_direct_validation.csv`
independently compares both extractors with direct density-matrix
measurements on interior, boundary-adjacent, and forbidden wraparound bonds.

`hartree_truncation_sweep_2d.jl` compares the current untruncated
`U*A_nn*n` application with truncation applied during that operation. It
sweeps cutoffs `1e-12,1e-10,1e-8` and maximum bond dimensions
`64,128,256,384,512`, recording runtime, allocation, resulting rank, relative
MPO error, and five direct Hartree probes.

Submit both independent CPU jobs:

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks

PROFILE_JOB="$(
  sbatch --parsable \
    --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
    benchmark/mean_field_profile_2d/component_profile_2d.slurm
)"

SWEEP_JOB="$(
  sbatch --parsable \
    --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
    benchmark/mean_field_profile_2d/hartree_truncation_sweep_2d.slurm
)"

echo "component profile job: $PROFILE_JOB"
echo "Hartree sweep job: $SWEEP_JOB"
```

The jobs request CPU memory through 20 allocated cores but deliberately keep
Julia single-threaded. This measures the current production CPU kernels rather
than thread scaling.

## Fock Hadamard interpretation

For the current real nearest-neighbour exchange convention, the complete
square Fock matrix can be written

```math
V_F = -U\,A_{\mathrm{nn}}\odot\operatorname{Re}(\rho),
```

where `odot` is the matrix-element Hadamard product. A local TT/MPO Hadamard
kernel could therefore replace both directional TCI reconstructions and their
translation-MPO assembly. Its uncompressed rank is bounded by the product of
the adjacency and density ranks, so it must be paired with controlled
compression. This benchmark first compares the already implemented
binary-carry Fock path against TCI; a Hadamard prototype should use those two
independent implementations as scientific references.
