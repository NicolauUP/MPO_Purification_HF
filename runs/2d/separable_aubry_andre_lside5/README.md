# Separable 2D Aubry--Andre HF, side level 5

This campaign is the square-lattice counterpart of the validated 1D
quasiperiodic reference cases. `side_level=5` means a `32 x 32` open square:
`ParametersSquare(L=10)`, `N=1024`, and `N_e=512` at half filling.

The hopping is separable by direction:

\[
t_x(x,y)=-1-V_2\cos[2\pi\tau(x+1/2)],\qquad
t_y(x,y)=-1-V_2\cos[2\pi\tau(y+1/2)],
\]

where \(\tau=\sqrt2-5/6\). There is no external onsite potential. The
temporary seed `S(x,y)=0.1(-1)^{x+y}` appears only in the first SCF
Hamiltonian and selects a checkerboard branch; it is replaced by the
self-consistent Hartree field.

| Array task | `V2` | `U` | Conservative SP2 bounds |
| --- | ---: | ---: | --- |
| 1 | 0.0 | 1.0 | `(-12.5, 12.5)` |
| 2 | 0.5 | 1.0 | `(-14.5, 14.5)` |
| 3 | 2.0 | 0.2 | `(-14.1, 14.1)` |

Submit the MPO calculation from the package root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch runs/2d/separable_aubry_andre_lside5/submit.slurm
```

Then submit the independent dense one-particle HF reference:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_dense_hf.slurm
```

The dense driver directly reconstructs the open-square Hamiltonian and
diagonalizes it at every SCF iteration; it does not form its reference from an
MPO. It writes under the separate campaign name
`separable_aubry_andre_lside5_seed0p1_dense_hf`, so both calculation paths
retain their densities, bond orders, SCF histories, energy components, and
provenance independently.

## P1.2: fixed-Hamiltonian SP2 cap diagnostic

Before retrying the MPO SCF calculation, run the initial-SP2 diagnostic:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_2_fixed_sp2.slurm
```

This submits one job per input Hamiltonian. Each job runs the same first SCF
Hamiltonian, `H0 + S`, at `maxdim = 128, 256, 384, 512`, without updating any
Hartree or Fock field. Every run writes `iterations.csv` with the trace,
candidate polynomial traces, branch, rank, cap-contact flag, idempotency, and
per-step resources. `summary.toml` records the production SP2 stopping result.

The result root is
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside5_seed0p1_fixed_sp2/`.
Interpret this experiment before changing a field extractor: if SP2 fails in
this isolated problem after reaching the cap, the limitation is purification
compression rather than SCF feedback.
