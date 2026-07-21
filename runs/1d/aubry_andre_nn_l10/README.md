# Aubry--André nearest-neighbour HF, `L=10`

Three independent half-filled open-chain SCF calculations for comparison with
external ED results. The physical chain has `N=2^10=1024` sites and
`N_e=512` spinless particles. There is no onsite potential or seed.

The hopping is

\[
t(x)=-1-V_2\cos\left(2\pi\tau(x+\tfrac12)\right),
\qquad \tau=\sqrt2-5/6\approx0.580880229,
\]

with the package's existing nearest-neighbour real-exchange HF convention.

| Array task | `V2` | `U` | Conservative spectral bounds |
| --- | ---: | ---: | --- |
| 1 | 0.0 | 1.0 | `(-6.5, 6.5)` |
| 2 | 0.5 | 1.0 | `(-7.5, 7.5)` |
| 3 | 2.0 | 0.2 | `(-7.3, 7.3)` |

The bounds use the open-chain hopping row-sum bound plus `2|U|` each for
Hartree and real-exchange Fock contributions, with `0.5` additional padding.
They are deliberately conservative rather than ED-derived.

Submit from `MPO_Purification_HF`:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch runs/1d/aubry_andre_nn_l10/submit.slurm
```

Each array task writes an immutable external result directory containing the
exact `input.jl`, selection/provenance metadata, full SCF history, final
`site_density.csv`, nearest-neighbour `bond_order.csv`, and final scalar
observables. Compare ED and MPO results using the same open boundary, hopping
phase convention, half filling, and real-exchange approximation.
