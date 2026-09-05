# QTT charge analysis

This directory validates compressed charge-density Fourier and correlation-
dimension calculations independently of the SCF solver.

`validate_qtt_fft_d2_ed.jl` builds a small open square Hamiltonian, obtains its
half-filled occupied projector by direct diagonalization, compresses that exact
projector as an MPO, and extracts its density diagonal as a QTT/MPS. It then
compares:

- every real-space density value against the exact projector diagonal;
- every normalized two-dimensional Fourier coefficient against a direct DFT;
- Parseval's identity;
- every dyadic correlation sum `Z2(s)` in real and Fourier space;
- the fitted `D2 = tau(2)` values.

Run the default `16 x 16` validation with:

```bash
julia --startup-file=no --project=. \
  benchmark/qtt_charge_analysis/validate_qtt_fft_d2_ed.jl \
  /tmp/qtt_fft_d2_ed
```

The optional second argument is `log2(side)`. Output is written atomically and
contains `metadata.toml`, `summary.toml`, and `scales.csv`.

The analyzed positive measure is

```math
\mu_{xy}=\frac{|n_{xy}-\bar n|^2}{\sum_{xy}|n_{xy}-\bar n|^2}.
```

Thus `D2` describes the spatial or momentum-space participation of the charge
modulation. It is not an orbital inverse participation ratio of the full
density matrix.

## Raw HDF5 export for plotting

`export_charge_hdf5.jl` reads a saved `converged_state.h5` (or
`final_fixed_sp2_state.h5`) and exports the diagonal charge density without
restarting SCF or constructing the dense density matrix. The exporter accepts
both public and legacy campaign files:

```bash
julia --startup-file=no --project=. \
  benchmark/qtt_charge_analysis/export_charge_hdf5.jl \
  CAMPAIGN.jl TASK RESULT_DIRECTORY OUTPUT.h5 auto 16777216 64
```

The final three arguments are optional: `mode` (`auto`, `dense`, or `sampled`),
the maximum number of sites allowed for a dense raster, and the sampling stride.
In dense mode the file contains a `density[side,side]` `Float64` dataset in the
public `(x,y)` ordering. In sampled mode it contains `sample_x`, `sample_y`,
and `sample_density`, evaluated directly from the QTT on a regular grid. Root
HDF5 attributes record the lattice size, QTT level count, source checkpoint,
SCF state kind, density definition, and extraction policy.

Use dense mode for publication plots when the raster is manageable (for example
`2048^2`), and sampled mode for very large lattices. The exporter never
overwrites an existing HDF5 file.

For cluster post-processing, set all paths explicitly and submit:

```bash
export CAMPAIGN=/gpfs/projects/epor78/MPO_HartreeFock_Purification/MPO_Purification_HF/runs/2d/c8_quasicrystal/campaign_double_c8_lside11_supermoire_mpo_sp2.jl
export TASK_INDEX=1
export RESULT_DIRECTORY=/gpfs/projects/epor78/MPO_HF_results/DOUBLE_C8_RESULT
export OUTPUT_H5=/gpfs/projects/epor78/MPO_HF_results/charge_analysis/double_c8.h5
sbatch --export=ALL,CAMPAIGN,TASK_INDEX,RESULT_DIRECTORY,OUTPUT_H5,MODE=auto \
  benchmark/qtt_charge_analysis/submit_charge_hdf5_export.slurm
```
