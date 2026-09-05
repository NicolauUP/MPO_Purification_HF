# Local GLMakie charge-density plots

This is intentionally a separate plotting environment. It does not modify the
solver `Project.toml` and it never reruns SCF.

Instantiate it on a machine with package access:

```bash
julia --project=analysis/charge_plots -e 'using Pkg; Pkg.instantiate()'
```

After exporting a dense HDF5 raster on the cluster, render locally:

```bash
julia --project=analysis/charge_plots \
  analysis/charge_plots/render_charge_report.jl \
  density_export.h5 figures/c8_double 10 128
```

The renderer writes `charge_report.png`, `charge_report.pdf`, `fft_peaks.csv`,
and `metadata.txt`. The full FFT uses the entire raster with a Hann window;
edge trimming masks the outer margin but retains the full lattice in the plot.
Sampled HDF5 files produce a sampled-density figure and do not pretend to
provide a full FFT.

For the double-C8 model, the publication layout also shows the centered
hopping modulation (with bond-center offsets) alongside the density and FFT:

```bash
julia --project=analysis/charge_plots \
  analysis/charge_plots/render_double_c8_publication.jl \
  density_export.h5 figures/double_c8_publication 10
```

This writes a high-resolution GLMakie PNG and `fft_peaks.csv`.

For presentation figures, use the raster `image!` renderer. It keeps the two
bare bond directions separate, writes full and 128×128 density views, and
computes all-bond correlations between the initial hopping and each endpoint
density separately (the endpoint average is retained only as a reference):

```bash
julia --project=analysis/charge_plots render_double_c8_powerpoint.jl \
  density_export.h5 figures/double_c8_powerpoint 10
```

The output includes raw and checkerboard-demodulated density images, separate
`δtₓ`/`δtᵧ` 128×128 crops, and `density_vs_initial_hopping_scatter.png`.

The scatter diagnostic reports endpoint densities separately.  An endpoint
average can cancel a staggered checkerboard response, so the figure also shows
`(-1)^(x+y)(n-1/2)` at each endpoint.  The direct response of a hopping field is
the bond order `2 Re ρᵢⱼ` (and hence the Fock field), which requires exporting
the converged bond order in addition to the density for a separate diagnostic.
