# New run checklist

Copy this checklist into the README of a new run directory before writing a
submission script.

- [ ] State the physical question and the expected qualitative result.
- [ ] Fix all model and numerical parameters, including spectral bounds.
- [ ] Choose `:sp2`, `:palser_manolopoulos`, or `:mcweeny_mu` and justify it.
- [ ] Define the smallest validation calculation.
- [ ] State the external cluster result root and a unique run identifier.
- [ ] Save Git revision, Julia version, and Slurm resource request with output.
- [ ] Define which final diagnostics make the run acceptable or inconclusive.
