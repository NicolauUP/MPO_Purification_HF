# Validation matrix

This is an evidence ledger, not a claim that every regime is validated.

| Path / regime | Independent comparison | Status | Remaining requirement |
|---|---|---|---|
| Dense HF, small open 1D/equal squares | direct eigensystem and invariants | reference | retain deterministic regressions |
| MPO-SP2, gapped 32x32 controls | dense projector and local observables | supported at recorded settings | report cap/cutoff and SP2 residual |
| MPO-SP2, low-gap 64x64 separable | fixed-H compression and cap ladders | difficult | trajectory-level convergence study |
| KPM-SCF, 32x32 separable | dense HF | validated high-resolution reference | retain moment/probe/seed checks |
| KPM-SCF, 64x64--256x256 | nested M/R and fixed-field audits | resolution-dependent | select M/R only after audit/seed agreement |
| KPM-SCF, 1024x512 and 2048x2048 | higher-resolution fixed-field audit | scalability evidence | report audit differences; no dense ED |
| Rectangular MPO/dense HF | none | unsupported | complete rectangular MPO validation R1--R6 |
| Non-separable quasiperiodic hopping | none | future model | dense/KPM checks and M/R/seed convergence |

## Required reporting

MPO results require spectral bounds, cutoff, cap, SP2 steps, trace, idempotency, Hermiticity, cap hits, SCF residuals, stationarity, and seed dependence. KPM results require `M`, deterministic probe design and `R`, chemical-potential/trace solve, trace-idempotency defect, SCF residuals, and independent audit differences. Compare nested moments at fixed probes and nested probes at fixed moments separately.

Use errors appropriate to the observable: a global order parameter can agree while local-density errors remain visible in a localized regime.
