# Current status and next work

## Objective

The conference target is the two-dimensional quasiperiodic nearest-neighbour Hartree--Fock model, initially with separable hopping modulation. The project uses QTT/MPO purification for tensor-network research and GPU KPM for large-scale physical exploration and independent references.

## Established evidence

- Small equal-square dense references validate Hartree/Fock conventions, purification diagnostics, and the public MPO/dense wrappers.
- A 32x32 KPM-SCF versus dense-HF reference at \(V_2=0.5,U=1\) reached density maximum error about \(1.7\times10^{-6}\), density RMS about \(4.8\times10^{-7}\), and bond maximum error about \(7.0\times10^{-7}\).
- Full CUDA KPM-SCF has been exercised at 1024x512 and 2048x2048 sites with deterministic coordinate-interleaved Hadamard probes and fixed-field audits.
- Low-gap 64x64 separable MPO-SP2 is hard: it reaches the cap and can stagnate. Exact-projector compression is useful, but does not validate a truncated SP2 trajectory.

## Active MPO risk

Repeated `apply(A,A; cutoff,maxdim)` drives SP2. Zip-up multiplication makes local contractions/factorizations followed by canonical truncation. In the low-gap case, accumulated trajectory compression changes the purification map before the desired projector is reached. Current H100 profiling of physical and synthetic MPO squaring is a performance diagnostic, not physics evidence.

An alternative feasibility path is documented in
[`plans/qtt_kpm_local_fields_prototype.md`](plans/qtt_kpm_local_fields_prototype.md):
Chebyshev MPO--MPS propagation of rank-one Hadamard probes followed by QTT
local density/bond fields. It remains pending until the physical SP2 profile
is interpreted.

## Recommended next tasks

1. Finish the `apply(A,A)` profile and identify factorization versus data movement as the measured bottleneck.
2. Perform a controlled 64x64 MPO SP2 cap/cutoff/trajectory study, retaining best measured iterates and fixed-H references.
3. Complete separable KPM moment/probe/seed convergence before mapping it.
4. Introduce non-separable hopping only after small dense/KPM validation and an explicit error budget.

## Do not infer

KPM audit agreement does not validate MPO purification; SCF convergence does not prove KPM resolution convergence; exact-projector compressibility does not prove a truncated SP2 trajectory; and public rectangle specification does not mean MPO/dense rectangular HF is implemented.
