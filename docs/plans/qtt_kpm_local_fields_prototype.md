# Prototype: QTT-KPM local-field Hartree--Fock

## Motivation

The MPO-SP2 path repeatedly compresses a full density-projector MPO. In a
low-gap quasiperiodic regime, its nonlinear squaring trajectory can saturate
the MPO cap before the desired projector is reached. This prototype explores
a distinct route: approximate only the local density and nearest-neighbour
bond fields needed by Hartree--Fock, using Chebyshev propagation of MPS probe
states and QTT representations of the resulting local fields.

This is a feasibility plan, not an SCF solver or a claim of equivalence to
current KPM or MPO paths.

## Initial fixed-H result (2026-08-25)

`benchmark/qtt_kpm_mps_probe/diagnose_qtt_kpm_mps_probe.jl` now implements
the first diagnostic. It propagates *rank-one, deterministic* Walsh--Hadamard
MPS probes by the MPO--MPS recurrence and separately evaluates the identical
dense-vector recurrence in QTT order.

For the 32x32 separable `V2=0.5, U=1` seeded fixed Hamiltonian, at
`M=512`, one probe, `maxdim=256`, and cutoff `1e-14`, the propagated MPS
never exceeded chi=32 and its final state differed from the vector recurrence
by `6.1e-7` in relative 2-norm. This is evidence that the *single-probe
recurrence* is compressible in this control case. It is not yet evidence that
the route is economical: it ran on CPU, materialised every amplitude for the
comparison, and a single probe has large estimator error.

The production KPM Hadamard probes additionally multiply every site by a
random sign. That diagonal random sign field is generally not rank-one in
QTT, so this first diagnostic intentionally uses unscrambled Walsh rows. A
probe-variance comparison against the production nested randomized estimator
is required before changing the solver.

## Required local fields

For the real nearest-neighbour functional, SCF requires only

\[
n_i=\rho_{ii},\qquad b_x(i)=\rho_{i,i+\hat{x}},\qquad
b_y(i)=\rho_{i,i+\hat{y}},
\]

not the complete density matrix. These determine the Hartree potential and
the real nearest-neighbour Fock field.

## Deterministic selected-site oracle

For binary site \(i\), form the product MPS \(\lvert i\rangle\), then
propagate with scaled Hamiltonian \(\widetilde H\):

\[
\lvert v_0^{(i)}\rangle=\lvert i\rangle,
\quad \lvert v_1^{(i)}\rangle=\widetilde H\lvert i\rangle,
\quad
\lvert v_{n+1}^{(i)}\rangle=2\widetilde H\lvert v_n^{(i)}\rangle-
\lvert v_{n-1}^{(i)}\rangle.
\]

After accumulating Jackson--Chebyshev coefficients,
\(\lvert y^{(i)}\rangle\approx P_M(H)\lvert i\rangle\). The oracle returns

\[
n_i=\langle i\vert y^{(i)}\rangle,
\qquad \rho_{i,i+\hat\delta}=\langle i\vert P_M(H)\vert i+\hat\delta\rangle.
\]

TCI/QTT cross interpolation can query selected points of \(n,b_x,b_y\).
Each query still performs a full MPO--MPS Chebyshev recurrence; TCI reduces
queries but cannot share their Krylov histories.

## Hadamard-probe QTT estimator

A Walsh--Hadamard row is rank-one in the binary QTT representation:

\[
r_a(i)=(-1)^{a\cdot i}.
\]

For each probe, propagate \(\lvert y_\alpha\rangle=P_M(H)\lvert r_\alpha\rangle\).
Its local estimators are

\[
n_i^{(\alpha)}=r_{\alpha,i}y_{\alpha,i},\qquad
\rho_{ij}^{(\alpha)}=y_{\alpha,i}r_{\alpha,j}.
\]

Average diagonal and nearest-neighbour bands directly as QTT/MPS fields. A
complete orthogonal Hadamard basis is exact in exact arithmetic; a smaller
randomized/nested subset is an approximation requiring independent probe
convergence checks.

## Feasibility sequence

1. Fixed-H single-probe 32x32 test at \(M=256\): rank trajectory and local
   contribution versus dense/vector KPM.
2. Fixed-H probe ladder \(R=8,16,32,\ldots\): compress averaged QTT density
   and bonds; compare ranks and local errors with dense observables.
3. QTT-compress exact field arrays independently, separating target-field
   rank from propagation-rank limitations.
4. Test selected-site TCI only if probes have high rank but fields remain
   QTT-compressible.
5. Validate a complete local-field QTT-KPM SCF against dense HF before any
   64x64 calculation.

## Required checks

Record Chebyshev order, interval, chemical-potential trace solve, MPS/QTT
cutoff and cap, rank trajectory, trace/idempotency diagnostics, local
density/bond errors, and wall time. Do not infer accurate propagation merely
from a compressible final fitted field.
