# Scientific conventions

## Model

The maintained model is real, spinless, open-boundary, nearest-neighbour tight-binding Hartree--Fock. The density matrix is

\[
\rho_{ij}=\langle c_j^\dagger c_i\rangle.
\]

Each undirected nearest-neighbour bond \(\langle ij\rangle\) is counted once:

\[
E[\rho]=\operatorname{Tr}(H_0\rho)+U\sum_{\langle ij\rangle}
\left[n_i n_j-\left(\operatorname{Re}\rho_{ij}\right)^2\right],
\qquad n_i=\rho_{ii}.
\]

`H0` contains hopping and external potential. There is no periodic wrapping, on-site Hubbard interaction, spin, or complex-exchange functional in the maintained MPO path. Fock intentionally uses `real(rho[i,j])`.

## Geometry and QTT encoding

- A chain with `L` binary bits has \(N=2^L\) open sites.
- A legacy equal-square MPO model with even `L` has side \(2^{L/2}\).
- Public `SquareModel(size=(Nx,Ny))` describes open power-of-two rectangles; public KPM supports them, but MPO/dense rectangles remain unvalidated.
- Julia site `i` is zero-based binary integer `i-1`.
- Equal-square bits interleave in integer-bit order \((y_0,x_0,y_1,x_1,\ldots)\); MPO tensors run most- to least-significant.

For separable Aubry--André hopping, campaign sources define \(t_x(x)=-1-V_2\cos[2\pi\tau(x+1/2)]\) and its y analogue. The campaign file is the executable specification of \(\tau\), seed, and spectral bounds.

## Numerical controls and diagnostics

`QTTSettings.tci_tol` controls function-to-QTT interpolation; `cutoff` and `maxdim` control MPO compression. Purification tolerance and maximum steps are inner-solver controls; SCF tolerance, mixing, Pulay safeguards, and two-cycle detection are outer-solver controls. The MPO spectral interval must enclose every SCF effective Hamiltonian and should be analytically justified.

For an explicit MPO projector, the zero-temperature diagnostic is a relative idempotency residual of the form

\[
\|\rho^2-\rho\|/\|\rho\|,
\]

distinct from trace error and outer SCF field residuals. KPM does not form \(\rho\) explicitly: it records a trace idempotency defect estimated from \(\operatorname{Tr}(P^2)-\operatorname{Tr}(P)\), normalized by occupied trace, and requires an independent higher-resolution fixed-field audit.

At zero temperature verify Hermiticity, filling, near-idempotency, stationarity, seed stability, and dense/reference agreement where possible.
