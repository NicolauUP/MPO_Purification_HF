# Atomistic TBLG circular-flake prototype

This directory is deliberately separate from the production square/rectangular
MPO and KPM workflows. It tests whether a finite atomistic twisted-bilayer
graph remains compressible after a quantics site ordering.

The current stage generates geometry and an on-demand hopping oracle. It does
not yet construct a TCI/MPO and it does not perform Hartree--Fock.

## Hopping convention

There is one carbon p_z orbital per site and no onsite term:

\[
H=\sum_{i,j}t_{ij}c_i^\dagger c_j,\qquad t_{ij}=t_{ji}\in\mathbb R.
\]

Intralayer hopping is restricted to the three honeycomb nearest neighbours and
uses intralayer_hopping=-2.70 in the prototype units. Interlayer hopping is
retained only for r_ij <= 6 and is

\[
t_{ij}=V_{pp\pi}(r_{ij})(1-\cos^2\gamma)
      +V_{pp\sigma}(r_{ij})\cos^2\gamma,
\qquad
\cos\gamma=\frac{z_i-z_j}{r_{ij}},
\]

\[
V_{pp\pi}(r)=V^0_{pp\pi}e^{-(r-a_0)/\delta},
\qquad
V_{pp\sigma}(r)=V^0_{pp\sigma}e^{-(r-d_0)/\delta}.
\]

The default prototype values are a0=1.42, d0=3.35,
Vpp_pi_0=-2.70, Vpp_sigma_0=0.48, and decay_length=0.319*a0. They are
recorded explicitly because they are exploratory values and should not silently
inherit the square-model parameters.

## First validation command

From the repository root:

~~~
julia --startup-file=no --project=. \
  benchmark/tblg_qtt/validate_tblg_circular_flake.jl 5.0 20.0
~~~

The output reports the physical site count, whether it is exactly 2^L,
layer/sublattice counts, padding needed for the next power of two, Hermiticity,
and the number of hopping entries.

For the centered disk at 5 degrees, a radius of 20 gives 962 physical sites,
not a power of two. This is expected. A shifted centre can be used to find a
discrete radius with exactly 1024 sites:

~~~
julia --startup-file=no -e '
include("benchmark/tblg_qtt/tblg_geometry.jl");
using .TBLGCircularFlake;
p=TBLGParameters(angle_deg=5.0);
r=find_exact_radius(1024,p; radius_max=40.0, center=(0.1,0.0));
println("radius=", r, " count=", length(circular_flake(r,p; center=(0.1,0.0))))
'
~~~

Run the deterministic unit checks with:

~~~
julia --startup-file=no benchmark/tblg_qtt/test_tblg_geometry.jl
~~~

## Binary-address validation

The next prototype step is a deterministic spatial binary address. The report
script finds the exact radius, applies a Morton ordering, verifies that the
hopping operator is unchanged under the resulting permutation, and compares
the hopping bandwidth with the existing coordinate sort:

~~~
julia --startup-file=no benchmark/tblg_qtt/report_tblg_binary_order.jl \
  1024 5.0 0.1 0.0 benchmark/tblg_qtt/results/tblg_binary_address_1024.csv
~~~

For the 5-degree, 1024-site test, the exact radius is
`20.818048992441863`. The permutation error is zero and Hermiticity is
preserved. In this first test the Morton bandwidth is 527, compared with 170
for the current coordinate sort, so Morton ordering is **not** yet a better
ordering for MPO compression. It remains optional; the production geometry
ordering has not been changed.

## Direct QTT-MPO construction prototype

`construct_tblg_hopping_mpo.jl` is the first representation test. It takes an
exact \(2^L\)-site sparse hopping matrix, writes it as an ITensor with \(L\)
output and \(L\) input binary indices, then uses the ITensorMPS SVD
factorization to construct an MPO. It compares coordinate order with a
deterministic reverse Cuthill--McKee (RCM) graph order, records the bond
dimensions, and validates 4096 nonzero plus 4096 random matrix elements.

~~~
julia --startup-file=no --project=. \\
  benchmark/tblg_qtt/construct_tblg_hopping_mpo.jl \\
  1024 5.0 1e-12 1024 benchmark/tblg_qtt/results/tblg_hopping_mpo_1024
~~~

The construction deliberately materializes a \(2^{2L}\)-entry tensor, so it
is restricted to 1024 sites and is only a **small-system MPO-compressibility
test**. It is not the scalable route: that must query the hopping oracle
through TCI/cross interpolation without materializing the matrix. In the
256-site smoke test, both coordinate and RCM orderings reproduced sampled
hopping elements to roughly \(3.3\times10^{-14}\); RCM reduced the peak bond
from 142 to 140. This is encouraging but does not establish the asymptotic
bond dimension.
