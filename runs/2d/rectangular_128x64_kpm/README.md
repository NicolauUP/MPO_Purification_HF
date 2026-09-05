# All-local KPM versus ED at N=8192

An exactly square power-of-two lattice cannot contain 8192 sites, so this
diagnostic uses an open `128 x 64` rectangle. It retains the separable
`V2=0.5` hopping, `tau=sqrt(2)-5/6`, checkerboard seed `0.1`, half filling,
the analytically safe interval `[-6.35,6.35]`, and 800 Jackson--Chebyshev
moments.

The driver diagonalizes the dense `8192 x 8192` Hamiltonian once. It then runs
Hadamard probe counts 256, 512, 1024, 2048, and 4096; counts 512, 1024, and
2048 are repeated with another randomized sign seed. Every one of the 8192
densities and all 16,192 open-boundary nearest-neighbour bonds are compared
with the occupied eigenspace.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/rectangular_128x64_kpm/submit_all_local_ed_cuda.slurm
```

The single job reuses its ED reference across the full probe ladder. Results
are stored below
`$MPO_RESULTS_ROOT/rectangular_128x64_kpm_all_local_ed/`.

After that reference finishes, compare a physically balanced Hadamard
hierarchy without repeating ED:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/rectangular_128x64_kpm/submit_coordinate_hadamard_cuda.slurm
```

This retains contiguous x-fastest row storage but assigns each row a
Hadamard code with interleaved low `(x,y)` bits. It compares
`R=256,512,1024,2048,4096,8192`, with second seeds at 512, 1024, and 2048.
The full `R=N=8192` result isolates the 800-moment polynomial error.
