# AGENTS.md — Repository Operating Guidelines

This repository (`MPO_Purification_HF`) implements high-performance Hartree–Fock and purification algorithms in Julia, targeting large-scale 1D and 2D lattice models using Matrix Product Operators (MPO), the Kernel Polynomial Method (KPM), custom GPU Block-SpAMM, and cuSPARSE CSR baselines.

Every agent working in this repository must strictly observe the architectural boundaries, CUDA safety invariants, and cluster orchestration rules documented below.

---

## 1. Directory Structure & File Placement

- `src/`: Core library code only.
  - `src/solvers/spamm/block_sp2.jl`: Module `BlockSp2Engine` containing GPU Block-SpAMM data structures (`CuBlockCSR`, `HostBlockCSR`) and CUDA kernels.
  - All solver engines must be proper Julia modules. Never use `SPAMM_LIBRARY_ONLY=true` or script-level include hacks.
- `ClusterRuns/`: The **sole** location for Slurm job submission scripts.
  - All cluster jobs must use the 5 canonical, parameterized templates:
    1. `block_spamm.slurm` — Block-SpAMM scaling, diagnostics, and SCF.
    2. `cusparse_csr.slurm` — cuSPARSE CSR scaling and SCF.
    3. `mpohf_cuda.slurm` — Unified dispatcher for `bin/mpohf.jl`.
    4. `kpm.slurm` — 2D square and rectangular KPM SCF and probing.
    5. `square_diagnostic.slurm` — Fixed-H, McWeeny, and projector compression diagnostics.
  - **Rule:** Never create one-off `.slurm` launchers in `runs/`, `benchmark/`, or subdirectories. Parameterize existing `ClusterRuns/` templates instead.
- `runs/`: Physical campaigns and scientific drivers.
  - Subdirectories (`runs/1d/`, `runs/2d/`) contain campaign specifications (`campaign.jl`) and model parameters.
  - `runs/common/`: Reusable drivers, validation tools, and shared runtimes (`slurm_cuda_runtime.sh`).
  - **Rule:** Never commit incident-specific debug scripts, snapshot replays, or transient CUDA reproducer scripts here.
- `benchmark/`: Standalone microbenchmarks and profiling harnesses.
- `logs/`: Slurm stdout and stderr dumps (`logs/%x_%j.out`). Always ignored by Git.

---

## 2. CUDA & Numerical Invariants

- **`CUDA.sortperm!` Safety Invariant:**
  `CUDA.sortperm!` mutates key vectors in-place on the device. **Always make an explicit copy of candidate keys before calling `CUDA.sortperm!`** to prevent silent memory corruption in SpAMM radix sorting.
- **Inner Loop Synchronization:**
  Never introduce blocking host-device synchronizations (e.g., `Array(cu_tensor)`, scalar indexing, or synchronous error logging) inside SP2 iterations or inner SCF loops.
- **Physical Sanity Checks:**
  - Hermiticity: $H = H^\dagger$, $\rho = \rho^\dagger$.
  - Particle number conservation: $|\text{Tr}(\rho) - N_e| / N_e \le 10^{-6}$.
  - Idempotency: $\|\rho^2 - \rho\| \le 10^{-8}$ at zero temperature.
  - No unphysical spectral intruder states with fractional occupation in the band gap.

---

## 3. Slurm & Execution Rules

- **Pre-Execution Directory Creation:**
  Slurm opens output log streams before executing the bash script body. Ensure `logs/` exists prior to dispatch, and always include `mkdir -p logs "$OUTPUT_DIR"` inside launchers.
- **Julia Output Path Collision:**
  Julia campaign scripts reject pre-existing target directories. Job templates must route runs into subdirectories per task (e.g., `"$OUTPUT_DIR/$TASK"`).
- **CLI Argument Inspection:**
  Before writing or adjusting any runner, inspect the target Julia script's argument parser (`ARGS`) directly. Never guess positional argument orders.
- **Script Validation:**
  Every shell/Slurm script created or modified must pass `bash -n` validation before concluding a task.

---

## 4. Deletion & Modification Guardrails

- **Campaign Files are Sacred:**
  Never delete `campaign.jl` or files defining physical models (e.g., C8 quasicrystal, Aubry–André, strong CDW) without explicit human confirmation. They encode irreplaceable scientific provenance.
- **Refactoring Strategy:**
  When replacing legacy scripts with canonical templates, first verify the new template with `bash -n` and check argument parity before pruning legacy files.