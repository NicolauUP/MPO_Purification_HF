raw"""
    nearest_neighbor_hf_energy_1d(sys)

Return the components of the zero-temperature, open-chain nearest-neighbour
Hartree--Fock energy implemented by the current 1D field extractors:

```math
E = \mathrm{Tr}(H_0\rho)
  + U\sum_{i=1}^{N-1}\left(n_i n_{i+1} - [\operatorname{Re}\rho_{i,i+1}]^2\right).
```

The returned named tuple has `kinetic`, `hartree`, `fock`, `interaction`, and
`total` fields. `H0` includes both hopping and the external potential `W`.
Each physical bond is counted once, so there is no factor of two in the
explicit interaction sum. There is no on-site interaction term and hence no
on-site self-interaction contribution.

The present Fock implementation uses `real(rho[i,i+1])`; this routine follows
that implemented real-density functional. It is therefore not an energy
functional for a complex exchange implementation.
"""
function nearest_neighbor_hf_energy_1d(sys::System)
    sys.params isa Parameters1D || throw(ArgumentError(
        "nearest-neighbor HF energy is currently implemented only for Parameters1D",
    ))

    params = sys.params
    N = 2^params.L
    density = Vector{Float64}(undef, N)
    bond_order = Vector{Float64}(undef, N - 1)
    for i in 1:N
        density[i] = real(MatrixChecker(
            sys.ρ, sys.sites, i, i, sys.bra_states, sys.ket_states,
        ))
        if i < N
            bond_order[i] = real(MatrixChecker(
                sys.ρ, sys.sites, i, i + 1, sys.bra_states, sys.ket_states,
            ))
        end
    end

    hartree = params.U * sum(density[i] * density[i + 1] for i in 1:(N - 1))
    fock = -params.U * sum(abs2(bond_order[i]) for i in 1:(N - 1))
    kinetic = real(tr(sys.H0 * sys.ρ))
    interaction = hartree + fock
    return (
        kinetic=kinetic,
        hartree=hartree,
        fock=fock,
        interaction=interaction,
        total=kinetic + interaction,
    )
end
