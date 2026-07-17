using MPO_MeanField
using ITensors, ITensorMPS
using LinearAlgebra
using Printf

function tensorial_hartree_benchmark_system(L)
    params = Parameters1D(
        L=L,
        t=-0.7,
        U=0.3,
        W=x -> 0.17 * cospi(Int(x) / (2^L - 1)) + 0.05 * sinpi(3Int(x) / (2^L - 1)),
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=64,
        density=0.5,
        purification_steps=35,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=5,
    )
    sys = System(params)
    N = 2^L
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        x -> 0.45 + 0.12 * cospi(Int(x) / (N - 1)),
        Float64,
        sys.sites,
        params.tci_tol,
    )
    bond = MPO_MeanField.diagonal_mpo_from_function(
        x -> Int(x) < N - 1 ? 0.04 * sinpi(3Int(x) / (N - 1)) : 0.0,
        Float64,
        sys.sites,
        params.tci_tol,
    )
    T_R, T_L = sys.translations
    rho_right = apply(bond, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_left = apply(T_L, ITensors.dag(bond); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.ρ = +(diagonal, rho_right, rho_left;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return sys
end

function _benchmark_dense_matrix(mpo, sys)
    N = 2^sys.params.L
    return [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states)
            for i in 1:N, j in 1:N]
end

"""
    benchmark_tensorial_hartree(; L=3, repetitions=1)

Compare the cached-TCI and direct tensorial Hartree extractors on one small,
nonuniform Hermitian 1D density MPO. Returned allocation counts are cumulative
allocations from `@timed`, not peak memory. This is a local prototype benchmark,
not a production performance claim.
"""
function benchmark_tensorial_hartree(; L=3, repetitions=1)
    L > 0 || throw(ArgumentError("L must be positive"))
    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    sys = tensorial_hartree_benchmark_system(L)

    # Warm compilation and establish field equivalence before timing.
    tci = extract_hartree_mpo_1d(sys)
    tensorial = extract_hartree_mpo_tensorial_1d(sys)
    error_norm = opnorm(_benchmark_dense_matrix(tci, sys) - _benchmark_dense_matrix(tensorial, sys))

    tci_times = Float64[]
    tensorial_times = Float64[]
    tci_bytes = Int[]
    tensorial_bytes = Int[]
    tci_bond_dimensions = Int[]
    tensorial_bond_dimensions = Int[]
    for _ in 1:repetitions
        GC.gc(true)
        timed_tci = @timed extract_hartree_mpo_1d(sys)
        push!(tci_times, timed_tci.time)
        push!(tci_bytes, timed_tci.bytes)
        push!(tci_bond_dimensions, maxlinkdim(timed_tci.value))

        GC.gc(true)
        timed_tensorial = @timed extract_hartree_mpo_tensorial_1d(sys)
        push!(tensorial_times, timed_tensorial.time)
        push!(tensorial_bytes, timed_tensorial.bytes)
        push!(tensorial_bond_dimensions, maxlinkdim(timed_tensorial.value))
    end
    median(values) = sort(values)[cld(length(values), 2)]
    result = (
        L=L,
        field_error_opnorm=error_norm,
        tci=(time_s=median(tci_times), bytes=median(tci_bytes), max_bond_dimension=median(tci_bond_dimensions)),
        tensorial=(time_s=median(tensorial_times), bytes=median(tensorial_bytes), max_bond_dimension=median(tensorial_bond_dimensions)),
    )
    @printf("L=%d | ||VH_tci-VH_tensorial||₂=%.3e\n", result.L, result.field_error_opnorm)
    @printf("TCI:       median_s=%.6f bytes=%d χ=%d\n", result.tci.time_s, result.tci.bytes, result.tci.max_bond_dimension)
    @printf("tensorial: median_s=%.6f bytes=%d χ=%d\n", result.tensorial.time_s, result.tensorial.bytes, result.tensorial.max_bond_dimension)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_tensorial_hartree()
end
