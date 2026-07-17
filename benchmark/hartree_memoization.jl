using MPO_MeanField
using ITensors, ITensorMPS
using LinearAlgebra
using Printf

mutable struct UncachedHartreeEvaluator
    sys
    density_reads::Int
end

function (evaluator::UncachedHartreeEvaluator)(coordinate::Real)
    i = Int(coordinate) + 1
    value = 0.0
    for j in (i - 1, i + 1)
        if 1 <= j <= 2^evaluator.sys.params.L
            evaluator.density_reads += 1
            value += evaluator.sys.params.U * real(MatrixChecker(
                evaluator.sys.ρ, evaluator.sys.sites, j, j,
                evaluator.sys.bra_states, evaluator.sys.ket_states,
            ))
        end
    end
    return value
end

function memoization_benchmark_system()
    params = Parameters1D(
        L=3,
        t=-0.7,
        U=0.3,
        W=nothing,
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=32,
        density=0.5,
        purification_steps=20,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=5,
    )
    sys = System(params)
    sys.ρ = Identity_MPO(sys.sites) * 0.5
    return sys
end

function dense_matrix_for_hartree_benchmark(mpo, sys)
    N = 2^sys.params.L
    return [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states)
            for i in 1:N, j in 1:N]
end

function benchmark_hartree_memoization(; repetitions=3)
    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    sys = memoization_benchmark_system()
    params = sys.params

    # Compile both paths before timing. The closure counts TCI sample calls;
    # the cache size is the number of distinct diagonal matrix elements read.
    cached = MPO_MeanField.HartreeEvaluate1D(sys, Dict{Int,Float64}())
    cached_calls = Ref(0)
    cached_mpo = MPO_MeanField.diagonal_mpo_from_function(
        x -> (cached_calls[] += 1; cached(x)), Float64, sys.sites, params.tci_tol,
    )
    uncached = UncachedHartreeEvaluator(sys, 0)
    uncached_mpo = MPO_MeanField.diagonal_mpo_from_function(
        x -> uncached(x), Float64, sys.sites, params.tci_tol,
    )
    @assert opnorm(
        dense_matrix_for_hartree_benchmark(cached_mpo, sys) -
        dense_matrix_for_hartree_benchmark(uncached_mpo, sys),
    ) < 1e-10

    cached_times = Float64[]
    uncached_times = Float64[]
    cached_calls_per_run = Int[]
    cached_reads = Int[]
    uncached_reads = Int[]
    for _ in 1:repetitions
        cached = MPO_MeanField.HartreeEvaluate1D(sys, Dict{Int,Float64}())
        cached_calls = Ref(0)
        GC.gc(true)
        timed_cached = @timed MPO_MeanField.diagonal_mpo_from_function(
            x -> (cached_calls[] += 1; cached(x)), Float64, sys.sites, params.tci_tol,
        )
        push!(cached_times, timed_cached.time)
        push!(cached_calls_per_run, cached_calls[])
        push!(cached_reads, length(cached.density_cache))

        uncached = UncachedHartreeEvaluator(sys, 0)
        GC.gc(true)
        timed_uncached = @timed MPO_MeanField.diagonal_mpo_from_function(
            x -> uncached(x), Float64, sys.sites, params.tci_tol,
        )
        push!(uncached_times, timed_uncached.time)
        push!(uncached_reads, uncached.density_reads)
    end
    median(values) = sort(values)[cld(length(values), 2)]
    @printf("cached:   median_s=%.6f TCI_calls=%d unique_density_reads=%d\n",
        median(cached_times), median(cached_calls_per_run), median(cached_reads))
    @printf("uncached: median_s=%.6f density_reads=%d\n", median(uncached_times), median(uncached_reads))
    println("Both paths produced identical Hartree MPOs; cached evaluator calls are sampled by TCI.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_hartree_memoization()
end
