using MPO_MeanField
using LinearAlgebra
using Printf

function purification_benchmark_parameters(
    ; itensors_tol=1e-12,
    itensors_maxdim=32,
    potential=(-2.0, -0.5, 0.7, 2.0),
)
    return Parameters1D(
        L=2,
        t=0.0,
        U=0.0,
        W=x -> potential[Int(x) + 1],
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=itensors_tol,
        itensors_maxdim=itensors_maxdim,
        density=0.5,
        purification_steps=35,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=5,
    )
end

function dense_matrix_for_benchmark(mpo, sys)
    N = 2^sys.params.L
    return [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states)
            for i in 1:N, j in 1:N]
end

function exact_projector(H, Ne)
    decomposition = eigen(Hermitian((H + H') / 2))
    vectors = decomposition.vectors[:, 1:Ne]
    return vectors * vectors'
end

function purification_benchmark_case(method, params, bounds)
    sys = System(params)
    rho0 = construct_rho_0(sys, params, bounds...; method=method)
    result = perform_purification(
        rho0, params;
        method=method,
        verbose=0,
        spectral_bounds=bounds,
        gc_policy=:automatic,
    )
    H = dense_matrix_for_benchmark(sys.H0, sys)
    rho = dense_matrix_for_benchmark(result.rho, sys)
    Ne = round(Int, size(H, 1) * params.density)
    exact = exact_projector(H, Ne)
    return (
        result=result,
        projector_error=opnorm(rho - exact),
        commutator=opnorm(H * rho - rho * H),
        band_energy=real(tr(H * rho)),
    )
end

function median(values)
    ordered = sort(values)
    return ordered[cld(length(ordered), 2)]
end

"""
    benchmark_purification_methods(; repetitions=3)

Compare SP2 and PM on a four-orbital non-interacting, gapped reference. The
two spectral intervals test affine scaling sensitivity; the two MPO settings
record tolerance/max-dimension sensitivity. This is a CPU microbenchmark, not
a production-size performance estimate.
"""
function benchmark_purification_methods(; repetitions=3)
    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    cases = (
        ("tight", purification_benchmark_parameters(), (-2.5, 2.5)),
        ("wide", purification_benchmark_parameters(itensors_tol=1e-8, itensors_maxdim=16), (-7.0, 7.0)),
        ("gap", purification_benchmark_parameters(potential=(-2.0, -0.01, 0.01, 2.0)), (-2.5, 2.5)),
    )
    for (_, params, bounds) in cases, method in (:sp2, :palser_manolopoulos)
        purification_benchmark_case(method, params, bounds) # compile/warm-up
    end

    println("case method  status              median_s  alloc_B  gc_s  live_heap_B rss_increase_B process_peak_rss_B  squares cubes  maxχ meanχ  trace_err idem  proj_err comm")
    for (label, params, bounds) in cases, method in (:sp2, :palser_manolopoulos)
        times = Float64[]
        allocations = Int[]
        gc_times = Float64[]
        latest = nothing
        rss_before = Sys.maxrss()
        for _ in 1:repetitions
            GC.gc(true)
            timing = @timed purification_benchmark_case(method, params, bounds)
            latest = timing.value
            push!(times, timing.time)
            push!(allocations, timing.bytes)
            push!(gc_times, timing.gctime)
        end
        result = latest.result
        rss_after = Sys.maxrss()
        @printf("%-4s %-20s %-18s %.6f %8d %.6f %11d %14d %20d %7d %5d %5d %.2f %.2e %.2e %.2e %.2e\n",
            label, String(method), String(result.termination_reason), median(times), median(allocations), median(gc_times),
            Base.gc_live_bytes(), max(0, rss_after - rss_before), rss_after,
            result.work.squares, result.work.cubes, result.work.max_bond_dimension,
            result.work.mean_bond_dimension, result.trace_error, result.idempotency_residual,
            latest.projector_error, latest.commutator,
        )
        labels = method == :sp2 ? ("rho", "rho2", "-") : ("rho", "rho2", "rho3")
        println("  bond dimensions by iteration ($(join(labels, ", "))): ", result.work.bond_dimensions)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_purification_methods()
end
