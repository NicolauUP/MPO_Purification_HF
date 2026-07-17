using MPO_MeanField
using Printf

function cleanup_benchmark_parameters()
    return Parameters1D(
        L=2,
        t=-0.7,
        U=0.3,
        W=x -> (-0.2, 0.1, -0.05, 0.25)[Int(x) + 1],
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=32,
        density=0.5,
        purification_steps=35,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=20,
    )
end

function run_cleanup_benchmark_case(policy; period=10, threshold_bytes=1 << 30)
    sys = System(cleanup_benchmark_parameters())
    run_scf!(sys, -5.0, 5.0;
        purification_method=:sp2,
        gc_policy=policy,
        gc_period=period,
        gc_threshold_bytes=threshold_bytes,
    ) || error("benchmark SCF did not converge for policy $policy")
    return nothing
end

function median(values)
    ordered = sort(values)
    return ordered[cld(length(ordered), 2)]
end

"""
    benchmark_cleanup_policies(; repetitions=3, gpu_peak_memory_bytes=nothing)

Run comparable CPU cleanup-policy measurements. If `gpu_peak_memory_bytes` is
provided, it must return the device allocator's peak bytes observed so far;
the benchmark records the maximum value returned after each repetition. This
keeps GPU instrumentation opt-in and avoids loading CUDA on CPU-only runs.
"""
function benchmark_cleanup_policies(; repetitions=3, gpu_peak_memory_bytes=nothing)
    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))

    # Compile the numerical path before collecting timings.
    run_cleanup_benchmark_case(:automatic)
    threshold = Base.gc_live_bytes() + 32 * 1024^2
    cases = (
        (:automatic, 10, 1 << 30),
        (:forced, 10, 1 << 30),
        (:periodic, 2, 1 << 30),
        (:threshold, 10, threshold),
    )
    println("policy       median_s   median_alloc_B  median_gc_s  live_heap_B  rss_increase_B  process_peak_rss_B  gpu_peak_B")
    for (policy, period, threshold_bytes) in cases
        times = Float64[]
        allocations = Int[]
        gc_times = Float64[]
        gpu_peaks = Int[]
        rss_before = Sys.maxrss()
        for _ in 1:repetitions
            GC.gc(true)
            timing = @timed run_cleanup_benchmark_case(
                policy; period=period, threshold_bytes=threshold_bytes,
            )
            push!(times, timing.time)
            push!(allocations, timing.bytes)
            push!(gc_times, timing.gctime)
            !isnothing(gpu_peak_memory_bytes) && push!(gpu_peaks, gpu_peak_memory_bytes())
        end
        rss_after = Sys.maxrss()
        gpu_peak = isnothing(gpu_peak_memory_bytes) ? "unavailable" : string(maximum(gpu_peaks))
        @printf("%-12s %.6f   %14d  %.6f    %11d  %14d  %18d  %s\n",
            String(policy), median(times), median(allocations), median(gc_times),
            Base.gc_live_bytes(), max(0, rss_after - rss_before), rss_after, gpu_peak,
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_cleanup_policies()
end
