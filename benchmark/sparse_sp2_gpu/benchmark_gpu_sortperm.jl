"""Benchmark the exact CUDA.jl sortperm-and-gather route proposed for SpAMM.

Usage:
  julia --project=. benchmark_gpu_sortperm.jl OUTPUT [SIZES]

SIZES is a comma-separated list, defaulting to 1000000,3000000,10000000.
The timed operation includes sortperm(keys) and gathering one UInt64 payload,
which represents packed (left,right) block indices.
"""

using CUDA
using Random
using Statistics
using TOML

function main(arguments)
    length(arguments) in 1:2 || error("usage: benchmark_gpu_sortperm.jl OUTPUT [SIZES]")
    output = abspath(arguments[1])
    sizes = length(arguments) == 2 ? parse.(Int, split(arguments[2], ',')) :
            [1_000_000, 3_000_000, 10_000_000]
    CUDA.functional() || error("CUDA is not functional on this node")
    isdir(output) && error("refusing to overwrite existing output directory: $output")
    mkpath(output)

    Random.seed!(510578)
    open(joinpath(output, "sortperm_benchmark.csv"), "w") do io
    println(io, "records,minimum_time_s,median_time_s,maximum_time_s,sorted_correctly")
    for records in sizes
        host_keys = rand(UInt64, records)
        host_payload = rand(UInt64, records)
        keys = CuArray(host_keys)
        payload = CuArray(host_payload)
        key_work = similar(keys)
        permutation = CUDA.zeros(Int32, records)

        # Compile and populate the CUDA memory pool before timing.
        copyto!(key_work, keys)
        sortperm!(permutation, key_work; initialized=false)
        warm_payload = payload[permutation]
        CUDA.synchronize()
        warm_payload = nothing
        GC.gc(true)
        CUDA.reclaim()

        times = Float64[]
        final_sorted_keys = nothing
        for _ in 1:5
            CUDA.synchronize()
            elapsed = @elapsed begin
                copyto!(key_work, keys)
                sortperm!(permutation, key_work; initialized=false)
                sorted_keys = key_work
                sorted_payload = payload[permutation]
                CUDA.synchronize()
            end
            push!(times, elapsed)
            final_sorted_keys = sorted_keys
            sorted_payload = nothing
            GC.gc(true)
            CUDA.reclaim()
        end
        sorted_correctly = issorted(Array(final_sorted_keys))
        println(io, join((records, minimum(times), median(times), maximum(times),
                          sorted_correctly), ','))
        flush(io)
    end
    end

    open(joinpath(output, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
        "diagnostic" => "cuda_sortperm_preallocated_scratch_and_payload_gather",
            "device_name" => CUDA.name(CUDA.device()),
            "record_sizes" => sizes,
            "key_type" => "UInt64",
            "payload_type" => "UInt64",
        "repetitions" => 5,
        "input_keys_restored_before_each_sort" => true,
        "permutation_preallocated" => true,
        ))
    end
end

main(ARGS)
