"""GPU front-end diagnostic for block SpAMM and GPU symbolic scheduling.

Implemented stages:
  1. GPU-resident block-CSR storage and Frobenius block norms.
  2. GPU generation of all compatible block-product candidates (i,k,j).
  3. GPU screening with ||A_ik||_F ||B_kj||_F >= tau.
  4. GPU compaction of retained candidates.
  5. GPU leaf multiplication and collision-safe output-block accumulation.

Every GPU result is checked against an independent host calculation, including
the conservative Frobenius error bound for the approximate product.

Usage:
  julia --project=. diagnose_block_spamm_gpu_frontend.jl OUTPUT [N] [BLOCK_SIZE] [V2] [SP2_STEPS]

When V2 is supplied, the input is an actual intermediate from fixed-H SP2 for
the hopping Aubry--Andre chain, with a staggered diagonal seed of +0.5.
"""

using CUDA
using LinearAlgebra
using SparseArrays
using Statistics
using TOML

const SPAMM_LIBRARY_ONLY = lowercase(get(ENV, "SPAMM_LIBRARY_ONLY", "false")) in
                           ("1", "true", "yes")
(!SPAMM_LIBRARY_ONLY && !(length(ARGS) in 1:5)) && error(
    "usage: diagnose_block_spamm_gpu_frontend.jl OUTPUT [N] [BLOCK_SIZE] [V2] [SP2_STEPS]",
)
output = SPAMM_LIBRARY_ONLY ? "" : abspath(ARGS[1])
N = SPAMM_LIBRARY_ONLY ? 16 : (length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 16384)
block_size = SPAMM_LIBRARY_ONLY ? 16 : (length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 16)
V2 = SPAMM_LIBRARY_ONLY ? nothing : (length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : nothing)
snapshot_steps = SPAMM_LIBRARY_ONLY ? 0 : (length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 12)
N % block_size == 0 || error("N must be divisible by BLOCK_SIZE")
block_size <= 32 || error("this diagnostic currently requires BLOCK_SIZE <= 32")
CUDA.functional() || error("CUDA is not functional on this node")
if !SPAMM_LIBRARY_ONLY
    isdir(output) && error("refusing to overwrite existing output directory: $output")
    mkpath(output)
end

struct HostBlockCSR
    nblockrows::Int
    block_size::Int
    rowptr::Vector{Int32}
    rowind::Vector{Int32}
    colind::Vector{Int32}
    values::Array{Float64,3}
end

struct CuBlockCSR
    nblockrows::Int
    block_size::Int
    rowptr::CuVector{Int32}
    rowind::CuVector{Int32}
    colind::CuVector{Int32}
    values::CuArray{Float64,3}
    norms::CuVector{Float64}
end

"A symmetric matrix with density-matrix-like exponential spatial decay."
function density_like_matrix(N::Int)
    rows = Int[]
    columns = Int[]
    values = Float64[]
    radius = 48
    sizehint!(rows, (2radius + 1) * N)
    sizehint!(columns, (2radius + 1) * N)
    sizehint!(values, (2radius + 1) * N)
    for i in 1:N
        for j in max(1, i - radius):min(N, i + radius)
            distance = abs(i - j)
            value = i == j ? 0.5 :
                0.12 * cos(0.37 * (i + j)) * cospi(distance / 2) * exp(-distance / 9)
            abs(value) >= 1e-10 || continue
            push!(rows, i); push!(columns, j); push!(values, value)
        end
    end
    A = sparse(rows, columns, values, N, N)
    return sparse(0.5 * (A + A'))
end

function csr_trace_kernel!(out, rowptr, colval, nzval, N)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if row <= N
        value = 0.0
        for position in Int(rowptr[row]):(Int(rowptr[row + 1]) - 1)
            if colval[position] == row
                value = nzval[position]
                break
            end
        end
        CUDA.@atomic out[1] += value
    end
    return
end

function csr_trace(A)
    result = CUDA.zeros(Float64, 1)
    threads = 256
    @cuda threads=threads blocks=cld(size(A, 1), threads) csr_trace_kernel!(
        result, A.rowPtr, A.colVal, A.nzVal, size(A, 1),
    )
    return Array(result)[1]
end

function host_csr(A, N::Int)
    rowptr = Array(A.rowPtr)
    columns = Int.(Array(A.colVal))
    values = Array(A.nzVal)
    rows = Vector{Int}(undef, length(values))
    for row in 1:N
        rows[Int(rowptr[row]):(Int(rowptr[row + 1]) - 1)] .= row
    end
    return sparse(rows, columns, values, N, N)
end

function fixed_hamiltonian(N::Int, V2::Float64, hopping_radius::Int,
                           hopping_decay_length::Float64)
    tau_aa = sqrt(2.0) - 5.0 / 6.0
    diagonal = [iseven(site) ? 0.5 : -0.5 for site in 1:N]
    diagonals = Pair{Int,Vector{Float64}}[0 => diagonal]
    maximum_row_sum = 0.0
    for distance in 1:hopping_radius
        envelope = exp(-(distance - 1) / hopping_decay_length)
        hopping = [
            -envelope * (1 + V2 * cos(2π * tau_aa * (site + distance / 2)))
            for site in 1:(N - distance)
        ]
        push!(diagonals, -distance => hopping, distance => hopping)
        maximum_row_sum += 2maximum(abs, hopping)
    end
    H = spdiagm(diagonals...)
    radius = maximum(abs, diagonal) + maximum_row_sum + 0.25
    return H, radius
end

function fixed_h_sp2_snapshot(N::Int, V2::Float64, steps::Int,
                             hopping_radius::Int=1,
                             hopping_decay_length::Float64=2.0)
    H, radius = fixed_hamiltonian(N, V2, hopping_radius,
                                  hopping_decay_length)
    initial = sparse((radius * I - H) / (2radius))
    P = CUDA.CUSPARSE.CuSparseMatrixCSR(initial)
    particles = N ÷ 2
    cutoff = 1e-8
    for iteration in 1:steps
        branch = csr_trace(P) >= particles ? :square : :hole
        product = P * P
        product = CUDA.CUSPARSE.prune(product, cutoff, 'O')
        P = branch == :square ? product : CUDA.CUSPARSE.geam(2.0, P, -1.0, product, 'O')
        P = CUDA.CUSPARSE.prune(P, cutoff, 'O')
        CUDA.synchronize()
        println("snapshot SP2 $iteration/$steps branch=$branch nnz_per_row=$(nnz(P) / N)")
    end
    return host_csr(P, N)
end

function host_block_csr(A::SparseMatrixCSC{Float64,Int}, block_size::Int)
    nblockrows = cld(size(A, 1), block_size)
    dictionary = Dict{Tuple{Int,Int},Matrix{Float64}}()
    rows, columns, entries = findnz(A)
    for (row, column, entry) in zip(rows, columns, entries)
        block_row, local_row = divrem(row - 1, block_size)
        block_column, local_column = divrem(column - 1, block_size)
        block = get!(() -> zeros(block_size, block_size), dictionary,
                     (block_row + 1, block_column + 1))
        block[local_row + 1, local_column + 1] = entry
    end

    keys_sorted = sort!(collect(keys(dictionary)))
    row_counts = zeros(Int32, nblockrows)
    for (block_row, _) in keys_sorted
        row_counts[block_row] += 1
    end
    rowptr = Vector{Int32}(undef, nblockrows + 1)
    rowptr[1] = 1
    for row in 1:nblockrows
        rowptr[row + 1] = rowptr[row] + row_counts[row]
    end
    rowind = Int32[key[1] for key in keys_sorted]
    colind = Int32[key[2] for key in keys_sorted]
    values = zeros(Float64, block_size, block_size, length(keys_sorted))
    for (position, key) in enumerate(keys_sorted)
        values[:, :, position] .= dictionary[key]
    end
    return HostBlockCSR(nblockrows, block_size, rowptr, rowind, colind, values)
end

function block_norm_kernel!(norms, values, block_size, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        total = 0.0
        for column in 1:block_size, row in 1:block_size
            total += abs2(values[row, column, block])
        end
        norms[block] = sqrt(total)
    end
    return
end

function cu_block_csr(host::HostBlockCSR)
    rowptr = CuArray(host.rowptr)
    rowind = CuArray(host.rowind)
    colind = CuArray(host.colind)
    values = CuArray(host.values)
    norms = CUDA.zeros(Float64, length(host.colind))
    threads = 256
    @cuda threads=threads blocks=cld(length(host.colind), threads) block_norm_kernel!(
        norms, values, host.block_size, length(host.colind),
    )
    CUDA.synchronize()
    return CuBlockCSR(host.nblockrows, host.block_size, rowptr, rowind, colind, values, norms)
end

"One output segment per left block; offsets are exclusive and one-based."
function candidate_offsets(host::HostBlockCSR)
    counts = zeros(Int64, length(host.colind))
    for left in eachindex(host.colind)
        k = Int(host.colind[left])
        counts[left] = Int(host.rowptr[k + 1] - host.rowptr[k])
    end
    offsets = Vector{Int64}(undef, length(counts) + 1)
    offsets[1] = 1
    for position in eachindex(counts)
        offsets[position + 1] = offsets[position] + counts[position]
    end
    return offsets
end

function candidate_kernel!(left_index, right_index, output_row, output_column,
                           scores, keep, rowptr, rowind, colind, norms, offsets, tau,
                           number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        i = rowind[left]
        k = Int(colind[left])
        first_right = Int(rowptr[k])
        last_right = Int(rowptr[k + 1]) - 1
        destination = Int(offsets[left])
        for right in first_right:last_right
            score = norms[left] * norms[right]
            left_index[destination] = left
            right_index[destination] = right
            output_row[destination] = i
            output_column[destination] = colind[right]
            scores[destination] = score
            keep[destination] = score >= tau
            destination += 1
        end
    end
    return
end

function gpu_candidates(A::CuBlockCSR, offsets::Vector{Int64}, tau::Float64)
    number_of_candidates = offsets[end] - 1
    left_index = CUDA.zeros(Int32, number_of_candidates)
    right_index = similar(left_index)
    output_row = similar(left_index)
    output_column = similar(left_index)
    scores = CUDA.zeros(Float64, number_of_candidates)
    keep = CUDA.zeros(Bool, number_of_candidates)
    device_offsets = CuArray(offsets)
    threads = 256
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=threads blocks=cld(length(A.colind), threads) candidate_kernel!(
            left_index, right_index, output_row, output_column, scores, keep,
            A.rowptr, A.rowind, A.colind, A.norms, device_offsets, tau, length(A.colind),
        )
        CUDA.synchronize()
    end
    return (; left_index, right_index, output_row, output_column, scores, keep,
            elapsed, number_of_candidates)
end

function host_candidate_reference(A::HostBlockCSR, norms::Vector{Float64}, tau::Float64)
    candidates = Tuple{Int32,Int32,Int32,Int32,Float64,Bool}[]
    for i in 1:A.nblockrows
        for left in Int(A.rowptr[i]):(Int(A.rowptr[i + 1]) - 1)
            k = Int(A.colind[left])
            for right in Int(A.rowptr[k]):(Int(A.rowptr[k + 1]) - 1)
                score = norms[left] * norms[right]
                push!(candidates, (left, right, i, A.colind[right], score, score >= tau))
            end
        end
    end
    return candidates
end

function output_block_indices(output_rows::Vector{Int32}, output_columns::Vector{Int32})
    keys = sort!(unique!(collect(zip(output_rows, output_columns))))
    lookup = Dict(key => Int32(position) for (position, key) in enumerate(keys))
    indices = Int32[lookup[(row, column)] for (row, column) in zip(output_rows, output_columns)]
    return keys, indices
end

function compact_kernel!(compact_left, compact_right, compact_output, counter,
                         left, right, output_index, keep, number_of_candidates)
    candidate = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if candidate <= number_of_candidates && keep[candidate]
        destination = CUDA.atomic_add!(pointer(counter), Int32(1)) + Int32(1)
        compact_left[destination] = left[candidate]
        compact_right[destination] = right[candidate]
        compact_output[destination] = output_index[candidate]
    end
    return
end

function compact_candidates(gpu, output_index::Vector{Int32})
    compact_left = similar(gpu.left_index)
    compact_right = similar(gpu.right_index)
    compact_output = similar(gpu.left_index)
    counter = CUDA.zeros(Int32, 1)
    device_output_index = CuArray(output_index)
    threads = 256
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=threads blocks=cld(gpu.number_of_candidates, threads) compact_kernel!(
            compact_left, compact_right, compact_output, counter, gpu.left_index,
            gpu.right_index, device_output_index, gpu.keep, gpu.number_of_candidates,
        )
        CUDA.synchronize()
    end
    retained = Int(Array(counter)[1])
    return (; compact_left, compact_right, compact_output, retained, elapsed)
end

function multiply_accumulate_kernel!(output_values, values, left_indices,
                                     right_indices, output_indices, retained,
                                     block_size)
    product = blockIdx().x
    element = threadIdx().x
    elements_per_block = block_size * block_size
    if product <= retained && element <= elements_per_block
        row = (element - 1) % block_size + 1
        column = (element - 1) ÷ block_size + 1
        left = left_indices[product]
        right = right_indices[product]
        output = output_indices[product]
        value = 0.0
        for inner in 1:block_size
            value += values[row, inner, left] * values[inner, column, right]
        end
        linear_index = row + (column - 1) * block_size +
                       (output - 1) * elements_per_block
        CUDA.atomic_add!(pointer(output_values, linear_index), value)
    end
    return
end

function multiply_accumulate(device::CuBlockCSR, compact, number_of_output_blocks::Int)
    output_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                               number_of_output_blocks)
    CUDA.synchronize()
    elapsed = @elapsed begin
        if compact.retained > 0
            @cuda threads=device.block_size^2 blocks=compact.retained multiply_accumulate_kernel!(
                output_values, device.values, compact.compact_left,
                compact.compact_right, compact.compact_output, compact.retained,
                device.block_size,
            )
        end
        CUDA.synchronize()
    end
    return output_values, elapsed
end

function blocks_to_sparse(keys, values::Array{Float64,3}, N::Int, block_size::Int)
    rows = Int[]
    columns = Int[]
    entries = Float64[]
    for (position, (block_row, block_column)) in enumerate(keys)
        row_offset = (Int(block_row) - 1) * block_size
        column_offset = (Int(block_column) - 1) * block_size
        for column in 1:block_size, row in 1:block_size
            value = values[row, column, position]
            iszero(value) && continue
            push!(rows, row_offset + row)
            push!(columns, column_offset + column)
            push!(entries, value)
        end
    end
    return sparse(rows, columns, entries, N, N)
end

function conservative_error_bound(reference)
    skipped_by_output = Dict{Tuple{Int32,Int32},Float64}()
    for candidate in reference
        candidate[6] && continue
        key = (candidate[3], candidate[4])
        skipped_by_output[key] = get(skipped_by_output, key, 0.0) + candidate[5]
    end
    return sqrt(sum(abs2, values(skipped_by_output)))
end

function prune_sparse!(A::SparseMatrixCSC{Float64,Int}, cutoff::Float64)
    dropzeros!(A)
    droptol!(A, cutoff)
    return A
end

"Fixed-H SP2 continuation with phase-resolved end-to-end timings."
function validated_sp2_continuation(P0::SparseMatrixCSC{Float64,Int}, block_size::Int,
                                    particles::Int, tau::Float64, steps::Int,
                                    output_cutoff::Float64, output::String;
                                    validate_reference::Bool=true,
                                    dense_projector::Union{Nothing,Matrix{Float64}}=nothing)
    approximate = copy(P0)
    reference = validate_reference ? copy(P0) : nothing
    open(joinpath(output, "full_sp2_history.csv"), "w") do io
        println(io, "iteration,branch,trace,relative_trace_error,idempotency_defect,input_nnz,input_blocks,candidates,screened_fraction,output_blocks,block_build_time_s,host_to_device_time_s,schedule_time_s,screening_time_s,multiply_time_s,device_to_host_time_s,sparse_rebuild_prune_time_s,trace_branch_time_s,total_iteration_time_s,reference_time_s,relative_reference_error,density_max_abs_error,gpu_free_memory_bytes")
        for iteration in 1:steps
            total_start = time_ns()
            trace_start = time_ns()
            trace_before = tr(approximate)
            branch = trace_before >= particles ? :square : :hole
            trace_branch_time = (time_ns() - trace_start) / 1e9

            block_start = time_ns()
            host = host_block_csr(approximate, block_size)
            block_build_time = (time_ns() - block_start) / 1e9

            transfer_start = time_ns()
            device = cu_block_csr(host)
            host_to_device_time = (time_ns() - transfer_start) / 1e9

            schedule_start = time_ns()
            schedule = grouped_output_schedule(host)
            schedule_time = (time_ns() - schedule_start) / 1e9

            screening_start = time_ns()
            scores = [norm(view(host.values, :, :, left)) * norm(view(host.values, :, :, right))
                      for (left, right) in zip(schedule.left, schedule.right)]
            screened_fraction = count(score -> score < tau, scores) / max(length(scores), 1)
            screening_time = (time_ns() - screening_start) / 1e9

            square_values, spamm_time = grouped_multiply(device, schedule, tau)

            device_to_host_start = time_ns()
            host_square_values = Array(square_values)
            device_to_host_time = (time_ns() - device_to_host_start) / 1e9

            rebuild_start = time_ns()
            square = blocks_to_sparse(schedule.keys, host_square_values,
                                      size(approximate, 1), block_size)
            approximate = branch == :square ? square : sparse(2approximate - square)
            prune_sparse!(approximate, output_cutoff)
            rebuild_prune_time = (time_ns() - rebuild_start) / 1e9

            reference_time = 0.0
            if validate_reference
                reference_start = time_ns()
                reference_square = reference * reference
                reference = branch == :square ? reference_square : sparse(2reference - reference_square)
                prune_sparse!(reference, output_cutoff)
                reference_time = (time_ns() - reference_start) / 1e9
            end

            trace_value = tr(approximate)
            idempotency = abs(sum(abs2, nonzeros(approximate)) - trace_value) / particles
            relative_error = validate_reference ?
                norm(approximate - reference) / max(norm(reference), eps()) : NaN
            density_error = validate_reference ?
                maximum(abs, diag(approximate) - diag(reference)) : NaN
            total_iteration_time = (time_ns() - total_start) / 1e9 - reference_time
            println(io, join((
                iteration, branch, trace_value, abs(trace_value - particles) / particles,
                idempotency, nnz(approximate), length(host.colind), length(schedule.left),
                screened_fraction, length(schedule.keys), block_build_time,
                host_to_device_time, schedule_time, screening_time, spamm_time,
                device_to_host_time, rebuild_prune_time, trace_branch_time,
                total_iteration_time, reference_time, relative_error, density_error,
                CUDA.free_memory(),
            ), ','))
            flush(io)
        end
    end

    trace_value = tr(approximate)
    idempotency = abs(sum(abs2, nonzeros(approximate)) - trace_value) / particles
    dense_relative_error = isnothing(dense_projector) ? NaN :
        norm(Matrix(approximate) - dense_projector) / norm(dense_projector)
    dense_density_error = isnothing(dense_projector) ? NaN :
        maximum(abs, diag(approximate) - diag(dense_projector))
    open(joinpath(output, "full_sp2_summary.toml"), "w") do io
        TOML.print(io, Dict(
            "matrix_dimension" => size(approximate, 1),
            "particles" => particles,
            "iterations" => steps,
            "spamm_tau" => tau,
            "output_cutoff" => output_cutoff,
            "reference_enabled" => validate_reference,
            "dense_validation_enabled" => !isnothing(dense_projector),
            "trace" => trace_value,
            "relative_trace_error" => abs(trace_value - particles) / particles,
            "trace_idempotency_defect" => idempotency,
            "final_nnz" => nnz(approximate),
            "final_nnz_per_row" => nnz(approximate) / size(approximate, 1),
            "dense_projector_relative_error" => dense_relative_error,
            "dense_density_max_abs_error" => dense_density_error,
            "timing_excludes_reference" => true,
            "gpu_resident_numeric_multiplication" => true,
            "host_structural_rebuild" => true,
        ))
    end
    return approximate, reference
end

"Symbolic candidate schedule grouped by output block C[i,j]."
function grouped_output_schedule(A)
    contributions = Dict{Tuple{Int32,Int32},Vector{Tuple{Int32,Int32}}}()
    for i in 1:A.nblockrows
        for left in Int(A.rowptr[i]):(Int(A.rowptr[i + 1]) - 1)
            k = Int(A.colind[left])
            for right in Int(A.rowptr[k]):(Int(A.rowptr[k + 1]) - 1)
                key = (Int32(i), A.colind[right])
                push!(get!(() -> Tuple{Int32,Int32}[], contributions, key),
                      (Int32(left), Int32(right)))
            end
        end
    end
    output_keys = sort!(collect(Base.keys(contributions)))
    rowptr = Vector{Int64}(undef, length(output_keys) + 1)
    rowptr[1] = 1
    left = Int32[]
    right = Int32[]
    for (position, key) in enumerate(output_keys)
        for (left_index, right_index) in contributions[key]
            push!(left, left_index)
            push!(right, right_index)
        end
        rowptr[position + 1] = length(left) + 1
    end
    return (; keys=output_keys, rowptr, left, right)
end

"Count screened-compatible right blocks for every active left block."
function screened_candidate_count_kernel!(counts, rowptr, colind, norms, tau,
                                          number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        k = Int(colind[left])
        count = Int32(0)
        for right in Int(rowptr[k]):(Int(rowptr[k + 1]) - 1)
            count += ifelse(norms[left] * norms[right] >= tau, Int32(1), Int32(0))
        end
        counts[left] = count
    end
    return
end

function screened_candidate_offset_kernel!(offsets, inclusive_counts,
                                            number_of_left_blocks)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index == 1 && (offsets[1] = Int64(1))
    if index <= number_of_left_blocks
        offsets[index + 1] = Int64(inclusive_counts[index]) + Int64(1)
    end
    return
end

"Write only candidates that passed the Frobenius-norm SpAMM screen."
function screened_candidate_write_kernel!(keys, left_indices, right_indices,
                                          rowptr, rowind, colind, norms, offsets,
                                          tau, number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        i = UInt64(rowind[left])
        k = Int(colind[left])
        destination = Int(offsets[left])
        for right in Int(rowptr[k]):(Int(rowptr[k + 1]) - 1)
            if norms[left] * norms[right] >= tau
                j = UInt64(colind[right])
                keys[destination] = (i << 32) | j
                left_indices[destination] = Int32(left)
                right_indices[destination] = Int32(right)
                destination += 1
            end
        end
    end
    return
end

"GPU count-scan-write candidate generation without rejected-candidate storage."
function gpu_prescreened_candidates(device::CuBlockCSR, tau::Float64)
    number_of_left_blocks = length(device.colind)
    counts = CUDA.zeros(Int32, number_of_left_blocks)
    threads = 256
    CUDA.synchronize()
    count_time = @elapsed begin
        @cuda threads=threads blocks=cld(number_of_left_blocks, threads) screened_candidate_count_kernel!(
            counts, device.rowptr, device.colind, device.norms, tau,
            number_of_left_blocks,
        )
        CUDA.synchronize()
    end

    scan_time = @elapsed begin
        inclusive_counts = CUDA.cumsum(counts)
        offsets = CUDA.zeros(Int64, number_of_left_blocks + 1)
        @cuda threads=threads blocks=cld(number_of_left_blocks, threads) screened_candidate_offset_kernel!(
            offsets, inclusive_counts, number_of_left_blocks,
        )
        CUDA.synchronize()
    end
    retained = Int(Array(view(offsets, number_of_left_blocks + 1:
                                      number_of_left_blocks + 1))[1] - 1)
    keys = CUDA.zeros(UInt64, retained)
    left = CUDA.zeros(Int32, retained)
    right = CUDA.zeros(Int32, retained)
    write_time = @elapsed begin
        if retained > 0
            @cuda threads=threads blocks=cld(number_of_left_blocks, threads) screened_candidate_write_kernel!(
                keys, left, right, device.rowptr, device.rowind,
                device.colind, device.norms, offsets, tau,
                number_of_left_blocks,
            )
        end
        CUDA.synchronize()
    end
    return (; keys, left, right, retained, count_time, scan_time, write_time)
end

mutable struct CandidateSortWorkspace
    capacity::Int
    keys
    permutation
    left
    right
end

CandidateSortWorkspace() = CandidateSortWorkspace(
    0, CUDA.zeros(UInt64, 0), CUDA.zeros(Int32, 0),
    CUDA.zeros(Int32, 0), CUDA.zeros(Int32, 0),
)

function ensure_sort_capacity!(workspace::CandidateSortWorkspace, required::Int)
    required <= workspace.capacity && return workspace
    capacity = max(required, max(1024, 2 * workspace.capacity))
    workspace.keys = CUDA.zeros(UInt64, capacity)
    workspace.permutation = CUDA.zeros(Int32, capacity)
    workspace.left = CUDA.zeros(Int32, capacity)
    workspace.right = CUDA.zeros(Int32, capacity)
    workspace.capacity = capacity
    return workspace
end

function gather_sorted_candidate_payload_kernel!(sorted_left, sorted_right,
                                                 left, right, permutation,
                                                 number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates
        source = permutation[index]
        sorted_left[index] = left[source]
        sorted_right[index] = right[source]
    end
    return
end

"Sort copied keys and gather payloads without mutating the candidate records."
function gpu_sort_candidates(candidates,
                             workspace::CandidateSortWorkspace=CandidateSortWorkspace())
    number_of_candidates = length(candidates.keys)
    CUDA.synchronize()
    elapsed = @elapsed begin
        permutation = CUDA.sortperm(candidates.keys)
        sorted_keys = candidates.keys[permutation]
        sorted_left = candidates.left[permutation]
        sorted_right = candidates.right[permutation]
        CUDA.synchronize()
    end
    return (; keys=sorted_keys, left=sorted_left, right=sorted_right, elapsed)
end

function group_boundary_kernel!(boundaries, keys, number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates
        boundaries[index] = index == 1 || keys[index] != keys[index - 1]
    end
    return
end

function group_scatter_kernel!(unique_keys, group_rowptr, keys, boundaries,
                               group_ids, number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates && boundaries[index] != 0
        group = Int(group_ids[index])
        unique_keys[group] = keys[index]
        group_rowptr[group] = Int64(index)
    end
    return
end

function group_finalize_kernel!(group_rowptr, number_of_groups,
                                number_of_candidates)
    if blockIdx().x == 1 && threadIdx().x == 1
        group_rowptr[number_of_groups + 1] = Int64(number_of_candidates + 1)
    end
    return
end

"GPU run-length encoding of sorted output keys."
function gpu_group_sorted_candidates(sorted)
    number_of_candidates = length(sorted.keys)
    number_of_candidates > 0 || error("cannot group an empty candidate list")
    threads = 256
    boundaries = CUDA.zeros(Int32, number_of_candidates)
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=threads blocks=cld(number_of_candidates, threads) group_boundary_kernel!(
            boundaries, sorted.keys, number_of_candidates,
        )
        group_ids = CUDA.cumsum(boundaries)
        CUDA.synchronize()
        number_of_groups = Int(Array(view(group_ids, number_of_candidates:
                                                    number_of_candidates))[1])
        1 <= number_of_groups <= number_of_candidates || error(
            "invalid GPU run-length group count: groups=$(number_of_groups), " *
            "candidates=$(number_of_candidates)",
        )
        unique_keys = CUDA.zeros(UInt64, number_of_groups)
        group_rowptr = CUDA.zeros(Int64, number_of_groups + 1)
        @cuda threads=threads blocks=cld(number_of_candidates, threads) group_scatter_kernel!(
            unique_keys, group_rowptr, sorted.keys, boundaries, group_ids,
            number_of_candidates,
        )
        @cuda threads=1 blocks=1 group_finalize_kernel!(
            group_rowptr, number_of_groups, number_of_candidates,
        )
        CUDA.synchronize()
    end
    return (; keys=unique_keys, rowptr=group_rowptr, left=sorted.left,
            right=sorted.right, elapsed)
end

decode_output_key(key::UInt64) = (Int32(key >> 32), Int32(key & 0xffffffff))

function grouped_multiply_device(device::CuBlockCSR, schedule)
    output_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                               length(schedule.keys))
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=device.block_size^2 blocks=length(schedule.keys) shmem=(
            2 * sizeof(Float64) * device.block_size^2
        ) grouped_multiply_kernel!(
            output_values, device.values, device.norms, schedule.rowptr,
            schedule.left, schedule.right, 0.0, device.block_size,
            length(schedule.keys),
        )
        CUDA.synchronize()
    end
    return output_values, elapsed
end

function pack_block_keys_kernel!(keys, rowind, colind, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        keys[block] = (UInt64(rowind[block]) << 32) | UInt64(colind[block])
    end
    return
end

function device_block_keys(device::CuBlockCSR)
    keys = CUDA.zeros(UInt64, length(device.colind))
    threads = 256
    @cuda threads=threads blocks=cld(length(keys), threads) pack_block_keys_kernel!(
        keys, device.rowind, device.colind, length(keys),
    )
    return keys
end

"Return sorted distinct device keys using boundary marking and a scan."
function gpu_unique_sorted_keys(keys)
    length(keys) > 0 || error("cannot uniquify an empty key list")
    sorted_keys = CUDA.sort(keys)
    boundaries = CUDA.zeros(Int32, length(sorted_keys))
    threads = 256
    @cuda threads=threads blocks=cld(length(sorted_keys), threads) group_boundary_kernel!(
        boundaries, sorted_keys, length(sorted_keys),
    )
    group_ids = CUDA.cumsum(boundaries)
    CUDA.synchronize()
    number_of_groups = Int(Array(view(group_ids, length(group_ids):length(group_ids)))[1])
    1 <= number_of_groups <= length(keys) || error(
        "invalid GPU unique-key count: groups=$(number_of_groups), " *
        "inputs=$(length(keys))",
    )
    unique_keys = CUDA.zeros(UInt64, number_of_groups)
    dummy_rowptr = CUDA.zeros(Int64, number_of_groups + 1)
    @cuda threads=threads blocks=cld(length(sorted_keys), threads) group_scatter_kernel!(
        unique_keys, dummy_rowptr, sorted_keys, boundaries, group_ids,
        length(sorted_keys),
    )
    CUDA.synchronize()
    return unique_keys
end

function exact_key_lookup_kernel!(source_to_union, source_keys, union_keys,
                                  number_of_source_keys, number_of_union_keys)
    source = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if source <= number_of_source_keys
        key = source_keys[source]
        low = 1
        high = number_of_union_keys
        found = 0
        while low <= high
            middle = (low + high) >>> 1
            candidate = union_keys[middle]
            if candidate < key
                low = middle + 1
            elseif candidate > key
                high = middle - 1
            else
                found = middle
                break
            end
        end
        source_to_union[source] = Int32(found)
    end
    return
end

function scatter_inverse_map_kernel!(union_to_source, source_to_union,
                                     number_of_source_keys)
    source = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if source <= number_of_source_keys
        destination = source_to_union[source]
        destination > 0 && (union_to_source[destination] = Int32(source))
    end
    return
end

function identity_map_kernel!(mapping, number_of_entries)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index <= number_of_entries && (mapping[index] = Int32(index))
    return
end

function map_sorted_keys_to_union(source_keys, union_keys)
    source_to_union = CUDA.zeros(Int32, length(source_keys))
    union_to_source = CUDA.zeros(Int32, length(union_keys))
    threads = 256
    @cuda threads=threads blocks=cld(length(source_keys), threads) exact_key_lookup_kernel!(
        source_to_union, source_keys, union_keys, length(source_keys),
        length(union_keys),
    )
    @cuda threads=threads blocks=cld(length(source_keys), threads) scatter_inverse_map_kernel!(
        union_to_source, source_to_union, length(source_keys),
    )
    return union_to_source
end

"Construct the selected SP2 output structure and source maps entirely on CUDA."
function gpu_branch_union(device::CuBlockCSR, square_schedule, branch::Symbol)
    CUDA.synchronize()
    elapsed = @elapsed begin
        input_keys = device_block_keys(device)
        if branch == :square
            output_keys = square_schedule.keys
            output_to_square = CUDA.zeros(Int32, length(output_keys))
            @cuda threads=256 blocks=cld(length(output_keys), 256) identity_map_kernel!(
                output_to_square, length(output_keys),
            )
            output_to_input = CUDA.zeros(Int32, length(output_keys))
        elseif branch == :hole
            output_keys = gpu_unique_sorted_keys(vcat(square_schedule.keys, input_keys))
            output_to_square = map_sorted_keys_to_union(square_schedule.keys, output_keys)
            output_to_input = map_sorted_keys_to_union(input_keys, output_keys)
        else
            error("unsupported SP2 branch: $branch")
        end
        CUDA.synchronize()
    end
    return (; keys=output_keys, square=output_to_square, input=output_to_input,
            elapsed)
end

function fused_branch_multiply_kernel!(output_values, input_values, square_rowptr,
                                       square_left, square_right, output_to_square,
                                       output_to_input, branch_sign, input_factor,
                                       output_cutoff, block_size,
                                       number_of_output_blocks)
    output = blockIdx().x
    element = threadIdx().x
    elements_per_block = block_size * block_size
    if output <= number_of_output_blocks && element <= elements_per_block
        row = (element - 1) % block_size + 1
        column = (element - 1) ÷ block_size + 1
        left_tile = CuDynamicSharedArray(Float64, (block_size, block_size))
        right_tile = CuDynamicSharedArray(
            Float64, (block_size, block_size), sizeof(Float64) * elements_per_block,
        )
        input_block = Int(output_to_input[output])
        value = input_block == 0 ? 0.0 :
                input_factor * input_values[row, column, input_block]
        square_group = Int(output_to_square[output])
        if square_group != 0
            first_contribution = Int(square_rowptr[square_group])
            last_contribution = Int(square_rowptr[square_group + 1]) - 1
            for contribution in first_contribution:last_contribution
                left = square_left[contribution]
                right = square_right[contribution]
                left_tile[row, column] = input_values[row, column, left]
                right_tile[row, column] = input_values[row, column, right]
                sync_threads()
                product = 0.0
                for inner in 1:block_size
                    product += left_tile[row, inner] * right_tile[inner, column]
                end
                value += branch_sign * product
                sync_threads()
            end
        end
        output_values[row, column, output] =
            abs(value) < output_cutoff ? 0.0 : value
    end
    return
end

function fused_branch_multiply(device::CuBlockCSR, square_schedule, branch_union,
                               branch::Symbol, output_cutoff::Float64)
    output_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                               length(branch_union.keys))
    branch_sign = branch == :square ? 1.0 : -1.0
    input_factor = branch == :square ? 0.0 : 2.0
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=device.block_size^2 blocks=length(branch_union.keys) shmem=(
            2 * sizeof(Float64) * device.block_size^2
        ) fused_branch_multiply_kernel!(
            output_values, device.values, square_schedule.rowptr,
            square_schedule.left, square_schedule.right, branch_union.square,
            branch_union.input, branch_sign, input_factor, output_cutoff,
            device.block_size, length(branch_union.keys),
        )
        CUDA.synchronize()
    end
    return output_values, elapsed
end

function mark_retained_blocks_kernel!(keep, norms, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    block <= number_of_blocks &&
        (keep[block] = ifelse(norms[block] > 0.0, Int32(1), Int32(0)))
    return
end

function compact_block_metadata_kernel!(keys_out, norms_out, keys_in, norms_in,
                                        keep, destinations, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks && keep[block] != 0
        destination = Int(destinations[block])
        keys_out[destination] = keys_in[block]
        norms_out[destination] = norms_in[block]
    end
    return
end

function compact_block_values_kernel!(values_out, values_in, keep, destinations,
                                      elements_per_block, number_of_blocks)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    total = elements_per_block * number_of_blocks
    if index <= total
        source_block = (index - 1) ÷ elements_per_block + 1
        if keep[source_block] != 0
            local_element = (index - 1) % elements_per_block + 1
            destination_block = Int(destinations[source_block])
            destination = local_element + (destination_block - 1) * elements_per_block
            values_out[destination] = values_in[index]
        end
    end
    return
end

function decode_keys_and_count_rows_kernel!(rowind, colind, row_counts, keys,
                                            number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        key = keys[block]
        row = Int32(key >> 32)
        column = Int32(key & 0xffffffff)
        rowind[block] = row
        colind[block] = column
        CUDA.@atomic row_counts[row] += Int32(1)
    end
    return
end

function int32_rowptr_kernel!(rowptr, inclusive_counts, number_of_rows)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    row == 1 && (rowptr[1] = Int32(1))
    row <= number_of_rows &&
        (rowptr[row + 1] = Int32(inclusive_counts[row]) + Int32(1))
    return
end

"Prune and rebuild the next CuBlockCSR without copying structural arrays to host."
function gpu_compact_output(device::CuBlockCSR, output_keys, output_values)
    threads = 256
    number_of_blocks = length(output_keys)
    output_norms = CUDA.zeros(Float64, number_of_blocks)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) block_norm_kernel!(
        output_norms, output_values, device.block_size, number_of_blocks,
    )
    keep = CUDA.zeros(Int32, number_of_blocks)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) mark_retained_blocks_kernel!(
        keep, output_norms, number_of_blocks,
    )
    destinations = CUDA.cumsum(keep)
    CUDA.synchronize()
    retained = Int(Array(view(destinations, number_of_blocks:number_of_blocks))[1])

    retained_keys = CUDA.zeros(UInt64, retained)
    retained_norms = CUDA.zeros(Float64, retained)
    retained_values = CUDA.zeros(Float64, device.block_size, device.block_size, retained)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) compact_block_metadata_kernel!(
        retained_keys, retained_norms, output_keys, output_norms, keep,
        destinations, number_of_blocks,
    )
    @cuda threads=threads blocks=cld(length(output_values), threads) compact_block_values_kernel!(
        retained_values, output_values, keep, destinations, device.block_size^2,
        number_of_blocks,
    )

    rowind = CUDA.zeros(Int32, retained)
    colind = CUDA.zeros(Int32, retained)
    row_counts = CUDA.zeros(Int32, device.nblockrows)
    @cuda threads=threads blocks=cld(retained, threads) decode_keys_and_count_rows_kernel!(
        rowind, colind, row_counts, retained_keys, retained,
    )
    inclusive_counts = CUDA.cumsum(row_counts)
    rowptr = CUDA.zeros(Int32, device.nblockrows + 1)
    @cuda threads=threads blocks=cld(device.nblockrows, threads) int32_rowptr_kernel!(
        rowptr, inclusive_counts, device.nblockrows,
    )
    CUDA.synchronize()
    return CuBlockCSR(device.nblockrows, device.block_size, rowptr, rowind, colind,
                      retained_values, retained_norms)
end

function grouped_multiply_kernel!(output_values, values, norms, group_rowptr,
                                  left_indices, right_indices, tau, block_size,
                                  number_of_output_blocks)
    output = blockIdx().x
    element = threadIdx().x
    elements_per_block = block_size * block_size
    if output <= number_of_output_blocks && element <= elements_per_block
        row = (element - 1) % block_size + 1
        column = (element - 1) ÷ block_size + 1
        left_tile = CuDynamicSharedArray(Float64, (block_size, block_size))
        right_tile = CuDynamicSharedArray(
            Float64, (block_size, block_size), sizeof(Float64) * elements_per_block,
        )
        value = 0.0
        first_contribution = Int(group_rowptr[output])
        last_contribution = Int(group_rowptr[output + 1]) - 1
        for contribution in first_contribution:last_contribution
            left = left_indices[contribution]
            right = right_indices[contribution]
            if norms[left] * norms[right] >= tau
                left_tile[row, column] = values[row, column, left]
                right_tile[row, column] = values[row, column, right]
                sync_threads()
                for inner in 1:block_size
                    value += left_tile[row, inner] * right_tile[inner, column]
                end
                sync_threads()
            end
        end
        output_values[row, column, output] = value
    end
    return
end

function grouped_multiply(device::CuBlockCSR, schedule, tau::Float64)
    device_rowptr = CuArray(schedule.rowptr)
    device_left = CuArray(schedule.left)
    device_right = CuArray(schedule.right)
    output_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                               length(schedule.keys))
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=device.block_size^2 blocks=length(schedule.keys) shmem=(
            2 * sizeof(Float64) * device.block_size^2
        ) grouped_multiply_kernel!(
            output_values, device.values, device.norms, device_rowptr, device_left,
            device_right, tau, device.block_size, length(schedule.keys),
        )
        CUDA.synchronize()
    end
    return output_values, elapsed
end

function scale_blocks_kernel!(values, factor, number_of_elements)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index <= number_of_elements && (values[index] *= factor)
    return
end

function add_mapped_blocks_kernel!(output_values, input_values, input_to_output,
                                   factor, elements_per_block, number_of_input_blocks)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    total = elements_per_block * number_of_input_blocks
    if index <= total
        input_block = (index - 1) ÷ elements_per_block + 1
        local_element = (index - 1) % elements_per_block + 1
        output_block = input_to_output[input_block]
        output_index = local_element + (output_block - 1) * elements_per_block
        output_values[output_index] += factor * input_values[index]
    end
    return
end

function prune_elements_kernel!(values, cutoff, number_of_elements)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_elements && abs(values[index]) < cutoff
        values[index] = 0.0
    end
    return
end

function gather_blocks_kernel!(destination, source, retained_indices,
                               elements_per_block, retained_blocks)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    total = elements_per_block * retained_blocks
    if index <= total
        destination_block = (index - 1) ÷ elements_per_block + 1
        local_element = (index - 1) % elements_per_block + 1
        source_block = retained_indices[destination_block]
        source_index = local_element + (source_block - 1) * elements_per_block
        destination[index] = source[source_index]
    end
    return
end

function block_metrics_kernel!(metrics, values, rowind, colind, block_size,
                               number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        squared_norm = 0.0
        trace_value = 0.0
        diagonal_block = rowind[block] == colind[block]
        for column in 1:block_size, row in 1:block_size
            value = values[row, column, block]
            squared_norm += value * value
            if diagonal_block && row == column
                trace_value += value
            end
        end
        CUDA.@atomic metrics[2] += squared_norm
        if diagonal_block
            CUDA.@atomic metrics[1] += trace_value
        end
    end
    return
end

function device_block_metrics(device::CuBlockCSR)
    metrics = CUDA.zeros(Float64, 2)
    number_of_blocks = length(device.colind)
    threads = 256
    @cuda threads=threads blocks=cld(number_of_blocks, threads) block_metrics_kernel!(
        metrics, device.values, device.rowind, device.colind, device.block_size,
        number_of_blocks,
    )
    CUDA.synchronize()
    result = Array(metrics)
    return (; trace=result[1], squared_norm=result[2])
end

function host_structure(device::CuBlockCSR)
    return (;
        nblockrows=device.nblockrows,
        block_size=device.block_size,
        rowptr=Array(device.rowptr),
        rowind=Array(device.rowind),
        colind=Array(device.colind),
    )
end

function structure_from_keys(keys, nblockrows::Int)
    row_counts = zeros(Int32, nblockrows)
    for (row, _) in keys
        row_counts[row] += 1
    end
    rowptr = Vector{Int32}(undef, nblockrows + 1)
    rowptr[1] = 1
    for row in 1:nblockrows
        rowptr[row + 1] = rowptr[row] + row_counts[row]
    end
    return (;
        rowptr,
        rowind=Int32[key[1] for key in keys],
        colind=Int32[key[2] for key in keys],
    )
end

function branch_prune_and_compact(device::CuBlockCSR, square_values,
                                  output_keys, branch::Symbol, cutoff::Float64)
    elements_per_block = device.block_size^2
    threads = 256
    lookup = Dict(key => Int32(position) for (position, key) in enumerate(output_keys))
    input_keys = collect(zip(Array(device.rowind), Array(device.colind)))
    all(haskey(lookup, key) for key in input_keys) || error(
        "P^2 symbolic structure does not contain every input block required by 2P-P^2",
    )
    input_to_output = CuArray(Int32[lookup[key] for key in input_keys])

    CUDA.synchronize()
    branch_start = time_ns()
    if branch == :hole
        @cuda threads=threads blocks=cld(length(square_values), threads) scale_blocks_kernel!(
            square_values, -1.0, length(square_values),
        )
        @cuda threads=threads blocks=cld(length(device.values), threads) add_mapped_blocks_kernel!(
            square_values, device.values, input_to_output, 2.0, elements_per_block,
            length(device.colind),
        )
    end
    @cuda threads=threads blocks=cld(length(square_values), threads) prune_elements_kernel!(
        square_values, cutoff, length(square_values),
    )
    output_norms = CUDA.zeros(Float64, length(output_keys))
    @cuda threads=threads blocks=cld(length(output_keys), threads) block_norm_kernel!(
        output_norms, square_values, device.block_size, length(output_keys),
    )
    CUDA.synchronize()
    branch_prune_time = (time_ns() - branch_start) / 1e9

    metadata_start = time_ns()
    keep = Array(output_norms) .> 0.0
    retained = findall(keep)
    retained_keys = output_keys[retained]
    structure = structure_from_keys(retained_keys, device.nblockrows)
    metadata_time = (time_ns() - metadata_start) / 1e9

    gather_start = time_ns()
    retained_device = CuArray(Int32.(retained))
    retained_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                                  length(retained))
    @cuda threads=threads blocks=cld(length(retained_values), threads) gather_blocks_kernel!(
        retained_values, square_values, retained_device, elements_per_block,
        length(retained),
    )
    retained_norms = CUDA.zeros(Float64, length(retained))
    @cuda threads=threads blocks=cld(length(retained), threads) block_norm_kernel!(
        retained_norms, retained_values, device.block_size, length(retained),
    )
    next_device = CuBlockCSR(
        device.nblockrows, device.block_size, CuArray(structure.rowptr),
        CuArray(structure.rowind), CuArray(structure.colind), retained_values,
        retained_norms,
    )
    CUDA.synchronize()
    gather_time = (time_ns() - gather_start) / 1e9
    return next_device, retained_keys, branch_prune_time, metadata_time, gather_time
end

function device_to_sparse(device::CuBlockCSR, N::Int)
    keys = collect(zip(Array(device.rowind), Array(device.colind)))
    return blocks_to_sparse(keys, Array(device.values), N, device.block_size)
end

"Device-resident block SP2; only structural metadata is rebuilt on the host."
function resident_sp2_continuation(P0::SparseMatrixCSC{Float64,Int}, block_size::Int,
                                   particles::Int, tau::Float64, steps::Int,
                                   output_cutoff::Float64, output::String;
                                   dense_projector::Union{Nothing,Matrix{Float64}}=nothing)
    device = cu_block_csr(host_block_csr(P0, block_size))
    open(joinpath(output, "resident_sp2_history.csv"), "w") do io
        println(io, "iteration,branch,trace,relative_trace_error,idempotency_defect,input_blocks,candidates,screened_fraction,output_blocks,structure_to_host_time_s,schedule_build_time_s,schedule_and_multiply_time_s,multiply_kernel_time_s,gpu_branch_prune_time_s,metadata_compaction_time_s,gpu_gather_time_s,metrics_time_s,total_iteration_time_s,gpu_free_memory_bytes")
        for iteration in 1:steps
            total_start = time_ns()
            metrics_start = time_ns()
            before = device_block_metrics(device)
            metrics_time = (time_ns() - metrics_start) / 1e9
            branch = before.trace >= particles ? :square : :hole

            structure_start = time_ns()
            structure = host_structure(device)
            host_norms = Array(device.norms)
            structure_to_host_time = (time_ns() - structure_start) / 1e9

            schedule_start = time_ns()
            schedule = grouped_output_schedule(structure)
            scores = [host_norms[left] * host_norms[right]
                      for (left, right) in zip(schedule.left, schedule.right)]
            screened_fraction = count(<(tau), scores) / max(length(scores), 1)
            schedule_build_time = (time_ns() - schedule_start) / 1e9

            multiply_start = time_ns()
            square_values, multiply_kernel_time = grouped_multiply(device, schedule, tau)
            schedule_and_multiply_time = (time_ns() - multiply_start) / 1e9

            device, _, branch_prune_time, metadata_time, gather_time =
                branch_prune_and_compact(
                    device, square_values, schedule.keys, branch, output_cutoff,
                )
            metrics_start = time_ns()
            after = device_block_metrics(device)
            metrics_time += (time_ns() - metrics_start) / 1e9
            idempotency = abs(after.squared_norm - after.trace) / particles
            total_iteration_time = (time_ns() - total_start) / 1e9
            println(io, join((
                iteration, branch, after.trace, abs(after.trace - particles) / particles,
                idempotency, length(structure.colind), length(schedule.left),
                screened_fraction, length(device.colind), structure_to_host_time,
                schedule_build_time, schedule_and_multiply_time, multiply_kernel_time,
                branch_prune_time, metadata_time, gather_time, metrics_time,
                total_iteration_time, CUDA.free_memory(),
            ), ','))
            flush(io)
        end
    end

    final_metrics = device_block_metrics(device)
    final_sparse = isnothing(dense_projector) ? nothing : device_to_sparse(device, size(P0, 1))
    dense_relative_error = isnothing(dense_projector) ? NaN :
        norm(Matrix(final_sparse) - dense_projector) / norm(dense_projector)
    dense_density_error = isnothing(dense_projector) ? NaN :
        maximum(abs, diag(final_sparse) - diag(dense_projector))
    open(joinpath(output, "resident_sp2_summary.toml"), "w") do io
        TOML.print(io, Dict(
            "matrix_dimension" => size(P0, 1),
            "particles" => particles,
            "iterations" => steps,
            "spamm_tau" => tau,
            "output_cutoff" => output_cutoff,
            "trace" => final_metrics.trace,
            "relative_trace_error" => abs(final_metrics.trace - particles) / particles,
            "trace_idempotency_defect" =>
                abs(final_metrics.squared_norm - final_metrics.trace) / particles,
            "final_blocks" => length(device.colind),
            "dense_validation_enabled" => !isnothing(dense_projector),
            "dense_projector_relative_error" => dense_relative_error,
            "dense_density_max_abs_error" => dense_density_error,
            "gpu_resident_block_values" => true,
            "gpu_branch_and_element_pruning" => true,
            "host_structural_metadata_only" => true,
        ))
    end
    return device
end

"SP2 trajectory with candidate scheduling, branch formation, and CSR rebuild on CUDA."
function gpu_pool_bytes(function_name)
    value = function_name()
    return ismissing(value) ? -1 : Int(value)
end

function gpu_scheduled_sp2_continuation(
    P0::SparseMatrixCSC{Float64,Int}, block_size::Int, particles::Int,
    tau::Float64, steps::Int, output_cutoff::Float64, output::String;
    dense_projector::Union{Nothing,Matrix{Float64}}=nothing,
    sort_workspace::CandidateSortWorkspace=CandidateSortWorkspace(),
)
    device = cu_block_csr(host_block_csr(P0, block_size))
    open(joinpath(output, "gpu_scheduled_sp2_history.csv"), "w") do io
        println(io, "iteration,branch,trace,relative_trace_error,idempotency_defect,input_blocks,compatible_candidates,retained_candidates,output_union_blocks,retained_output_blocks,candidate_count_time_s,candidate_scan_time_s,candidate_write_time_s,sort_time_s,group_time_s,union_time_s,fused_multiply_time_s,gpu_compaction_time_s,metrics_time_s,total_iteration_time_s,gpu_pool_used_bytes,gpu_pool_reserved_bytes,gpu_iteration_peak_used_bytes,gpu_free_memory_bytes")
        for iteration in 1:steps
            total_start = time_ns()
            iteration_peak_used = gpu_pool_bytes(CUDA.used_memory)
            metrics_start = time_ns()
            before = device_block_metrics(device)
            metrics_time = (time_ns() - metrics_start) / 1e9
            branch = before.trace >= particles ? :square : :hole
            input_blocks = length(device.colind)

            candidates = gpu_prescreened_candidates(device, tau)
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))
            # Compatible count is evaluated from CSR degrees on-device without
            # materializing rejected records; retain the screened count directly.
            degree_counts = CUDA.zeros(Int64, input_blocks)
            @cuda threads=256 blocks=cld(input_blocks, 256) compatible_candidate_count_kernel!(
                degree_counts, device.rowptr, device.colind, input_blocks,
            )
            compatible_candidates = Int(sum(degree_counts))

            sorted = gpu_sort_candidates(candidates, sort_workspace)
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))
            square_schedule = gpu_group_sorted_candidates(sorted)
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))
            branch_union = gpu_branch_union(device, square_schedule, branch)
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))
            output_values, fused_multiply_time = fused_branch_multiply(
                device, square_schedule, branch_union, branch, output_cutoff,
            )
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))
            compact_start = time_ns()
            device = gpu_compact_output(device, branch_union.keys, output_values)
            gpu_compaction_time = (time_ns() - compact_start) / 1e9
            iteration_peak_used = max(iteration_peak_used,
                                      gpu_pool_bytes(CUDA.used_memory))

            metrics_start = time_ns()
            after = device_block_metrics(device)
            metrics_time += (time_ns() - metrics_start) / 1e9
            idempotency = abs(after.squared_norm - after.trace) / particles
            total_iteration_time = (time_ns() - total_start) / 1e9
            println(io, join((
                iteration, branch, after.trace,
                abs(after.trace - particles) / particles, idempotency,
                input_blocks, compatible_candidates, candidates.retained,
                length(branch_union.keys), length(device.colind),
                candidates.count_time, candidates.scan_time, candidates.write_time,
                sorted.elapsed, square_schedule.elapsed, branch_union.elapsed,
                fused_multiply_time, gpu_compaction_time, metrics_time,
                total_iteration_time, gpu_pool_bytes(CUDA.used_memory),
                gpu_pool_bytes(CUDA.cached_memory), iteration_peak_used,
                CUDA.free_memory(),
            ), ','))
            flush(io)
        end
    end

    final_metrics = device_block_metrics(device)
    final_sparse = isnothing(dense_projector) ? nothing :
                   device_to_sparse(device, size(P0, 1))
    dense_relative_error = isnothing(dense_projector) ? NaN :
        norm(Matrix(final_sparse) - dense_projector) / norm(dense_projector)
    dense_density_error = isnothing(dense_projector) ? NaN :
        maximum(abs, diag(final_sparse) - diag(dense_projector))
    eigenvalue_minimum = NaN
    eigenvalue_maximum = NaN
    intruder_eigenvalues = -1
    if !isnothing(final_sparse)
        eigenvalues = eigvals(Symmetric(Matrix(final_sparse + final_sparse') / 2))
        eigenvalue_minimum = minimum(eigenvalues)
        eigenvalue_maximum = maximum(eigenvalues)
        intruder_eigenvalues = count(value -> 0.1 < value < 0.9, eigenvalues)
        dense_relative_error <= 1e-6 || error(
            "dense projector relative error exceeds 1e-6: $dense_relative_error",
        )
        dense_density_error <= 1e-6 || error(
            "dense diagonal error exceeds 1e-6: $dense_density_error",
        )
        eigenvalue_minimum >= -1e-6 || error(
            "projector eigenvalue below tolerance: $eigenvalue_minimum",
        )
        eigenvalue_maximum <= 1 + 1e-6 || error(
            "projector eigenvalue above tolerance: $eigenvalue_maximum",
        )
        intruder_eigenvalues == 0 || error(
            "projector has $intruder_eigenvalues eigenvalues in (0.1,0.9)",
        )
    end
    open(joinpath(output, "gpu_scheduled_sp2_summary.toml"), "w") do io
        TOML.print(io, Dict(
            "matrix_dimension" => size(P0, 1),
            "particles" => particles,
            "iterations" => steps,
            "spamm_tau" => tau,
            "output_cutoff" => output_cutoff,
            "trace" => final_metrics.trace,
            "relative_trace_error" => abs(final_metrics.trace - particles) / particles,
            "trace_idempotency_defect" =>
                abs(final_metrics.squared_norm - final_metrics.trace) / particles,
            "final_blocks" => length(device.colind),
            "dense_validation_enabled" => !isnothing(dense_projector),
            "dense_projector_relative_error" => dense_relative_error,
            "dense_density_max_abs_error" => dense_density_error,
            "projector_eigenvalue_minimum" => eigenvalue_minimum,
            "projector_eigenvalue_maximum" => eigenvalue_maximum,
            "projector_intruder_eigenvalues_0p1_0p9" => intruder_eigenvalues,
            "gpu_prescreened_candidates" => true,
            "gpu_sorted_grouped_schedule" => true,
            "gpu_fused_sp2_branch" => true,
            "gpu_output_pruning_and_csr_rebuild" => true,
            "host_structural_metadata" => false,
            "sort_workspace_reused" => true,
            "memory_telemetry" => "CUDA pool used/reserved plus sampled iteration peak",
        ))
    end
    return device
end

function compatible_candidate_count_kernel!(counts, rowptr, colind,
                                            number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        k = Int(colind[left])
        counts[left] = Int64(rowptr[k + 1]) - Int64(rowptr[k])
    end
    return
end

if !SPAMM_LIBRARY_ONLY
hopping_radius = parse(Int, get(ENV, "SPAMM_HOPPING_RADIUS", "1"))
hopping_decay_length = parse(Float64, get(ENV, "SPAMM_HOPPING_DECAY_LENGTH", "2.0"))
hopping_radius >= 1 || error("SPAMM_HOPPING_RADIUS must be positive")
hopping_decay_length > 0 || error("SPAMM_HOPPING_DECAY_LENGTH must be positive")
A = isnothing(V2) ? density_like_matrix(N) : fixed_h_sp2_snapshot(
    N, V2, snapshot_steps, hopping_radius, hopping_decay_length,
)
host = host_block_csr(A, block_size)
device = cu_block_csr(host)
host_norms = [norm(view(host.values, :, :, block)) for block in axes(host.values, 3)]
device_norms = Array(device.norms)
norm_error = maximum(abs, host_norms .- device_norms)
norm_error <= 100eps(Float64) * maximum(host_norms) || error("GPU block norms disagree with host")
offsets = candidate_offsets(host)
thresholds = (0.0, 1e-6, 1e-4, 1e-3, 1e-2)
grouped_schedule = grouped_output_schedule(host)
length(grouped_schedule.left) == offsets[end] - 1 || error(
    "grouped candidate schedule does not contain every compatible product",
)

production_only = lowercase(get(ENV, "SPAMM_PRODUCTION_ONLY", "false")) in ("1", "true", "yes")

# Experimental all-GPU symbolic scheduler for P^2 only.  It deliberately remains
# outside resident_sp2_continuation until the hole-branch union has been added.
scheduler_thresholds = production_only ?
    (parse(Float64, get(ENV, "SPAMM_TAU", "1e-6")),) : thresholds
run_scheduler_diagnostic = lowercase(get(
    ENV, "SPAMM_RUN_SCHEDULER_DIAGNOSTIC", production_only ? "false" : "true",
)) in ("1", "true", "yes")
if run_scheduler_diagnostic
open(joinpath(output, "gpu_scheduler.csv"), "w") do scheduler_io
    println(scheduler_io, "tau,compatible_candidates,retained_candidates,output_blocks,count_time_s,scan_time_s,write_time_s,sort_time_s,group_time_s,multiply_time_s,candidates_match,groups_match,product_max_abs_difference,product_relative_difference")
    sort_workspace = CandidateSortWorkspace()
    for tau in scheduler_thresholds
        candidates = gpu_prescreened_candidates(device, tau)

        # Warm up sort/group before recording the representative run.
        warm_sorted = gpu_sort_candidates(candidates, sort_workspace)
        warm_grouped = gpu_group_sorted_candidates(warm_sorted)
        grouped_multiply_device(device, warm_grouped)
        # One post-warmup measurement here avoids retaining several potentially
        # gigabyte-sized output tensors.  The standalone sort benchmark records
        # repeated sort measurements without allocating product blocks.
        sorted = gpu_sort_candidates(candidates, sort_workspace)
        sort_time = sorted.elapsed
        device_schedule = gpu_group_sorted_candidates(sorted)
        group_time = device_schedule.elapsed
        device_product_values, multiply_time =
            grouped_multiply_device(device, device_schedule)

        candidates_match = true
        groups_match = true
        product_max_abs_difference = NaN
        product_relative_difference = NaN
        if !production_only
            reference = host_candidate_reference(host, host_norms, tau)
            reference_records = sort!([
                ((UInt64(candidate[3]) << 32) | UInt64(candidate[4]),
                 candidate[1], candidate[2])
                for candidate in reference if candidate[6]
            ])
            device_records = sort!(collect(zip(
                Array(candidates.keys), Array(candidates.left), Array(candidates.right),
            )))
            candidates_match = device_records == reference_records

            reference_schedule = grouped_output_schedule(host)
            reference_values, _ = grouped_multiply(device, reference_schedule, tau)
            reference_product = blocks_to_sparse(
                reference_schedule.keys, Array(reference_values), N, block_size,
            )
            device_keys = decode_output_key.(Array(device_schedule.keys))
            device_product = blocks_to_sparse(
                device_keys, Array(device_product_values), N, block_size,
            )
            difference = reference_product - device_product
            product_max_abs_difference = maximum(abs, difference)
            product_relative_difference = norm(difference) / max(norm(reference_product), eps())

            expected_keys = sort!(unique!(first.(reference_records)))
            expected_counts = Dict{UInt64,Int}()
            for record in reference_records
                expected_counts[record[1]] = get(expected_counts, record[1], 0) + 1
            end
            actual_keys = Array(device_schedule.keys)
            actual_rowptr = Array(device_schedule.rowptr)
            groups_match = actual_keys == expected_keys &&
                all(actual_rowptr[index + 1] - actual_rowptr[index] ==
                    expected_counts[actual_keys[index]] for index in eachindex(actual_keys))

            candidates_match || error("GPU pre-screened candidates disagree at tau=$tau")
            groups_match || error("GPU run-length groups disagree at tau=$tau")
            product_relative_difference <= 1e-12 || error(
                "GPU sorted/grouped product disagrees at tau=$tau: relative difference=$product_relative_difference",
            )
        end

        compatible_candidates = offsets[end] - 1
        println(scheduler_io, join((
            tau, compatible_candidates, candidates.retained,
            length(device_schedule.keys), candidates.count_time,
            candidates.scan_time, candidates.write_time, sort_time, group_time,
            multiply_time, candidates_match, groups_match,
            product_max_abs_difference, product_relative_difference,
        ), ','))
        flush(scheduler_io)
    end
end
end

if !production_only
exact = A * A
exact_norm = norm(exact)
open(joinpath(output, "summary.csv"), "w") do io
    println(io, "tau,N,block_size,input_nnz,input_blocks,candidates,kept,screened,screened_fraction,candidate_kernel_time_s,compaction_time_s,atomic_multiply_time_s,grouped_multiply_time_s,grouped_speedup,output_blocks,error_bound,actual_error,relative_error,atomic_grouped_max_abs_difference,bound_satisfied,max_block_norm_error,host_gpu_match")
    for tau in thresholds
        gpu = gpu_candidates(device, offsets, tau)
        reference = host_candidate_reference(host, host_norms, tau)
        gpu_tuples = collect(zip(
            Array(gpu.left_index), Array(gpu.right_index), Array(gpu.output_row),
            Array(gpu.output_column), Array(gpu.scores), Array(gpu.keep),
        ))
        indices_match = all(gpu_tuples[index][1:4] == reference[index][1:4]
                            for index in eachindex(reference))
        scores_match = all(isapprox(gpu_tuples[index][5], reference[index][5];
                                    rtol=10eps(Float64), atol=0.0)
                           for index in eachindex(reference))
        masks_match = all(gpu_tuples[index][6] == reference[index][6]
                          for index in eachindex(reference))
        host_gpu_match = indices_match && scores_match && masks_match
        host_gpu_match || error("GPU candidate construction/screening mismatch for tau=$tau")
        kept = count(Array(gpu.keep))
        screened = gpu.number_of_candidates - kept
        output_keys, output_index = output_block_indices(
            Array(gpu.output_row), Array(gpu.output_column),
        )
        compact = compact_candidates(gpu, output_index)
        compact.retained == kept || error(
            "GPU compaction count mismatch for tau=$tau: $(compact.retained) != $kept",
        )
        atomic_runs = [multiply_accumulate(device, compact, length(output_keys))
                       for _ in 1:7]
        output_values = atomic_runs[end][1]
        multiply_time = median(last.(atomic_runs))
        approximate = blocks_to_sparse(output_keys, Array(output_values), N, block_size)
        grouped_runs = [grouped_multiply(device, grouped_schedule, tau) for _ in 1:7]
        grouped_values = grouped_runs[end][1]
        grouped_time = median(last.(grouped_runs))
        grouped_approximate = blocks_to_sparse(
            grouped_schedule.keys, Array(grouped_values), N, block_size,
        )
        atomic_grouped_difference = maximum(abs, approximate - grouped_approximate)
        atomic_grouped_difference <= 1e-12 || error(
            "atomic and grouped GPU products disagree for tau=$tau: $atomic_grouped_difference",
        )
        actual_error = norm(exact - approximate)
        error_bound = conservative_error_bound(reference)
        tolerance_slack = 100eps(Float64) * max(exact_norm, error_bound, 1.0)
        bound_satisfied = actual_error <= error_bound + tolerance_slack
        bound_satisfied || error(
            "GPU SpAMM error bound violated for tau=$tau: $actual_error > $error_bound",
        )
        println(io, join((
            tau, N, block_size, nnz(A), length(host.colind), gpu.number_of_candidates,
            kept, screened, screened / max(gpu.number_of_candidates, 1), gpu.elapsed,
            compact.elapsed, multiply_time, grouped_time, multiply_time / grouped_time,
            length(output_keys), error_bound, actual_error, actual_error / exact_norm,
            atomic_grouped_difference, bound_satisfied, norm_error, host_gpu_match,
        ), ','))
        flush(io)
    end
end
end

open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "diagnostic" => "gpu_block_spamm_stages_1_to_5",
        "matrix_dimension" => N,
        "block_size" => block_size,
        "input_nnz" => nnz(A),
        "input_blocks" => length(host.colind),
        "input_source" => isnothing(V2) ? "synthetic_density_like" : "fixed_h_sp2_snapshot",
        "V2" => isnothing(V2) ? "not_applicable" : V2,
        "hopping_radius" => hopping_radius,
        "hopping_decay_length" => hopping_decay_length,
        "hopping_profile" => "-exp(-(r-1)/xi)*(1+V2*cos(2pi*tau*(i+r/2)))",
        "snapshot_sp2_steps" => isnothing(V2) ? 0 : snapshot_steps,
        "device_name" => CUDA.name(CUDA.device()),
        "implemented_stages" => [
            "GPU block-CSR storage and Frobenius norms",
            "GPU candidate triple construction",
            "GPU block-norm product screening",
            "GPU retained-candidate compaction",
            "GPU leaf multiplication and atomic output-block accumulation",
            "GPU output-grouped multiplication without atomics",
        ],
        "multiplication_implemented" => true,
        "production_only" => production_only,
    ))
end

full_sp2_steps = parse(Int, get(ENV, "SPAMM_FULL_SP2_STEPS", "0"))
if full_sp2_steps > 0
    isnothing(V2) && error("SPAMM_FULL_SP2_STEPS requires a fixed-H V2 input")
    spamm_tau = parse(Float64, get(ENV, "SPAMM_TAU", "1e-6"))
    output_cutoff = parse(Float64, get(ENV, "SPAMM_OUTPUT_CUTOFF", "1e-8"))
    validate_reference = lowercase(get(ENV, "SPAMM_VALIDATE_REFERENCE", N <= 1024 ? "true" : "false")) in
                         ("1", "true", "yes")
    dense_projector = nothing
    if lowercase(get(ENV, "SPAMM_VALIDATE_DENSE", N <= 1024 ? "true" : "false")) in
       ("1", "true", "yes")
        N <= 4096 || error("dense-projector validation is restricted to N <= 4096")
        H_sparse, _ = fixed_hamiltonian(
            N, V2, hopping_radius, hopping_decay_length,
        )
        H = Symmetric(Matrix(H_sparse))
        eig = eigen(H)
        occupied = view(eig.vectors, :, 1:(N ÷ 2))
        dense_projector = Matrix(occupied * occupied')
    end
    gpu_scheduler_mode = lowercase(get(ENV, "SPAMM_GPU_SCHEDULER_MODE", "false")) in
                         ("1", "true", "yes")
    resident_mode = lowercase(get(ENV, "SPAMM_RESIDENT_MODE", "false")) in
                    ("1", "true", "yes")
    if gpu_scheduler_mode
        shared_sort_workspace = CandidateSortWorkspace()
        warmup_gpu = lowercase(get(ENV, "SPAMM_WARMUP_GPU_SCHEDULED", "false")) in
                     ("1", "true", "yes")
        if warmup_gpu
            # Warm-up exists only to compile and initialize the scheduled GPU
            # kernels.  Running the complete production trajectory here both
            # duplicates the benchmark and can fill the CUDA pool before the
            # measured cycles begin.  Two SP2 updates exercise the initial
            # square/hole alternation for the fixed-H validation inputs.
            warmup_steps = parse(Int, get(ENV, "SPAMM_WARMUP_STEPS", "2"))
            1 <= warmup_steps <= full_sp2_steps ||
                error("SPAMM_WARMUP_STEPS must lie in 1:SPAMM_FULL_SP2_STEPS")
            warmup_output = joinpath(output, "warmup")
            mkpath(warmup_output)
            gpu_scheduled_sp2_continuation(
                A, block_size, N ÷ 2, spamm_tau, warmup_steps, output_cutoff,
                warmup_output; dense_projector=nothing,
                sort_workspace=shared_sort_workspace,
            )
            GC.gc(true)
            CUDA.reclaim()
        end
        trajectory_cycles = parse(Int, get(ENV, "SPAMM_TRAJECTORY_CYCLES", "1"))
        trajectory_cycles >= 1 || error("SPAMM_TRAJECTORY_CYCLES must be positive")
        open(joinpath(output, "trajectory_cycles.csv"), "w") do cycle_io
            println(cycle_io, "cycle,elapsed_time_s,pool_used_before_bytes,pool_used_after_bytes,pool_reserved_before_bytes,pool_reserved_after_bytes,free_before_bytes,free_after_bytes,sort_workspace_capacity")
            for cycle in 1:trajectory_cycles
                cycle_output = trajectory_cycles == 1 ? output :
                               joinpath(output, "cycle_$(lpad(cycle, 2, '0'))")
                mkpath(cycle_output)
                used_before = gpu_pool_bytes(CUDA.used_memory)
                reserved_before = gpu_pool_bytes(CUDA.cached_memory)
                free_before = CUDA.free_memory()
                result_device = nothing
                elapsed = @elapsed begin
                    result_device = gpu_scheduled_sp2_continuation(
                        A, block_size, N ÷ 2, spamm_tau, full_sp2_steps,
                        output_cutoff, cycle_output; dense_projector,
                        sort_workspace=shared_sort_workspace,
                    )
                end
                result_device = nothing
                GC.gc(true)
                used_after = gpu_pool_bytes(CUDA.used_memory)
                reserved_after = gpu_pool_bytes(CUDA.cached_memory)
                free_after = CUDA.free_memory()
                println(cycle_io, join((
                    cycle, elapsed, used_before, used_after, reserved_before,
                    reserved_after, free_before, free_after,
                    shared_sort_workspace.capacity,
                ), ','))
                flush(cycle_io)
            end
        end
    elseif resident_mode
        resident_sp2_continuation(
            A, block_size, N ÷ 2, spamm_tau, full_sp2_steps, output_cutoff, output;
            dense_projector,
        )
    else
        validated_sp2_continuation(
            A, block_size, N ÷ 2, spamm_tau, full_sp2_steps, output_cutoff, output;
            validate_reference, dense_projector,
        )
    end
end

println("GPU block-SpAMM stages 1--5 verified: $output")
end # !SPAMM_LIBRARY_ONLY
