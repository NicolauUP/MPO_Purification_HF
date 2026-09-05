module BlockSp2Engine

using CUDA
using LinearAlgebra
using SparseArrays
using TOML

export HostBlockCSR, CuBlockCSR
export host_to_cublock_csr, cublock_to_host_csr, purify_gpu_block_sp2!

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

function _sparse_to_hostblock(A::SparseMatrixCSC{Float64,Int}, block_size::Int)
    size(A, 1) == size(A, 2) || throw(ArgumentError("matrix must be square"))
    size(A, 1) % block_size == 0 ||
        throw(ArgumentError("matrix dimension must be divisible by block_size"))
    1 <= block_size <= 32 || throw(ArgumentError("block_size must lie in 1:32"))

    nblockrows = size(A, 1) ÷ block_size
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

function _block_norm_kernel!(norms, values, block_size, number_of_blocks)
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

function host_to_cublock_csr(host::HostBlockCSR)
    rowptr = CuArray(host.rowptr)
    rowind = CuArray(host.rowind)
    colind = CuArray(host.colind)
    values = CuArray(host.values)
    norms = CUDA.zeros(Float64, length(host.colind))
    if !isempty(host.colind)
        threads = 256
        @cuda threads=threads blocks=cld(length(host.colind), threads) _block_norm_kernel!(
            norms, values, host.block_size, length(host.colind),
        )
        CUDA.synchronize()
    end
    return CuBlockCSR(host.nblockrows, host.block_size, rowptr, rowind, colind,
                      values, norms)
end

host_to_cublock_csr(A::SparseMatrixCSC{Float64,Int}, block_size::Int) =
    host_to_cublock_csr(_sparse_to_hostblock(A, block_size))

function cublock_to_host_csr(device::CuBlockCSR)
    return HostBlockCSR(
        device.nblockrows,
        device.block_size,
        Array(device.rowptr),
        Array(device.rowind),
        Array(device.colind),
        Array(device.values),
    )
end

function cublock_to_host_csr(device::CuBlockCSR, matrix_dimension::Int)
    matrix_dimension <= device.nblockrows * device.block_size ||
        throw(ArgumentError("matrix_dimension exceeds the block storage extent"))
    host = cublock_to_host_csr(device)
    keys = collect(zip(host.rowind, host.colind))
    rows = Int[]
    columns = Int[]
    entries = Float64[]
    for (position, (block_row, block_column)) in enumerate(keys)
        row_offset = (Int(block_row) - 1) * host.block_size
        column_offset = (Int(block_column) - 1) * host.block_size
        for column in 1:host.block_size, row in 1:host.block_size
            global_row = row_offset + row
            global_column = column_offset + column
            global_row <= matrix_dimension && global_column <= matrix_dimension || continue
            value = host.values[row, column, position]
            iszero(value) && continue
            push!(rows, global_row)
            push!(columns, global_column)
            push!(entries, value)
        end
    end
    return sparse(rows, columns, entries, matrix_dimension, matrix_dimension)
end

function _screened_candidate_count_kernel!(counts, rowptr, colind, norms, tau,
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

function _screened_candidate_offset_kernel!(offsets, inclusive_counts,
                                             number_of_left_blocks)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index == 1 && (offsets[1] = Int64(1))
    index <= number_of_left_blocks &&
        (offsets[index + 1] = Int64(inclusive_counts[index]) + Int64(1))
    return
end

function _screened_candidate_write_kernel!(keys, left_indices, right_indices,
                                           rowptr, rowind, colind, norms, offsets,
                                           tau, number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        i = UInt64(rowind[left])
        k = Int(colind[left])
        destination = Int(offsets[left])
        for right in Int(rowptr[k]):(Int(rowptr[k + 1]) - 1)
            if norms[left] * norms[right] >= tau
                keys[destination] = (i << 32) | UInt64(colind[right])
                left_indices[destination] = Int32(left)
                right_indices[destination] = Int32(right)
                destination += 1
            end
        end
    end
    return
end

function _gpu_prescreened_candidates(device::CuBlockCSR, tau::Float64)
    number_of_left_blocks = length(device.colind)
    number_of_left_blocks > 0 || error("cannot square an empty block matrix")
    counts = CUDA.zeros(Int32, number_of_left_blocks)
    threads = 256
    CUDA.synchronize()
    count_time = @elapsed begin
        @cuda threads=threads blocks=cld(number_of_left_blocks, threads) _screened_candidate_count_kernel!(
            counts, device.rowptr, device.colind, device.norms, tau,
            number_of_left_blocks,
        )
        CUDA.synchronize()
    end
    scan_time = @elapsed begin
        inclusive_counts = CUDA.cumsum(counts)
        offsets = CUDA.zeros(Int64, number_of_left_blocks + 1)
        @cuda threads=threads blocks=cld(number_of_left_blocks, threads) _screened_candidate_offset_kernel!(
            offsets, inclusive_counts, number_of_left_blocks,
        )
        CUDA.synchronize()
    end
    retained = Int(Array(view(offsets, number_of_left_blocks + 1:
                                      number_of_left_blocks + 1))[1] - 1)
    retained > 0 || error("SpAMM threshold removed every compatible block product")
    keys = CUDA.zeros(UInt64, retained)
    left = CUDA.zeros(Int32, retained)
    right = CUDA.zeros(Int32, retained)
    write_time = @elapsed begin
        @cuda threads=threads blocks=cld(number_of_left_blocks, threads) _screened_candidate_write_kernel!(
            keys, left, right, device.rowptr, device.rowind, device.colind,
            device.norms, offsets, tau, number_of_left_blocks,
        )
        CUDA.synchronize()
    end
    return (; keys, left, right, retained, count_time, scan_time, write_time)
end

mutable struct CandidateSortWorkspace
    size::Int
    keys::CuVector{UInt64}
    permutation::CuVector{Int32}
    left::CuVector{Int32}
    right::CuVector{Int32}
end

CandidateSortWorkspace() = CandidateSortWorkspace(
    0, CUDA.zeros(UInt64, 0), CUDA.zeros(Int32, 0),
    CUDA.zeros(Int32, 0), CUDA.zeros(Int32, 0),
)

function _ensure_sort_size!(workspace::CandidateSortWorkspace, required::Int)
    workspace.size == required && return workspace
    workspace.keys = CUDA.zeros(UInt64, required)
    workspace.permutation = CUDA.zeros(Int32, required)
    workspace.left = CUDA.zeros(Int32, required)
    workspace.right = CUDA.zeros(Int32, required)
    workspace.size = required
    return workspace
end

function _gather_sorted_payload_kernel!(sorted_left, sorted_right, left, right,
                                        permutation, number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates
        source = permutation[index]
        sorted_left[index] = left[source]
        sorted_right[index] = right[source]
    end
    return
end

function _gpu_sort_candidates(candidates, workspace::CandidateSortWorkspace)
    number_of_candidates = length(candidates.keys)
    _ensure_sort_size!(workspace, number_of_candidates)
    threads = 256
    CUDA.synchronize()
    elapsed = @elapsed begin
        # CUDA.sortperm! sorts its key input in-place. Never pass candidate keys:
        # later branch construction still relies on the original record ordering.
        copyto!(workspace.keys, candidates.keys)
        CUDA.sortperm!(workspace.permutation, workspace.keys; initialized=false)
        @cuda threads=threads blocks=cld(number_of_candidates, threads) _gather_sorted_payload_kernel!(
            workspace.left, workspace.right, candidates.left, candidates.right,
            workspace.permutation, number_of_candidates,
        )
        CUDA.synchronize()
    end
    return (; keys=workspace.keys, left=workspace.left, right=workspace.right, elapsed)
end

function _group_boundary_kernel!(boundaries, keys, number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates
        boundaries[index] = index == 1 || keys[index] != keys[index - 1]
    end
    return
end

function _group_scatter_kernel!(unique_keys, group_rowptr, keys, boundaries,
                                group_ids, number_of_candidates)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if index <= number_of_candidates && boundaries[index] != 0
        group = Int(group_ids[index])
        unique_keys[group] = keys[index]
        group_rowptr[group] = Int64(index)
    end
    return
end

function _group_finalize_kernel!(group_rowptr, number_of_groups,
                                 number_of_candidates)
    if blockIdx().x == 1 && threadIdx().x == 1
        group_rowptr[number_of_groups + 1] = Int64(number_of_candidates + 1)
    end
    return
end

function _gpu_group_sorted_candidates(sorted)
    number_of_candidates = length(sorted.keys)
    boundaries = CUDA.zeros(Int32, number_of_candidates)
    threads = 256
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=threads blocks=cld(number_of_candidates, threads) _group_boundary_kernel!(
            boundaries, sorted.keys, number_of_candidates,
        )
        group_ids = CUDA.cumsum(boundaries)
        CUDA.synchronize()
        number_of_groups = Int(Array(view(group_ids, number_of_candidates:
                                                    number_of_candidates))[1])
        unique_keys = CUDA.zeros(UInt64, number_of_groups)
        group_rowptr = CUDA.zeros(Int64, number_of_groups + 1)
        @cuda threads=threads blocks=cld(number_of_candidates, threads) _group_scatter_kernel!(
            unique_keys, group_rowptr, sorted.keys, boundaries, group_ids,
            number_of_candidates,
        )
        @cuda threads=1 blocks=1 _group_finalize_kernel!(
            group_rowptr, number_of_groups, number_of_candidates,
        )
        CUDA.synchronize()
    end
    return (; keys=unique_keys, rowptr=group_rowptr, left=sorted.left,
            right=sorted.right, elapsed)
end

function _pack_block_keys_kernel!(keys, rowind, colind, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    block <= number_of_blocks &&
        (keys[block] = (UInt64(rowind[block]) << 32) | UInt64(colind[block]))
    return
end

function _device_block_keys(device::CuBlockCSR)
    keys = CUDA.zeros(UInt64, length(device.colind))
    threads = 256
    @cuda threads=threads blocks=cld(length(keys), threads) _pack_block_keys_kernel!(
        keys, device.rowind, device.colind, length(keys),
    )
    return keys
end

function _gpu_unique_sorted_keys(keys)
    isempty(keys) && error("cannot uniquify an empty key list")
    sorted_keys = sort(copy(keys))
    boundaries = CUDA.zeros(Int32, length(sorted_keys))
    threads = 256
    @cuda threads=threads blocks=cld(length(sorted_keys), threads) _group_boundary_kernel!(
        boundaries, sorted_keys, length(sorted_keys),
    )
    group_ids = CUDA.cumsum(boundaries)
    CUDA.synchronize()
    number_of_groups = Int(Array(view(group_ids, length(group_ids):length(group_ids)))[1])
    unique_keys = CUDA.zeros(UInt64, number_of_groups)
    dummy_rowptr = CUDA.zeros(Int64, number_of_groups + 1)
    @cuda threads=threads blocks=cld(length(sorted_keys), threads) _group_scatter_kernel!(
        unique_keys, dummy_rowptr, sorted_keys, boundaries, group_ids,
        length(sorted_keys),
    )
    CUDA.synchronize()
    return unique_keys
end

function _exact_key_lookup_kernel!(source_to_union, source_keys, union_keys,
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

function _scatter_inverse_map_kernel!(union_to_source, source_to_union,
                                      number_of_source_keys)
    source = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if source <= number_of_source_keys
        destination = source_to_union[source]
        destination > 0 && (union_to_source[destination] = Int32(source))
    end
    return
end

function _identity_map_kernel!(mapping, number_of_entries)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    index <= number_of_entries && (mapping[index] = Int32(index))
    return
end

function _map_sorted_keys_to_union(source_keys, union_keys)
    source_to_union = CUDA.zeros(Int32, length(source_keys))
    union_to_source = CUDA.zeros(Int32, length(union_keys))
    threads = 256
    @cuda threads=threads blocks=cld(length(source_keys), threads) _exact_key_lookup_kernel!(
        source_to_union, source_keys, union_keys, length(source_keys),
        length(union_keys),
    )
    @cuda threads=threads blocks=cld(length(source_keys), threads) _scatter_inverse_map_kernel!(
        union_to_source, source_to_union, length(source_keys),
    )
    return union_to_source
end

function _gpu_branch_union(device::CuBlockCSR, square_schedule, branch::Symbol)
    CUDA.synchronize()
    elapsed = @elapsed begin
        input_keys = _device_block_keys(device)
        if branch == :square
            output_keys = square_schedule.keys
            output_to_square = CUDA.zeros(Int32, length(output_keys))
            @cuda threads=256 blocks=cld(length(output_keys), 256) _identity_map_kernel!(
                output_to_square, length(output_keys),
            )
            output_to_input = CUDA.zeros(Int32, length(output_keys))
        elseif branch == :hole
            output_keys = _gpu_unique_sorted_keys(vcat(square_schedule.keys, input_keys))
            output_to_square = _map_sorted_keys_to_union(square_schedule.keys, output_keys)
            output_to_input = _map_sorted_keys_to_union(input_keys, output_keys)
        else
            error("unsupported SP2 branch: $branch")
        end
        CUDA.synchronize()
    end
    return (; keys=output_keys, square=output_to_square, input=output_to_input,
            elapsed)
end

function _fused_branch_multiply_kernel!(output_values, input_values,
                                        square_rowptr, square_left, square_right,
                                        output_to_square, output_to_input,
                                        branch_sign, input_factor, output_cutoff,
                                        block_size, number_of_output_blocks)
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

function _fused_branch_multiply(device::CuBlockCSR, square_schedule,
                                branch_union, branch::Symbol,
                                output_cutoff::Float64)
    output_values = CUDA.zeros(Float64, device.block_size, device.block_size,
                               length(branch_union.keys))
    branch_sign = branch == :square ? 1.0 : -1.0
    input_factor = branch == :square ? 0.0 : 2.0
    CUDA.synchronize()
    elapsed = @elapsed begin
        @cuda threads=device.block_size^2 blocks=length(branch_union.keys) shmem=(
            2 * sizeof(Float64) * device.block_size^2
        ) _fused_branch_multiply_kernel!(
            output_values, device.values, square_schedule.rowptr,
            square_schedule.left, square_schedule.right, branch_union.square,
            branch_union.input, branch_sign, input_factor, output_cutoff,
            device.block_size, length(branch_union.keys),
        )
        CUDA.synchronize()
    end
    return output_values, elapsed
end

function _mark_retained_blocks_kernel!(keep, norms, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    block <= number_of_blocks &&
        (keep[block] = ifelse(norms[block] > 0.0, Int32(1), Int32(0)))
    return
end

function _compact_block_metadata_kernel!(keys_out, norms_out, keys_in, norms_in,
                                         keep, destinations, number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks && keep[block] != 0
        destination = Int(destinations[block])
        keys_out[destination] = keys_in[block]
        norms_out[destination] = norms_in[block]
    end
    return
end

function _compact_block_values_kernel!(values_out, values_in, keep, destinations,
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

function _decode_keys_and_count_rows_kernel!(rowind, colind, row_counts, keys,
                                             number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        key = keys[block]
        row = Int32(key >> 32)
        rowind[block] = row
        colind[block] = Int32(key & 0xffffffff)
        CUDA.@atomic row_counts[row] += Int32(1)
    end
    return
end

function _int32_rowptr_kernel!(rowptr, inclusive_counts, number_of_rows)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    row == 1 && (rowptr[1] = Int32(1))
    row <= number_of_rows &&
        (rowptr[row + 1] = Int32(inclusive_counts[row]) + Int32(1))
    return
end

function _gpu_compact_output(device::CuBlockCSR, output_keys, output_values)
    threads = 256
    number_of_blocks = length(output_keys)
    output_norms = CUDA.zeros(Float64, number_of_blocks)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) _block_norm_kernel!(
        output_norms, output_values, device.block_size, number_of_blocks,
    )
    keep = CUDA.zeros(Int32, number_of_blocks)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) _mark_retained_blocks_kernel!(
        keep, output_norms, number_of_blocks,
    )
    destinations = CUDA.cumsum(keep)
    CUDA.synchronize()
    retained = Int(Array(view(destinations, number_of_blocks:number_of_blocks))[1])
    retained > 0 || error("SP2 output contains no retained blocks")

    retained_keys = CUDA.zeros(UInt64, retained)
    retained_norms = CUDA.zeros(Float64, retained)
    retained_values = CUDA.zeros(Float64, device.block_size, device.block_size, retained)
    @cuda threads=threads blocks=cld(number_of_blocks, threads) _compact_block_metadata_kernel!(
        retained_keys, retained_norms, output_keys, output_norms, keep,
        destinations, number_of_blocks,
    )
    @cuda threads=threads blocks=cld(length(output_values), threads) _compact_block_values_kernel!(
        retained_values, output_values, keep, destinations, device.block_size^2,
        number_of_blocks,
    )
    rowind = CUDA.zeros(Int32, retained)
    colind = CUDA.zeros(Int32, retained)
    row_counts = CUDA.zeros(Int32, device.nblockrows)
    @cuda threads=threads blocks=cld(retained, threads) _decode_keys_and_count_rows_kernel!(
        rowind, colind, row_counts, retained_keys, retained,
    )
    inclusive_counts = CUDA.cumsum(row_counts)
    rowptr = CUDA.zeros(Int32, device.nblockrows + 1)
    @cuda threads=threads blocks=cld(device.nblockrows, threads) _int32_rowptr_kernel!(
        rowptr, inclusive_counts, device.nblockrows,
    )
    CUDA.synchronize()
    return CuBlockCSR(device.nblockrows, device.block_size, rowptr, rowind, colind,
                      retained_values, retained_norms)
end

function _block_metrics_kernel!(metrics, values, rowind, colind, block_size,
                                number_of_blocks)
    block = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if block <= number_of_blocks
        squared_norm = 0.0
        trace_value = 0.0
        diagonal_block = rowind[block] == colind[block]
        for column in 1:block_size, row in 1:block_size
            value = values[row, column, block]
            squared_norm += value * value
            diagonal_block && row == column && (trace_value += value)
        end
        CUDA.@atomic metrics[2] += squared_norm
        if diagonal_block
            CUDA.@atomic metrics[1] += trace_value
        end
    end
    return
end

function _device_block_metrics(device::CuBlockCSR)
    metrics = CUDA.zeros(Float64, 2)
    number_of_blocks = length(device.colind)
    threads = 256
    @cuda threads=threads blocks=cld(number_of_blocks, threads) _block_metrics_kernel!(
        metrics, device.values, device.rowind, device.colind, device.block_size,
        number_of_blocks,
    )
    CUDA.synchronize()
    result = Array(metrics)
    return (; trace=result[1], squared_norm=result[2])
end

function _compatible_candidate_count_kernel!(counts, rowptr, colind,
                                             number_of_left_blocks)
    left = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if left <= number_of_left_blocks
        k = Int(colind[left])
        counts[left] = Int64(rowptr[k + 1]) - Int64(rowptr[k])
    end
    return
end

function _gpu_pool_bytes(function_name)
    value = function_name()
    return ismissing(value) ? -1 : Int(value)
end

function _sp2_step(device::CuBlockCSR, particles::Int, tau::Float64,
                   output_cutoff::Float64, sort_workspace::CandidateSortWorkspace)
    total_start = time_ns()
    peak_used = _gpu_pool_bytes(CUDA.used_memory)
    metrics_start = time_ns()
    before = _device_block_metrics(device)
    metrics_time = (time_ns() - metrics_start) / 1e9
    branch = before.trace >= particles ? :square : :hole
    input_blocks = length(device.colind)

    candidates = _gpu_prescreened_candidates(device, tau)
    peak_used = max(peak_used, _gpu_pool_bytes(CUDA.used_memory))
    degree_counts = CUDA.zeros(Int64, input_blocks)
    @cuda threads=256 blocks=cld(input_blocks, 256) _compatible_candidate_count_kernel!(
        degree_counts, device.rowptr, device.colind, input_blocks,
    )
    compatible_candidates = Int(sum(degree_counts))
    sorted = _gpu_sort_candidates(candidates, sort_workspace)
    schedule = _gpu_group_sorted_candidates(sorted)
    branch_union = _gpu_branch_union(device, schedule, branch)
    output_values, multiply_time = _fused_branch_multiply(
        device, schedule, branch_union, branch, output_cutoff,
    )
    compact_start = time_ns()
    next_device = _gpu_compact_output(device, branch_union.keys, output_values)
    compact_time = (time_ns() - compact_start) / 1e9
    peak_used = max(peak_used, _gpu_pool_bytes(CUDA.used_memory))
    metrics_start = time_ns()
    after = _device_block_metrics(next_device)
    metrics_time += (time_ns() - metrics_start) / 1e9
    total_time = (time_ns() - total_start) / 1e9
    stats = (;
        branch, trace=after.trace,
        relative_trace_error=abs(after.trace - particles) / particles,
        idempotency_defect=abs(after.squared_norm - after.trace) / particles,
        input_blocks, compatible_candidates, retained_candidates=candidates.retained,
        output_union_blocks=length(branch_union.keys),
        retained_output_blocks=length(next_device.colind),
        candidate_count_time=candidates.count_time,
        candidate_scan_time=candidates.scan_time,
        candidate_write_time=candidates.write_time,
        sort_time=sorted.elapsed, group_time=schedule.elapsed,
        union_time=branch_union.elapsed, multiply_time, compact_time, metrics_time,
        total_time, pool_used=_gpu_pool_bytes(CUDA.used_memory),
        pool_reserved=_gpu_pool_bytes(CUDA.cached_memory), peak_used,
        free_memory=CUDA.free_memory(),
    )
    return next_device, stats
end

"""Purify `P0` with GPU-scheduled block-SpAMM SP2 and return a `CuBlockCSR`.

When `output` is supplied, iteration diagnostics and a final TOML summary are
written there. `P0` itself is not modified; the `!` denotes mutation of the
GPU-resident iterate and reusable sorting workspace during purification.
"""
function purify_gpu_block_sp2!(
    P0::SparseMatrixCSC{Float64,Int}, block_size::Int, particles::Int,
    tau::Float64, steps::Int, output_cutoff::Float64;
    output::Union{Nothing,AbstractString}=nothing,
)
    particles > 0 || throw(ArgumentError("particles must be positive"))
    steps >= 0 || throw(ArgumentError("steps must be nonnegative"))
    tau >= 0 || throw(ArgumentError("tau must be nonnegative"))
    output_cutoff >= 0 || throw(ArgumentError("output_cutoff must be nonnegative"))
    device = host_to_cublock_csr(P0, block_size)
    workspace = CandidateSortWorkspace()
    history_io = nothing
    if !isnothing(output)
        mkpath(output)
        history_io = open(joinpath(output, "gpu_scheduled_sp2_history.csv"), "w")
        println(history_io, "iteration,branch,trace,relative_trace_error,idempotency_defect,input_blocks,compatible_candidates,retained_candidates,output_union_blocks,retained_output_blocks,candidate_count_time_s,candidate_scan_time_s,candidate_write_time_s,sort_time_s,group_time_s,union_time_s,fused_multiply_time_s,gpu_compaction_time_s,metrics_time_s,total_iteration_time_s,gpu_pool_used_bytes,gpu_pool_reserved_bytes,gpu_iteration_peak_used_bytes,gpu_free_memory_bytes")
    end
    try
        for iteration in 1:steps
            device, stats = _sp2_step(device, particles, tau, output_cutoff, workspace)
            if !isnothing(history_io)
                println(history_io, join((
                    iteration, stats.branch, stats.trace, stats.relative_trace_error,
                    stats.idempotency_defect, stats.input_blocks,
                    stats.compatible_candidates, stats.retained_candidates,
                    stats.output_union_blocks, stats.retained_output_blocks,
                    stats.candidate_count_time, stats.candidate_scan_time,
                    stats.candidate_write_time, stats.sort_time, stats.group_time,
                    stats.union_time, stats.multiply_time, stats.compact_time,
                    stats.metrics_time, stats.total_time, stats.pool_used,
                    stats.pool_reserved, stats.peak_used, stats.free_memory,
                ), ','))
                flush(history_io)
            end
        end
    finally
        !isnothing(history_io) && close(history_io)
    end

    if !isnothing(output)
        metrics = _device_block_metrics(device)
        open(joinpath(output, "gpu_scheduled_sp2_summary.toml"), "w") do io
            TOML.print(io, Dict(
                "matrix_dimension" => size(P0, 1),
                "particles" => particles,
                "iterations" => steps,
                "spamm_tau" => tau,
                "output_cutoff" => output_cutoff,
                "trace" => metrics.trace,
                "relative_trace_error" => abs(metrics.trace - particles) / particles,
                "trace_idempotency_defect" =>
                    abs(metrics.squared_norm - metrics.trace) / particles,
                "final_blocks" => length(device.colind),
                "candidate_keys_copied_before_sortperm" => true,
            ))
        end
    end
    return device
end

end # module BlockSp2Engine
