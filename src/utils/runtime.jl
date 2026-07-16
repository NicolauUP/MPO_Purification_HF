const _GC_POLICIES = (:automatic, :forced, :periodic, :threshold)

function _validate_gc_policy(
    gc_policy::Symbol,
    gc_period::Integer,
    gc_threshold_bytes::Integer,
)
    gc_policy in _GC_POLICIES || throw(ArgumentError(
        "gc_policy must be one of $(_GC_POLICIES), got $gc_policy",
    ))
    gc_period > 0 || throw(ArgumentError("gc_period must be positive, got $gc_period"))
    gc_threshold_bytes > 0 || throw(ArgumentError(
        "gc_threshold_bytes must be positive, got $gc_threshold_bytes",
    ))
    return nothing
end

function _should_collect_garbage(
    iteration::Integer,
    gc_policy::Symbol,
    gc_period::Integer,
    gc_threshold_bytes::Integer,
)
    if gc_policy == :automatic
        return false
    elseif gc_policy == :forced
        return true
    elseif gc_policy == :periodic
        return iteration % gc_period == 0
    end
    return Base.gc_live_bytes() >= gc_threshold_bytes
end

"""
    maybe_collect_garbage!(iteration; gc_policy=:automatic, gc_period=10,
                            gc_threshold_bytes=1 << 30)

Apply the requested host-GC policy at an iteration boundary. `:automatic`
performs no explicit collection and leaves the decision to Julia. `:forced`
collects at every boundary, `:periodic` every `gc_period` boundaries, and
`:threshold` when Julia reports at least `gc_threshold_bytes` live heap bytes.
The threshold uses Julia's live-heap counter; it is a host-memory trigger, not
a GPU-memory measurement.
"""
function maybe_collect_garbage!(
    iteration::Integer;
    gc_policy::Symbol=:automatic,
    gc_period::Integer=10,
    gc_threshold_bytes::Integer=1 << 30,
)
    _validate_gc_policy(gc_policy, gc_period, gc_threshold_bytes)
    _should_collect_garbage(iteration, gc_policy, gc_period, gc_threshold_bytes) && GC.gc(true)
    return nothing
end
