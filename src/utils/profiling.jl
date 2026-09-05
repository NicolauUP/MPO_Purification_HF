"""Run one operation with optional synchronized timing instrumentation.

When `callback` is `nothing`, this performs only the operation: no clock,
allocation measurement, device synchronization, or event construction occurs.
The helper is private because persistence and presentation belong in benchmark
scripts rather than the scientific kernels.
"""
function _profiled_operation(
    operation,
    phase::Symbol;
    callback::Union{Nothing,Function}=nothing,
    synchronize::Function=() -> nothing,
    metadata::NamedTuple=(;),
)
    isnothing(callback) && return operation()
    synchronize()
    timing = @timed begin
        value = operation()
        synchronize()
        value
    end
    callback(merge(metadata, (
        phase=phase,
        elapsed_time_s=timing.time,
        allocations_bytes=timing.bytes,
        gc_time_s=timing.gctime,
    )))
    return timing.value
end
