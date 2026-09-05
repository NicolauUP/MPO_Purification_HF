"""Public physical-model and numerical-setting types.

These types deliberately sit in front of the established `Parameters1D` and
`ParametersSquare` implementation. The conversion functions at the bottom of
this file are lossless for the fields represented by the legacy parameter
types; numerical kernels continue to consume those legacy types unchanged.
"""

abstract type AbstractPhysicalModel end

function _validated_power_of_two_size(size::Integer, label::AbstractString)
    size >= 2 || throw(ArgumentError("$label must be at least 2, got $size"))
    size <= typemax(Int) || throw(ArgumentError("$label is too large for Int: $size"))
    value = Int(size)
    value & (value - 1) == 0 || throw(ArgumentError(
        "$label must be an exact power of two for the current QTT representation, got $size",
    ))
    return value
end

function _validated_filling(filling::Real, N::Integer)
    isfinite(filling) && 0.0 < filling < 1.0 || throw(ArgumentError(
        "filling must be finite and lie strictly between 0 and 1, got $filling",
    ))
    value = Float64(filling)
    particles = round(Int, N * value)
    0 < particles < N || throw(ArgumentError(
        "filling=$filling rounds to unsupported occupation Ne=$particles for N=$N",
    ))
    return value
end

function _validated_interaction(interaction)
    interaction isa Real || throw(ArgumentError(
        "interaction must be a real scalar; spatially dependent and complex interactions are not implemented",
    ))
    isfinite(interaction) || throw(ArgumentError("interaction must be finite, got $interaction"))
    return interaction
end

function _validated_optional_field(field, label::AbstractString)
    (isnothing(field) || field isa Function) || throw(ArgumentError(
        "$label must be nothing or a function",
    ))
    return field
end

function _validated_hopping_component(component, label::AbstractString)
    (component isa Number || component isa Function) || throw(ArgumentError(
        "$label must be a number or function",
    ))
    component isa Number && !isfinite(component) && throw(ArgumentError(
        "numeric $label must be finite, got $component",
    ))
    return component
end

"""
    ChainModel(; size, hopping, interaction, potential=nothing, seed=nothing,
               filling=0.5, boundary=:open)

Describe the physical open chain independently of QTT/MPO and SCF settings.
`size` is the physical number of sites and must be a power of two in the
currently supported QTT implementation.
"""
struct ChainModel{Tt,Tu,Tw,Ts} <: AbstractPhysicalModel
    size::Int
    hopping::Tt
    interaction::Tu
    potential::Tw
    seed::Ts
    filling::Float64
    boundary::Symbol
end

function ChainModel(
    ; size::Integer,
    hopping,
    interaction,
    potential=nothing,
    seed=nothing,
    filling::Real=0.5,
    boundary::Symbol=:open,
)
    boundary == :open || throw(ArgumentError(
        "the current ChainModel compatibility layer supports only boundary=:open, got $boundary",
    ))
    physical_size = _validated_power_of_two_size(size, "chain size")
    _validated_hopping_component(hopping, "chain hopping")
    _validated_interaction(interaction)
    _validated_optional_field(potential, "potential")
    _validated_optional_field(seed, "seed")
    physical_filling = _validated_filling(filling, physical_size)
    return ChainModel{typeof(hopping),typeof(interaction),typeof(potential),typeof(seed)}(
        physical_size, hopping, interaction, potential, seed, physical_filling, boundary,
    )
end

"""
    SquareModel(; size, hopping, interaction, potential=nothing, seed=nothing,
                filling=0.5, boundary=:open)

Describe an open two-dimensional lattice with power-of-two dimensions. Its
historical name is retained for compatibility: `size=(Nx, Ny)` may describe a
square when `Nx == Ny` or a rectangle otherwise. The public model is geometry
complete, whereas conversion to the current legacy MPO/dense core remains
available only for equal squares; call `legacy_parameters` to receive that
explicit capability error before an expensive calculation is started.
"""
struct SquareModel{Tt,Tu,Tw,Ts} <: AbstractPhysicalModel
    size::Tuple{Int,Int}
    hopping::Tt
    interaction::Tu
    potential::Tw
    seed::Ts
    filling::Float64
    boundary::Symbol
end

"""
    GraphModel(; hopping, interaction, potential=nothing, seed=nothing,
               filling=0.5, probe_codes=nothing)

Describe a real, symmetric, finite tight-binding graph for `method=:kpm`.
`hopping` is a sparse, zero-diagonal matrix whose nonzero off-diagonal entries
define both the hopping graph and the nearest-neighbour interaction graph.
`potential` and `seed` are optional site vectors. The current nested Hadamard
estimator requires a power-of-two number of sites. `probe_codes`, when given,
is a zero-based permutation that defines the nested probe hierarchy; its
default is the physical storage order.

This is intentionally KPM-only. It does not imply that the legacy MPO or
dense compatibility kernels support arbitrary graphs.
"""
struct GraphModel <: AbstractPhysicalModel
    hopping::SparseMatrixCSC{Float64,Int}
    interaction::Float64
    potential::Vector{Float64}
    seed::Vector{Float64}
    filling::Float64
    probe_codes::Vector{Int}
end

function _graph_site_vector(values, N::Int, label::AbstractString)
    isnothing(values) && return zeros(Float64, N)
    values isa AbstractVector || throw(ArgumentError("$label must be nothing or a vector of length $N"))
    length(values) == N || throw(ArgumentError("$label must have length $N, got $(length(values))"))
    result = Float64.(values)
    all(isfinite, result) || throw(ArgumentError("$label must contain only finite real values"))
    return result
end

function GraphModel(
    ; hopping::SparseMatrixCSC,
    interaction,
    potential=nothing,
    seed=nothing,
    filling::Real=0.5,
    probe_codes=nothing,
)
    size(hopping, 1) == size(hopping, 2) || throw(ArgumentError("graph hopping matrix must be square"))
    N = _validated_power_of_two_size(size(hopping, 1), "graph site count")
    eltype(hopping) <: Real || throw(ArgumentError("graph hopping matrix must be real"))
    graph_hopping = SparseMatrixCSC{Float64,Int}(hopping)
    all(iszero, diag(graph_hopping)) || throw(ArgumentError(
        "graph hopping must have a zero diagonal; place onsite terms in potential",
    ))
    isapprox(graph_hopping, graph_hopping'; rtol=0, atol=0) || throw(ArgumentError(
        "graph hopping matrix must be exactly symmetric",
    ))
    codes = isnothing(probe_codes) ? collect(0:(N - 1)) : Int.(probe_codes)
    length(codes) == N && sort(codes) == collect(0:(N - 1)) || throw(ArgumentError(
        "probe_codes must be a zero-based permutation of 0:$(N - 1)",
    ))
    return GraphModel(
        graph_hopping, Float64(_validated_interaction(interaction)),
        _graph_site_vector(potential, N, "potential"), _graph_site_vector(seed, N, "seed"),
        _validated_filling(filling, N), codes,
    )
end

function SquareModel(
    ; size::Tuple{<:Integer,<:Integer},
    hopping,
    interaction,
    potential=nothing,
    seed=nothing,
    filling::Real=0.5,
    boundary::Symbol=:open,
)
    boundary == :open || throw(ArgumentError(
        "the current SquareModel compatibility layer supports only boundary=:open, got $boundary",
    ))
    nx = _validated_power_of_two_size(size[1], "square size along x")
    ny = _validated_power_of_two_size(size[2], "square size along y")
    hopping isa Tuple && length(hopping) == 2 || throw(ArgumentError(
        "square hopping must contain exactly (t_x, t_y)",
    ))
    _validated_hopping_component(hopping[1], "square hopping t_x")
    _validated_hopping_component(hopping[2], "square hopping t_y")
    _validated_interaction(interaction)
    _validated_optional_field(potential, "potential")
    _validated_optional_field(seed, "seed")
    physical_filling = _validated_filling(filling, Base.checked_mul(nx, ny))
    return SquareModel{typeof(hopping),typeof(interaction),typeof(potential),typeof(seed)}(
        (nx, ny), hopping, interaction, potential, seed, physical_filling, boundary,
    )
end

"""
    QTTSettings(; encoding=:auto, tci_tol=1e-10, cutoff=1e-10, maxdim=256)

Representation controls for QTT construction and MPO compression. `:auto`
resolves to the legacy binary chain or interleaved square encoding during
conversion; no alternate encoding is implemented by this compatibility layer.
"""
struct QTTSettings
    encoding::Symbol
    tci_tol::Float64
    cutoff::Float64
    maxdim::Int
end

function QTTSettings(
    ; encoding::Symbol=:auto,
    tci_tol::Real=1e-10,
    cutoff::Real=1e-10,
    maxdim::Integer=256,
)
    encoding in (:auto, :binary, :interleaved) || throw(ArgumentError(
        "encoding must be :auto, :binary, or :interleaved, got $encoding",
    ))
    isfinite(tci_tol) && tci_tol > 0 || throw(ArgumentError(
        "tci_tol must be finite and positive, got $tci_tol",
    ))
    isfinite(cutoff) && cutoff > 0 || throw(ArgumentError(
        "cutoff must be finite and positive, got $cutoff",
    ))
    0 < maxdim <= typemax(Int) || throw(ArgumentError(
        "maxdim must be a positive Int-compatible integer, got $maxdim",
    ))
    return QTTSettings(encoding, Float64(tci_tol), Float64(cutoff), Int(maxdim))
end

"""
    SCFSettings(; purification=:sp2, mixing=0.5, tolerance=0.1,
                maxiter=30, purification_maxiter=50,
                square_fock_method=:binary_carry,
                sp2_idempotency_tolerance=1e-3,
                sp2_relative_trace_tolerance=1e-6,
                record_energy=true, stable_iterations=2, require_stationarity=true,
                measure_stationarity=true,
                detect_two_cycles=true, mixing_method=:linear,
                pulay_history=4, pulay_warmup=4,
                pulay_regularization=1e-12,
                pulay_coefficient_limit=8.0, pulay_step_limit=20.0)

Algorithm controls for the public SCF solvers. The SP2 trace tolerance is
stored relative to the canonical particle count and converted to the absolute
trace tolerance required by the legacy MPO kernel at solve time.
"""
struct SCFSettings
    purification::Symbol
    mixing::Float64
    tolerance::Float64
    maxiter::Int
    purification_maxiter::Int
    square_fock_method::Symbol
    sp2_idempotency_tolerance::Float64
    sp2_relative_trace_tolerance::Float64
    record_energy::Bool
    stable_iterations::Int
    require_stationarity::Bool
    measure_stationarity::Bool
    detect_two_cycles::Bool
    mixing_method::Symbol
    pulay_history::Int
    pulay_warmup::Int
    pulay_regularization::Float64
    pulay_coefficient_limit::Float64
    pulay_step_limit::Float64
end

"""Numerical controls for the public sparse-vector KPM solver.

`moments` controls the Jackson-smoothed Chebyshev projector. `probes` is the
size of a deterministic, nested Hadamard hierarchy used for both the trace and
local density/bond estimators. The audit is an independent, fixed-field KPM
calculation; it is not another SCF solve.
"""
struct KPMSettings
    moments::Int
    probes::Int
    probe_seed::Int
    audit_moments::Int
    audit_probes::Int
    audit_seed::Int
end

function KPMSettings(; moments::Integer=1200, probes::Integer=1024,
    probe_seed::Integer=510578, audit_moments::Integer=1600,
    audit_probes::Integer=1024, audit_seed::Integer=20260730)
    moments >= 2 || throw(ArgumentError("KPM moments must be at least 2"))
    audit_moments >= 2 || throw(ArgumentError("KPM audit_moments must be at least 2"))
    probes > 0 || throw(ArgumentError("KPM probes must be positive"))
    audit_probes > 0 || throw(ArgumentError("KPM audit_probes must be positive"))
    return KPMSettings(Int(moments), Int(probes), Int(probe_seed), Int(audit_moments),
        Int(audit_probes), Int(audit_seed))
end

"""
    RuntimeSettings(; backend=:cpu, device_scalar_type=Float64)

Execution choices separate from the physical model and numerical
representation. The current public MPO CUDA path keeps purification and
effective-Hamiltonian work on the GPU, but moves the density to the host for
Hartree/Fock extraction each SCF iteration. Dense public HF is CPU-only.
"""
struct RuntimeSettings
    backend::Symbol
    device_scalar_type::DataType
end

function RuntimeSettings(; backend::Symbol=:cpu, device_scalar_type::DataType=Float64)
    backend in (:cpu, :cuda) || throw(ArgumentError("backend must be :cpu or :cuda, got $backend"))
    device_scalar_type == Float64 || throw(ArgumentError(
        "the current public MPO CUDA adapter supports only device_scalar_type=Float64, got $device_scalar_type",
    ))
    return RuntimeSettings(backend, device_scalar_type)
end

"""Actual runtime capability and execution-path report from [`runtime_preflight`](@ref)."""
struct RuntimeReport
    requested_backend::Symbol
    active_backend::Union{Nothing,Symbol}
    execution_path::Symbol
    runnable::Bool
    device_name::Union{Nothing,String}
    device_total_memory_bytes::Union{Nothing,Int}
    device_free_memory_before_bytes::Union{Nothing,Int}
    message::String
end

const _CUDA_PACKAGE_ID = Base.PkgId(
    Base.UUID("052768ef-5323-5732-b1bb-66c8b64840ba"), "CUDA",
)

function _cuda_module_or_error()
    try
        # `Base.require` returns the already-loaded module when the CLI has
        # preloaded CUDA. Avoid creating a new global binding with `eval` from
        # inside an active solver call: Julia 1.12 world-age semantics can then
        # hide CUDA extension methods from that call.
        return Base.require(_CUDA_PACKAGE_ID), nothing
    catch error
        return nothing, sprint(showerror, error)
    end
end

_cuda_call(cuda, name::Symbol, arguments...) =
    Base.invokelatest(getproperty(cuda, name), arguments...)

"""
    runtime_preflight(runtime=RuntimeSettings(); method=:mpo)

Check the requested runtime without constructing a `System` or an MPO. CUDA
is loaded lazily only when requested, so ordinary CPU users and login-node
imports remain unaffected.
"""
function runtime_preflight(runtime::RuntimeSettings=RuntimeSettings(); method::Symbol=:mpo)
    method in (:mpo, :dense, :kpm) || throw(ArgumentError("runtime preflight supports :mpo, :dense, or :kpm"))
    if runtime.backend == :cpu
        return RuntimeReport(:cpu, :cpu, :cpu_end_to_end, true, nothing, nothing, nothing,
            "all public solver work executes on the CPU")
    elseif method == :dense
        return RuntimeReport(:cuda, nothing, :unsupported, false, nothing, nothing, nothing,
            "method=:dense is CPU-only; request RuntimeSettings(backend=:cpu)")
    end
    cuda, load_error = _cuda_module_or_error()
    !isnothing(load_error) && return RuntimeReport(:cuda, nothing, :unsupported, false,
        nothing, nothing, nothing, "CUDA import failed: $load_error")
    functional = try
        Bool(_cuda_call(cuda, :functional))
    catch error
        return RuntimeReport(:cuda, nothing, :unsupported, false, nothing, nothing, nothing,
            "CUDA functional check failed: $(sprint(showerror, error))")
    end
    functional || return RuntimeReport(:cuda, nothing, :unsupported, false, nothing, nothing,
        nothing, "CUDA.jl is installed but not functional on this node")
    try
        device = _cuda_call(cuda, :device)
        return RuntimeReport(
            :cuda, :cuda, method == :kpm ? :cuda_kpm_sparse_vectors : :cuda_purification_host_fields, true,
            String(_cuda_call(cuda, :name, device)),
            Int(_cuda_call(cuda, :total_memory)),
            Int(_cuda_call(cuda, :free_memory)),
            method == :kpm ?
            "sparse-vector KPM recurrence on CUDA; filtered probes return to the host for local estimators" :
            "purification/effective Hamiltonian on CUDA; density and mean-field extraction cross the host-device boundary each SCF iteration",
        )
    catch error
        return RuntimeReport(:cuda, nothing, :unsupported, false, nothing, nothing, nothing,
            "CUDA device query failed: $(sprint(showerror, error))")
    end
end

"""
    PreflightIssue(severity, code, message)

A non-mutating diagnostic emitted by [`preflight`](@ref). `severity` is either
`:error` (the requested compatibility path cannot run) or `:warning` (the
configuration is accepted, but requires an explicit choice before solving).
"""
struct PreflightIssue
    severity::Symbol
    code::Symbol
    message::String
end

"""
    PreflightReport

Pure summary of the model geometry, QTT representation, canonical particle
count, requested legacy solver path, and any capability or spectral-bound
issues. Constructing this report never allocates an MPO, initializes CUDA, or
creates an output directory.
"""
struct PreflightReport
    method::Symbol
    model_kind::Symbol
    geometry::Symbol
    physical_size::Tuple{Int,Int}
    physical_sites::Int
    qtt_bit_counts::Tuple{Int,Int}
    qtt_levels::Int
    encoding::Union{Nothing,Symbol}
    target_particles::Int
    representation::QTTSettings
    solver::SCFSettings
    kpm::KPMSettings
    runtime::RuntimeSettings
    runtime_report::RuntimeReport
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}
    spectral_bounds_status::Symbol
    runnable::Bool
    issues::Vector{PreflightIssue}
end

function Base.show(io::IO, report::PreflightReport)
    print(io, "Preflight report\n")
    print(io, "  method: ", report.method, "\n")
    print(io, "  geometry: ", report.geometry)
    report.model_kind == :graph ? print(io, ", explicit edges\n") : print(io, ", open\n")
    print(io, "  physical lattice: ", report.physical_size[1])
    report.model_kind == :square && print(io, " × ", report.physical_size[2])
    print(io, "\n  physical sites: ", report.physical_sites, "\n")
    print(io, "  QTT bits: ", report.qtt_bit_counts, "; levels: ", report.qtt_levels, "\n")
    print(io, "  encoding: ", something(report.encoding, :unresolved), "\n")
    print(io, "  target particles: ", report.target_particles, "\n")
    print(io, "  runtime: requested=", report.runtime.backend,
        ", active=", something(report.runtime_report.active_backend, :unavailable),
        ", path=", report.runtime_report.execution_path, "\n")
    print(io, "  spectral bounds: ", report.spectral_bounds_status)
    !isnothing(report.spectral_bounds) && print(io, " ", report.spectral_bounds)
    print(io, "\n  runnable through current compatibility path: ", report.runnable)
    for issue in report.issues
        print(io, "\n  ", uppercase(String(issue.severity)), " [", issue.code, "]: ", issue.message)
    end
end

_preflight_geometry(model::ChainModel) = (
    model_kind=:chain,
    geometry=:chain,
    physical_size=(model.size, 1),
    physical_sites=model.size,
    bit_counts=(qtt_bits(model), 0),
)

function _preflight_geometry(model::SquareModel)
    nx, ny = model.size
    return (
        model_kind=:square,
        geometry=nx == ny ? :square : :rectangle,
        physical_size=(nx, ny),
        physical_sites=Base.checked_mul(nx, ny),
        bit_counts=qtt_bits(model),
    )
end

_preflight_geometry(model::GraphModel) = (
    model_kind=:graph,
    geometry=:graph,
    physical_size=(size(model.hopping, 1), 1),
    physical_sites=size(model.hopping, 1),
    bit_counts=(qtt_bits(model), 0),
)

function _preflight_encoding(model, representation, issues)
    try
        return _resolved_encoding(model, representation)
    catch error
        error isa ArgumentError || rethrow()
        push!(issues, PreflightIssue(:error, :unsupported_encoding, sprint(showerror, error)))
        return nothing
    end
end

function _preflight_spectral_bounds(
    solver::SCFSettings,
    spectral_bounds,
    issues;
    required::Bool=true,
)
    if isnothing(spectral_bounds)
        if required && solver.purification in (:sp2, :mcweeny_mu)
            push!(issues, PreflightIssue(
                :error,
                :spectral_bounds_required,
                "purification=$(solver.purification) requires enclosing spectral bounds; " *
                "the public bound-policy wrapper is not implemented yet, so supply " *
                "spectral_bounds=(H_min, H_max) to preflight and run_scf!.",
            ))
            return nothing, :required_not_supplied
        end
        return nothing, :not_required_by_selected_purification
    end
    spectral_bounds isa Tuple && length(spectral_bounds) == 2 || begin
        push!(issues, PreflightIssue(
            :error, :invalid_spectral_bounds,
            "spectral_bounds must be a two-tuple (H_min, H_max)",
        ))
        return nothing, :invalid
    end
    try
        bounds = validate_spectral_bounds(spectral_bounds...)
        return bounds, :supplied_unverified
    catch error
        error isa ArgumentError || rethrow()
        push!(issues, PreflightIssue(:error, :invalid_spectral_bounds, sprint(showerror, error)))
        return nothing, :invalid
    end
end

"""
    preflight(model, representation=QTTSettings(), solver=SCFSettings(),
              runtime=RuntimeSettings(); method=:mpo, spectral_bounds=nothing,
              kpm=KPMSettings())

Inspect the public configuration before any MPO, GPU, or output allocation.
`method` presently describes a compatibility path, not a new solver: only
`:mpo`, `:dense`, and a sparse graph-based `:kpm` adapter are available.
The public model accepts rectangles, but reports them as unavailable until the
rectangular MPO/dense kernel plan is completed. The public KPM adapter accepts
open power-of-two `SquareModel` rectangles through a row-major sparse graph;
its probe hierarchy uses interleaved coordinate-bit colors.
"""
function preflight(
    model::AbstractPhysicalModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings();
    runtime::RuntimeSettings=RuntimeSettings(),
    method::Symbol=:mpo,
    spectral_bounds=nothing,
    kpm::KPMSettings=KPMSettings(),
)
    method in (:mpo, :dense, :kpm) || throw(ArgumentError(
        "method must be :mpo, :dense, or :kpm, got $method",
    ))
    geometry = _preflight_geometry(model)
    issues = PreflightIssue[]
    encoding = _preflight_encoding(model, representation, issues)
    runtime_report = runtime_preflight(runtime; method=method)
    !runtime_report.runnable && push!(issues, PreflightIssue(
        :error, :runtime_unavailable, runtime_report.message,
    ))
    bounds, bounds_status = _preflight_spectral_bounds(
        solver, spectral_bounds, issues; required=method == :mpo,
    )

    if method == :dense && !isnothing(spectral_bounds)
        push!(issues, PreflightIssue(
            :error,
            :spectral_bounds_not_used_by_dense,
            "method=:dense diagonalizes the effective Hamiltonian directly; omit spectral_bounds.",
        ))
    end

    if geometry.geometry == :rectangle && method in (:mpo, :dense)
        push!(issues, PreflightIssue(
            :error,
            :rectangular_backend_unavailable,
            "SquareModel(size=$(geometry.physical_size)) is a valid public geometry, " *
            "but the current legacy $method backend assumes equal x/y QTT bit counts. " *
            "See rectangular stages R1--R6 in docs/plans/public_api_redesign.md.",
        ))
    elseif method in (:mpo, :dense) && model isa GraphModel
        push!(issues, PreflightIssue(
            :error, :graph_backend_unavailable,
            "GraphModel is currently supported only by method=:kpm; legacy MPO and dense " *
            "compatibility kernels require ChainModel or equal-square SquareModel.",
        ))
    elseif method == :kpm && !(model isa Union{SquareModel,GraphModel})
        push!(issues, PreflightIssue(
            :error, :kpm_geometry_unavailable,
            "the public KPM adapter currently supports SquareModel and GraphModel; " *
            "ChainModel KPM support has not been added yet.",
        ))
    elseif method == :kpm && (kpm.probes > geometry.physical_sites || kpm.audit_probes > geometry.physical_sites)
        push!(issues, PreflightIssue(
            :error, :kpm_probe_count_invalid,
            "KPM probe counts must not exceed the physical site count $(geometry.physical_sites)",
        ))
    end

    runnable = isempty(filter(issue -> issue.severity == :error, issues))
    return PreflightReport(
        method, geometry.model_kind, geometry.geometry, geometry.physical_size,
        geometry.physical_sites, geometry.bit_counts, qtt_levels(model), encoding,
        round(Int, geometry.physical_sites * model.filling), representation, solver, kpm, runtime,
        runtime_report,
        bounds, bounds_status, runnable, issues,
    )
end

"""
    SolveResult

Normalized public result from the current compatibility solver. It stores
physical model/settings provenance, scalar SCF history, and measured local
observables, but deliberately does not expose mutable `System` state or raw
MPO tensors. `observables` is `nothing` only when purification fails before a
candidate density can be installed in the legacy system.
"""
struct SolveResult{M}
    model::M
    representation::QTTSettings
    solver::SCFSettings
    kpm::Union{Nothing,KPMSettings}
    runtime::RuntimeSettings
    runtime_report::RuntimeReport
    method::Symbol
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}
    target_particles::Int
    converged::Bool
    termination_reason::Symbol
    elapsed_time_s::Float64
    diagnostics::Any
    observables::Any
end

function Base.show(io::IO, result::SolveResult)
    print(io, "SolveResult(method=", result.method)
    print(io, ", backend=", something(result.runtime_report.active_backend, :unavailable))
    print(io, ", converged=", result.converged)
    print(io, ", termination=", result.termination_reason)
    print(io, ", elapsed_s=", result.elapsed_time_s)
    print(io, ", particles=", result.target_particles)
    print(io, ", spectral_bounds=", something(result.spectral_bounds, :not_used), ")")
end

function _preflight_error_message(report::PreflightReport)
    errors = filter(issue -> issue.severity == :error, report.issues)
    return join(("[$(issue.code)] $(issue.message)" for issue in errors), "\n")
end

function _mpo_runtime_callbacks(report::RuntimeReport)
    report.active_backend == :cpu && return (identity, identity, () -> nothing)
    report.active_backend == :cuda || throw(ArgumentError(
        "MPO runtime is not active: $(report.message)",
    ))
    cuda, load_error = _cuda_module_or_error()
    isnothing(load_error) || throw(ArgumentError("CUDA import failed after preflight: $load_error"))
    to_device = value -> ITensors.adapt(getproperty(cuda, :CuArray), value)
    to_host = ITensors.cpu
    synchronize = () -> _cuda_call(cuda, :synchronize)
    return to_device, to_host, synchronize
end

"""
    solve(model, representation=QTTSettings(), solver=SCFSettings();
          method=:mpo, spectral_bounds, verbose=:nothing, kwargs...)

Run the existing MPO Hartree--Fock implementation through the public
configuration types. This is a compatibility wrapper: it constructs legacy
parameters and `System`, calls the existing `run_scf!`, then measures existing
geometry-appropriate observables. It does not alter equations, QTT ordering,
spectral scaling, truncation, or stopping criteria.

`method=:mpo` is the legacy purification path and requires
`spectral_bounds=(H_min,H_max)` whenever its purification method needs them.
`method=:dense` is an independent canonical, real dense Hartree--Fock
reference and does not use purification or spectral bounds. Remaining MPO
keywords are forwarded unchanged to `run_scf!`, except `purification_method`,
which is owned by `SCFSettings.purification`.
"""
function solve(
    model::AbstractPhysicalModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings();
    runtime::RuntimeSettings=RuntimeSettings(),
    method::Symbol=:mpo,
    kpm::KPMSettings=KPMSettings(),
    spectral_bounds=nothing,
    verbose::Symbol=:nothing,
    checkpoint_path::Union{Nothing,AbstractString}=nothing,
    measure_observables::Bool=true,
    kwargs...,
)
    method in (:mpo, :dense, :kpm) || throw(ArgumentError(
        "public solve currently supports method=:mpo, :dense, or :kpm; requested method=$method",
    ))
    method != :mpo && !isnothing(checkpoint_path) && throw(ArgumentError(
        "checkpoint_path is currently supported only by method=:mpo",
    ))
    method != :mpo && !measure_observables && throw(ArgumentError(
        "measure_observables=false is currently supported only by method=:mpo",
    ))
    owned_solver_keywords = (
        :purification_method, :square_fock_method, :sp2_idempotency_tolerance,
        :sp2_trace_tolerance, :record_energy, :stable_iterations, :require_stationarity,
        :measure_stationarity, :detect_two_cycles,
        :mixing_method, :pulay_history, :pulay_warmup, :pulay_regularization,
        :pulay_coefficient_limit, :pulay_step_limit,
    )
    conflicting = filter(keyword -> haskey(kwargs, keyword), owned_solver_keywords)
    isempty(conflicting) || throw(ArgumentError(
        "$(join(string.(conflicting), ", ")) is controlled by SCFSettings; do not override it in solve",
    ))
    haskey(kwargs, :to_gpu) && throw(ArgumentError(
        "to_gpu is controlled by RuntimeSettings; pass runtime=RuntimeSettings(backend=:cuda) instead",
    ))
    haskey(kwargs, :to_cpu) && throw(ArgumentError(
        "to_cpu is controlled by RuntimeSettings; host/device transfers are selected by the runtime",
    ))

    report = preflight(
        model, representation, solver;
        runtime=runtime, method=method, spectral_bounds=spectral_bounds, kpm=kpm,
    )
    report.runnable || throw(ArgumentError(
        "public solve cannot run this configuration:\n" * _preflight_error_message(report),
    ))
    started = time_ns()
    if method == :kpm
        isempty(kwargs) || throw(ArgumentError(
            "method=:kpm accepts controls through KPMSettings and SCFSettings; got $(collect(keys(kwargs)))",
        ))
        model isa Union{SquareModel,GraphModel} || throw(ArgumentError(
            "public KPM currently requires SquareModel or GraphModel",
        ))
        kpm_result = _solve_kpm(model, solver, kpm, report.runtime_report)
        return SolveResult(
            model, representation, solver, kpm, runtime, report.runtime_report, method, nothing,
            report.target_particles, kpm_result.converged, kpm_result.termination_reason,
            (time_ns() - started) / 1e9, kpm_result.diagnostics, kpm_result.observables,
        )
    end
    params = legacy_parameters(model, representation, solver)
    if method == :dense
        isempty(kwargs) || throw(ArgumentError(
            "method=:dense accepts no MPO-specific run_scf! keywords in this stage; got $(collect(keys(kwargs)))",
        ))
        dense = _solve_dense_hf(params)
        return SolveResult(
            model, representation, solver, nothing, runtime, report.runtime_report, method, nothing, report.target_particles,
            dense.converged, dense.termination_reason, (time_ns() - started) / 1e9,
            dense.diagnostics, dense.observables,
        )
    end
    bounds = isnothing(report.spectral_bounds) ? throw(ArgumentError(
        "public solve requires explicit spectral_bounds for this purification method",
    )) : report.spectral_bounds
    to_device, to_host, synchronize = _mpo_runtime_callbacks(report.runtime_report)
    static_phase_callback = verbose == :nothing ? nothing : event -> begin
        @printf(stdout, "Initialization %-32s %.3f s\n", string(event.phase), event.elapsed_time_s)
        flush(stdout)
    end
    verbose != :nothing && begin
        println("Initializing static QTT/MPO operators on $(report.runtime_report.active_backend)...")
        flush(stdout)
    end
    sys = System(
        params;
        static_to_backend=to_device,
        static_phase_callback=static_phase_callback,
        static_phase_synchronize=synchronize,
    )
    converged = run_scf!(
        sys, bounds...;
        verbose=verbose,
        purification_method=solver.purification,
        square_fock_method=solver.square_fock_method,
        sp2_idempotency_tolerance=solver.sp2_idempotency_tolerance,
        sp2_trace_tolerance=solver.sp2_relative_trace_tolerance * report.target_particles,
        record_energy=solver.record_energy,
        stable_iterations=solver.stable_iterations,
        require_stationarity=solver.require_stationarity,
        measure_stationarity=solver.measure_stationarity,
        detect_two_cycles=solver.detect_two_cycles,
        mixing_method=solver.mixing_method,
        pulay_history=solver.pulay_history,
        pulay_warmup=solver.pulay_warmup,
        pulay_regularization=solver.pulay_regularization,
        pulay_coefficient_limit=solver.pulay_coefficient_limit,
        pulay_step_limit=solver.pulay_step_limit,
        to_gpu=to_device,
        to_cpu=to_host,
        phase_synchronize=synchronize,
        kwargs...,
    )
    if report.runtime_report.active_backend == :cuda
        synchronize()
        sys.H0 = to_host(sys.H0)
        sys.VH = to_host(sys.VH)
        sys.VF = to_host(sys.VF)
        synchronize()
    end
    diagnostics = scf_diagnostics(sys)
    if converged && !isnothing(checkpoint_path)
        write_mpo_checkpoint(sys, checkpoint_path)
    end
    measured_observables = (
        measure_observables && diagnostics.termination_reason != :purification_failed
    ) ? observables(sys; measure_stationarity=solver.measure_stationarity) : nothing
    return SolveResult(
        model, representation, solver, nothing, runtime, report.runtime_report, method, bounds, report.target_particles,
        converged, diagnostics.termination_reason, (time_ns() - started) / 1e9,
        diagnostics, measured_observables,
    )
end

_dense_relative_change(candidate, reference) =
    norm(candidate - reference) / max(norm(reference), sqrt(eps(Float64)))

function _dense_chain_terms(params::Parameters1D)
    N = 2^params.L
    onsite = isnothing(params.W) ? zeros(Float64, N) :
        [Float64(params.W(site - 1)) for site in 1:N]
    hopping = [Float64(params.t isa Number ? params.t : params.t(bond - 1)) for bond in 1:(N - 1)]
    seed = isnothing(params.S) ? zeros(Float64, N) :
        [Float64(params.S(site - 1)) for site in 1:N]
    return onsite, hopping, seed
end

function _dense_chain_mean_fields(rho::Matrix{Float64}, U::Float64)
    N = size(rho, 1)
    hartree = zeros(Float64, N)
    for site in 1:N
        site > 1 && (hartree[site] += U * rho[site - 1, site - 1])
        site < N && (hartree[site] += U * rho[site + 1, site + 1])
    end
    return hartree, [-U * rho[site, site + 1] for site in 1:(N - 1)]
end

function _dense_chain_energy(onsite, hopping, rho::Matrix{Float64}, U::Float64)
    density = diag(rho)
    bonds = [rho[site, site + 1] for site in 1:(size(rho, 1) - 1)]
    kinetic = sum(onsite .* density) + 2sum(hopping .* bonds)
    hartree = U * sum(density[site] * density[site + 1] for site in 1:(length(density) - 1))
    fock = -U * sum(abs2, bonds)
    return (; kinetic, hartree, fock, interaction=hartree + fock, total=kinetic + hartree + fock)
end

function _dense_square_h0(params::ParametersSquare)
    L = params.L
    side = 2^div(L, 2)
    H0 = zeros(Float64, side^2, side^2)
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) : Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) : Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, L)
        H0[site, site] = isnothing(params.W) ? 0.0 : Float64(params.W(x, y))
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, L)
            H0[site, neighbour] = H0[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, L)
            H0[site, neighbour] = H0[neighbour, site] = ty(x, y)
        end
    end
    return H0
end

function _dense_square_seed(params::ParametersSquare)
    seed = zeros(Float64, 2^params.L, 2^params.L)
    isnothing(params.S) && return seed
    for site in axes(seed, 1)
        x, y = square_lattice_decoder(site - 1, params.L)
        seed[site, site] = Float64(params.S(x, y))
    end
    return seed
end

function _dense_square_mean_fields(rho::Matrix{Float64}, U::Float64, L::Int)
    N = size(rho, 1)
    hartree = zeros(Float64, N, N)
    fock = zeros(Float64, N, N)
    for site in 1:N
        for neighbour in values(square_neighbours(site, L))
            !isnothing(neighbour) && (hartree[site, site] += U * rho[neighbour, neighbour])
        end
    end
    for (site, neighbour, _) in square_undirected_bonds(L)
        fock[site, neighbour] = fock[neighbour, site] = -U * rho[site, neighbour]
    end
    return hartree, fock
end

function _dense_square_energy(H0, rho::Matrix{Float64}, U::Float64, L::Int)
    density = diag(rho)
    hartree = 0.0
    fock = 0.0
    for (site, neighbour, _) in square_undirected_bonds(L)
        hartree += U * density[site] * density[neighbour]
        fock -= U * abs2(rho[site, neighbour])
    end
    kinetic = real(tr(H0 * rho))
    return (; kinetic, hartree, fock, interaction=hartree + fock, total=kinetic + hartree + fock)
end

function _dense_diagnostics_record(
    iteration,
    rho,
    new_hartree,
    old_hartree,
    new_fock,
    old_fock,
    rho_previous,
    rho_two_steps_ago,
    H,
    energy_total,
)
    vh_residual = iteration == 1 ? Inf : _dense_relative_change(new_hartree, old_hartree)
    vf_residual = iteration == 1 ? Inf : _dense_relative_change(new_fock, old_fock)
    rho_residual = iteration == 1 ? Inf : _dense_relative_change(rho, rho_previous)
    two_cycle_residual = iteration < 3 ? Inf : _dense_relative_change(rho, rho_two_steps_ago)
    commutator = _dense_relative_change(H * rho, rho * H)
    return SCFIterationRecord(
        iteration, real(tr(rho)), vh_residual, vf_residual, rho_residual, commutator,
        two_cycle_residual, true, :not_applicable, 0, 0, nothing, nothing, nothing,
        nothing, energy_total,
    )
end

function _dense_observables_chain(onsite, hopping, rho, H, U)
    energy = _dense_chain_energy(onsite, hopping, rho, U)
    return (
        site_density=Vector{Float64}(diag(rho)),
        bond_order=ComplexF64[rho[site, site + 1] for site in 1:(size(rho, 1) - 1)],
        particle_number=real(tr(rho)), energy,
        hermiticity_residual=_dense_relative_change(rho, rho'),
        idempotency_residual=_dense_relative_change(rho * rho, rho),
        stationarity_residual=_dense_relative_change(H * rho, rho * H),
    )
end

function _dense_observables_square(H0, rho, H, U, L)
    horizontal_bonds = Tuple{Int,Int}[]
    vertical_bonds = Tuple{Int,Int}[]
    horizontal_bond_order = ComplexF64[]
    vertical_bond_order = ComplexF64[]
    for (site, neighbour, orientation) in square_undirected_bonds(L)
        if orientation == :horizontal
            push!(horizontal_bonds, (site, neighbour))
            push!(horizontal_bond_order, rho[site, neighbour])
        else
            push!(vertical_bonds, (site, neighbour))
            push!(vertical_bond_order, rho[site, neighbour])
        end
    end
    energy = _dense_square_energy(H0, rho, U, L)
    return (
        site_density=Vector{Float64}(diag(rho)), horizontal_bonds, vertical_bonds,
        horizontal_bond_order, vertical_bond_order, particle_number=real(tr(rho)), energy,
        hermiticity_residual=_dense_relative_change(rho, rho'),
        idempotency_residual=_dense_relative_change(rho * rho, rho),
        stationarity_residual=_dense_relative_change(H * rho, rho * H),
    )
end

function _solve_dense_hf(params::Parameters1D)
    N = 2^params.L
    Ne = round(Int, params.density * N)
    onsite, hopping, hartree = _dense_chain_terms(params)
    fock = zeros(Float64, N - 1)
    rho_previous = zeros(Float64, N, N)
    rho_two_steps_ago = nothing
    history = SCFIterationRecord[]
    tolerance = params.scf_tol / 100
    for iteration in 1:params.scf_max_iterations
        H = Matrix(SymTridiagonal(onsite .+ hartree, hopping .+ fock))
        eigenpairs = eigen(Symmetric(H))
        rho = Matrix(@view(eigenpairs.vectors[:, 1:Ne]) * @view(eigenpairs.vectors[:, 1:Ne])')
        new_hartree, new_fock = _dense_chain_mean_fields(rho, Float64(params.U))
        energy = _dense_chain_energy(onsite, hopping, rho, Float64(params.U))
        record = _dense_diagnostics_record(
            iteration, rho, new_hartree, hartree, new_fock, fock, rho_previous,
            rho_two_steps_ago, H, energy.total,
        )
        push!(history, record)
        stable = iteration >= 2 && _scf_record_within_tolerance(record, tolerance)
        previous_stable = length(history) >= 2 && _scf_record_within_tolerance(history[end - 1], tolerance)
        diagnostics = if stable && previous_stable
            SCFDiagnostics(history, true, :converged)
        elseif _scf_is_two_cycle(record, tolerance)
            SCFDiagnostics(history, false, :two_cycle_detected)
        else
            nothing
        end
        if !isnothing(diagnostics)
            return (; converged=diagnostics.converged, termination_reason=diagnostics.termination_reason,
                diagnostics, observables=_dense_observables_chain(onsite, hopping, rho, H, Float64(params.U)))
        end
        rho_two_steps_ago = rho_previous
        if iteration == 1
            hartree, fock = new_hartree, new_fock
        else
            hartree = params.scf_mixing .* new_hartree .+ (1 - params.scf_mixing) .* hartree
            fock = params.scf_mixing .* new_fock .+ (1 - params.scf_mixing) .* fock
        end
        rho_previous = rho
    end
    H = Matrix(SymTridiagonal(onsite .+ hartree, hopping .+ fock))
    return (; converged=false, termination_reason=:max_iterations,
        diagnostics=SCFDiagnostics(history, false, :max_iterations),
        observables=_dense_observables_chain(onsite, hopping, rho_previous, H, Float64(params.U)))
end

function _solve_dense_hf(params::ParametersSquare)
    N = 2^params.L
    Ne = round(Int, params.density * N)
    H0 = _dense_square_h0(params)
    hartree = _dense_square_seed(params)
    fock = zeros(Float64, N, N)
    rho_previous = zeros(Float64, N, N)
    rho_two_steps_ago = nothing
    history = SCFIterationRecord[]
    tolerance = params.scf_tol / 100
    for iteration in 1:params.scf_max_iterations
        H = H0 + hartree + fock
        eigenpairs = eigen(Symmetric(H))
        rho = Matrix(@view(eigenpairs.vectors[:, 1:Ne]) * @view(eigenpairs.vectors[:, 1:Ne])')
        new_hartree, new_fock = _dense_square_mean_fields(rho, Float64(params.U), params.L)
        energy = _dense_square_energy(H0, rho, Float64(params.U), params.L)
        record = _dense_diagnostics_record(
            iteration, rho, new_hartree, hartree, new_fock, fock, rho_previous,
            rho_two_steps_ago, H, energy.total,
        )
        push!(history, record)
        stable = iteration >= 2 && _scf_record_within_tolerance(record, tolerance)
        previous_stable = length(history) >= 2 && _scf_record_within_tolerance(history[end - 1], tolerance)
        diagnostics = if stable && previous_stable
            SCFDiagnostics(history, true, :converged)
        elseif _scf_is_two_cycle(record, tolerance)
            SCFDiagnostics(history, false, :two_cycle_detected)
        else
            nothing
        end
        if !isnothing(diagnostics)
            return (; converged=diagnostics.converged, termination_reason=diagnostics.termination_reason,
                diagnostics, observables=_dense_observables_square(H0, rho, H, Float64(params.U), params.L))
        end
        rho_two_steps_ago = rho_previous
        if iteration == 1
            hartree, fock = new_hartree, new_fock
        else
            hartree = params.scf_mixing * new_hartree + (1 - params.scf_mixing) * hartree
            fock = params.scf_mixing * new_fock + (1 - params.scf_mixing) * fock
        end
        rho_previous = rho
    end
    H = H0 + hartree + fock
    return (; converged=false, termination_reason=:max_iterations,
        diagnostics=SCFDiagnostics(history, false, :max_iterations),
        observables=_dense_observables_square(H0, rho_previous, H, Float64(params.U), params.L))
end

const _PUBLIC_RESULT_FORMAT_VERSION = 1

function _toml_escape(value::AbstractString)
    return replace(value, '\\' => "\\\\", '"' => "\\\"", '\n' => "\\n")
end

function _toml_value(value)
    if value isa AbstractString
        return '"' * _toml_escape(value) * '"'
    elseif value isa Bool || value isa Integer || value isa AbstractFloat
        return repr(value)
    elseif value isa AbstractVector
        return "[" * join(_toml_value.(value), ", ") * "]"
    end
    throw(ArgumentError("unsupported public-result TOML value type: $(typeof(value))"))
end

function _write_public_toml(path::AbstractString, values::AbstractDict)
    open(path, "w") do io
        for key in sort!(collect(keys(values)))
            println(io, key, " = ", _toml_value(values[key]))
        end
    end
end

function _split_toml_vector(text::AbstractString)
    parts = String[]
    current = IOBuffer()
    quoted = false
    escaped = false
    for character in text
        if escaped
            write(current, character)
            escaped = false
        elseif character == '\\' && quoted
            write(current, character)
            escaped = true
        elseif character == '"'
            write(current, character)
            quoted = !quoted
        elseif character == ',' && !quoted
            push!(parts, strip(String(take!(current))))
        else
            write(current, character)
        end
    end
    quoted && throw(ArgumentError("unterminated TOML string in vector: $text"))
    tail = strip(String(take!(current)))
    isempty(tail) || push!(parts, tail)
    return parts
end

function _parse_public_toml_value(text::AbstractString)
    value = strip(text)
    if startswith(value, '"') && endswith(value, '"')
        interior = value[nextind(value, firstindex(value)):prevind(value, lastindex(value))]
        return replace(replace(replace(interior, "\\n" => "\n"), "\\\"" => "\""), "\\\\" => "\\")
    elseif value == "true" || value == "false"
        return value == "true"
    elseif startswith(value, '[') && endswith(value, ']')
        interior = strip(value[nextind(value, firstindex(value)):prevind(value, lastindex(value))])
        return isempty(interior) ? Any[] : _parse_public_toml_value.(_split_toml_vector(interior))
    end
    integer = tryparse(Int, value)
    !isnothing(integer) && return integer
    floating = tryparse(Float64, value)
    !isnothing(floating) && return floating
    throw(ArgumentError("unsupported public-result TOML value: $value"))
end

function _read_public_toml(path::AbstractString)
    values = Dict{String,Any}()
    for (line_number, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        isempty(line) && continue
        startswith(line, '#') && continue
        key_value = split(line, '='; limit=2)
        length(key_value) == 2 || throw(ArgumentError("invalid TOML line $line_number in $path"))
        key = strip(key_value[1])
        isempty(key) && throw(ArgumentError("empty TOML key at line $line_number in $path"))
        values[key] = _parse_public_toml_value(key_value[2])
    end
    return values
end

"""Portable, parsed contents of a directory written by [`write_result`](@ref)."""
struct StoredResult
    directory::String
    input::Dict{String,Any}
    metadata::Dict{String,Any}
    observables::Dict{String,Any}
    history::Vector{Dict{String,String}}
    site_density::Vector{Float64}
    bond_order::Vector{NamedTuple{(:site_left, :site_right, :orientation, :value),Tuple{Int,Int,Symbol,ComplexF64}}}
end

function _csv_escape(value)
    return '"' * replace(string(value), '"' => "\"\"") * '"'
end

_write_csv_row(io, values) = println(io, join(_csv_escape.(values), ','))

function _model_input_summary(model::ChainModel)
    return Dict(
        "model_kind" => "chain",
        "geometry" => "chain",
        "size" => [model.size],
        "hopping" => string(model.hopping),
        "hopping_representation" => model.hopping isa Function ? "function" : "scalar",
        "interaction" => Float64(model.interaction),
        "potential_representation" => isnothing(model.potential) ? "none" : "function",
        "seed_representation" => isnothing(model.seed) ? "none" : "function",
        "filling" => model.filling,
        "boundary" => string(model.boundary),
    )
end

function _model_input_summary(model::SquareModel)
    hopping_representation(component) = component isa Function ? "function" : "scalar"
    return Dict(
        "model_kind" => "square",
        "geometry" => model.size[1] == model.size[2] ? "square" : "rectangle",
        "size" => collect(model.size),
        "hopping_x" => string(model.hopping[1]),
        "hopping_y" => string(model.hopping[2]),
        "hopping_x_representation" => hopping_representation(model.hopping[1]),
        "hopping_y_representation" => hopping_representation(model.hopping[2]),
        "interaction" => Float64(model.interaction),
        "potential_representation" => isnothing(model.potential) ? "none" : "function",
        "seed_representation" => isnothing(model.seed) ? "none" : "function",
        "filling" => model.filling,
        "boundary" => string(model.boundary),
    )
end

function _model_input_summary(model::GraphModel)
    return Dict(
        "model_kind" => "graph",
        "geometry" => "graph",
        "size" => [size(model.hopping, 1)],
        "edge_count" => nnz(model.hopping) ÷ 2,
        "hopping_representation" => "sparse_symmetric_matrix",
        "interaction" => model.interaction,
        "potential_representation" => "site_vector",
        "seed_representation" => "site_vector",
        "probe_codes_representation" => "zero_based_permutation",
        "filling" => model.filling,
        "boundary" => "explicit_graph",
    )
end

function _result_input(result::SolveResult)
    input = _model_input_summary(result.model)
    input["format_version"] = _PUBLIC_RESULT_FORMAT_VERSION
    input["input_serialization"] = "descriptive_not_executable"
    input["function_serialization"] = "functions are described by type/string only; preserve source separately for reproduction"
    input["qtt_encoding"] = string(result.representation.encoding)
    input["qtt_tci_tol"] = result.representation.tci_tol
    input["mpo_cutoff"] = result.representation.cutoff
    input["mpo_maxdim"] = result.representation.maxdim
    input["runtime_requested_backend"] = string(result.runtime.backend)
    input["runtime_device_scalar_type"] = string(result.runtime.device_scalar_type)
    input["purification"] = string(result.solver.purification)
    input["purification_applied"] = result.method == :mpo
    if !isnothing(result.kpm)
        input["kpm_moments"] = result.kpm.moments
        input["kpm_probes"] = result.kpm.probes
        input["kpm_probe_seed"] = result.kpm.probe_seed
        input["kpm_audit_moments"] = result.kpm.audit_moments
        input["kpm_audit_probes"] = result.kpm.audit_probes
        input["kpm_audit_seed"] = result.kpm.audit_seed
        input["kpm_probe_method"] = "hadamard_nested_random_signs"
        input["kpm_probe_ordering"] = result.model isa SquareModel ?
            "interleaved_coordinate_bits" : "explicit_probe_codes"
        input["kpm_site_storage"] = result.model isa SquareModel ?
            "row_major_x_fastest" : "graph_vertex_order"
    end
    input["scf_mixing"] = result.solver.mixing
    input["scf_tolerance_percent"] = result.solver.tolerance
    input["scf_max_iterations"] = result.solver.maxiter
    input["purification_max_iterations"] = result.solver.purification_maxiter
    input["square_fock_method"] = string(result.solver.square_fock_method)
    input["sp2_idempotency_tolerance"] = result.solver.sp2_idempotency_tolerance
    input["sp2_relative_trace_tolerance"] = result.solver.sp2_relative_trace_tolerance
    input["sp2_absolute_trace_tolerance"] =
        result.solver.sp2_relative_trace_tolerance * result.target_particles
    input["record_energy"] = result.solver.record_energy
    input["stable_iterations"] = result.solver.stable_iterations
    input["require_stationarity"] = result.solver.require_stationarity
    input["measure_stationarity"] = result.solver.measure_stationarity
    input["detect_two_cycles"] = result.solver.detect_two_cycles
    input["mixing_method"] = string(result.solver.mixing_method)
    input["pulay_history"] = result.solver.pulay_history
    input["pulay_warmup"] = result.solver.pulay_warmup
    input["pulay_regularization"] = result.solver.pulay_regularization
    input["pulay_coefficient_limit"] = result.solver.pulay_coefficient_limit
    input["pulay_step_limit"] = result.solver.pulay_step_limit
    input["spectral_bounds_status"] = isnothing(result.spectral_bounds) ? "not_used" : "supplied"
    !isnothing(result.spectral_bounds) && (input["spectral_bounds"] = collect(result.spectral_bounds))
    return input
end

function _result_metadata(result::SolveResult)
    bits = qtt_bits(result.model)
    physical_size = result.model isa ChainModel ? (result.model.size, 1) :
        result.model isa SquareModel ? result.model.size : (Base.size(result.model.hopping, 1), 1)
    metadata = Dict(
        "format_version" => _PUBLIC_RESULT_FORMAT_VERSION,
        "created_at_unix_s" => time(),
        "solver" => result.method == :mpo ? "public_mpo_compatibility" :
            result.method == :dense ? "public_dense_hf_compatibility" : "public_sparse_vector_kpm_open_graph_hf",
        "method" => string(result.method),
        "runtime_requested_backend" => string(result.runtime_report.requested_backend),
        "runtime_active_backend" => string(something(result.runtime_report.active_backend, :unavailable)),
        "runtime_execution_path" => string(result.runtime_report.execution_path),
        "runtime_runnable" => result.runtime_report.runnable,
        "runtime_message" => result.runtime_report.message,
        "model_kind" => result.model isa ChainModel ? "chain" :
            result.model isa SquareModel ? "square" : "graph",
        "geometry" => result.model isa ChainModel ? "chain" :
            result.model isa SquareModel ? (physical_size[1] == physical_size[2] ? "square" : "rectangle") : "graph",
        "physical_size" => collect(physical_size),
        "matrix_dimension" => physical_size[1] * physical_size[2],
        "qtt_bit_counts" => result.model isa SquareModel ? collect(bits) : [bits],
        "qtt_levels" => qtt_levels(result.model),
        "target_particles" => result.target_particles,
        "scf_converged" => result.converged,
        "scf_termination_reason" => string(result.termination_reason),
        "scf_iterations" => length(result.diagnostics.history),
        "solve_elapsed_time_s" => result.elapsed_time_s,
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
    )
    if isnothing(result.spectral_bounds)
        metadata["spectral_bounds_status"] = "not_used"
    else
        metadata["spectral_bounds_status"] = "supplied"
        metadata["spectral_lower"] = result.spectral_bounds[1]
        metadata["spectral_upper"] = result.spectral_bounds[2]
    end
    if !isnothing(result.runtime_report.device_name)
        metadata["device_name"] = result.runtime_report.device_name
        metadata["device_total_memory_bytes"] = result.runtime_report.device_total_memory_bytes
        metadata["device_free_memory_before_bytes"] = result.runtime_report.device_free_memory_before_bytes
    end
    if result.diagnostics isa KPMDiagnostics
        audit = result.diagnostics.audit
        metadata["kpm_probe_ordering"] = result.model isa SquareModel ?
            "interleaved_coordinate_bits" : "explicit_probe_codes"
        metadata["kpm_site_storage"] = result.model isa SquareModel ?
            "row_major_x_fastest" : "graph_vertex_order"
        metadata["kpm_hermiticity_residual_status"] = "not_estimated_from_local_kpm_observables"
        metadata["kpm_stationarity_residual_status"] = "not_estimated_from_local_kpm_observables"
        metadata["kpm_audit_trace_error"] = audit.trace_error
        metadata["kpm_audit_trace_idempotency_defect"] = audit.trace_idempotency_defect
        metadata["kpm_audit_density_rms_difference"] = audit.density_rms_difference
        metadata["kpm_audit_bond_rms_difference"] = audit.bond_rms_difference
    end
    return metadata
end

function _result_observables_summary(result::SolveResult)
    isnothing(result.observables) && return Dict("available" => false)
    observables = result.observables
    return Dict(
        "available" => true,
        "particle_number" => observables.particle_number,
        "energy_kinetic" => observables.energy.kinetic,
        "energy_hartree" => observables.energy.hartree,
        "energy_fock" => observables.energy.fock,
        "energy_interaction" => observables.energy.interaction,
        "energy_total" => observables.energy.total,
        "hermiticity_residual" => observables.hermiticity_residual,
        "idempotency_residual" => observables.idempotency_residual,
        "stationarity_residual" => observables.stationarity_residual,
    )
end

function _write_result_history(path::AbstractString, diagnostics::SCFDiagnostics)
    open(path, "w") do io
        _write_csv_row(io, (
            "iteration", "trace", "vh_residual", "vf_residual", "rho_residual",
            "commutator_residual", "two_cycle_residual", "purification_converged",
            "purification_termination_reason", "purification_iterations",
            "purification_selected_iteration", "rho_bond_dimension",
            "hartree_bond_dimension", "fock_bond_dimension",
            "effective_hamiltonian_bond_dimension", "energy_total",
        ))
        for record in diagnostics.history
            _write_csv_row(io, (
                record.iteration, record.trace, record.vh_residual, record.vf_residual,
                record.rho_residual, record.commutator_residual, record.two_cycle_residual,
                record.purification_converged, record.purification_termination_reason,
                record.purification_iterations, record.purification_selected_iteration,
                record.rho_bond_dimension, record.hartree_bond_dimension,
                record.fock_bond_dimension, record.effective_hamiltonian_bond_dimension,
                record.energy_total,
            ))
        end
    end
end

function _write_result_history(path::AbstractString, diagnostics)
    open(path, "w") do io
        _write_csv_row(io, (
            "iteration", "spectral_lower", "spectral_upper", "chemical_potential",
            "trace", "trace_error", "trace_squared_estimate", "trace_idempotency_defect",
            "vh_residual", "vf_residual", "density_residual", "bond_residual",
            "two_cycle_residual", "energy_total", "checkerboard_order",
        ))
        for record in diagnostics.history
            _write_csv_row(io, (
                record.iteration, record.spectral_lower, record.spectral_upper,
                record.chemical_potential, record.trace, record.trace_error,
                record.trace_squared_estimate, record.trace_idempotency_defect,
                record.vh_residual, record.vf_residual, record.density_residual,
                record.bond_residual, record.two_cycle_residual, record.energy_total,
                record.checkerboard_order,
            ))
        end
    end
end

function _write_result_local_observables(directory::AbstractString, result::SolveResult)
    open(joinpath(directory, "site_density.csv"), "w") do io
        _write_csv_row(io, ("site", "density"))
        !isnothing(result.observables) || return
        for (site, density) in enumerate(result.observables.site_density)
            _write_csv_row(io, (site, density))
        end
    end
    open(joinpath(directory, "bond_order.csv"), "w") do io
        _write_csv_row(io, ("site_left", "site_right", "orientation", "real", "imag"))
        !isnothing(result.observables) || return
        observables = result.observables
        if result.model isa ChainModel
            for (site, value) in enumerate(observables.bond_order)
                _write_csv_row(io, (site, site + 1, "chain", real(value), imag(value)))
            end
        elseif result.method == :kpm
            for (bond, value) in zip(observables.bonds, observables.bond_order)
                _write_csv_row(io, (bond[1], bond[2], string(bond[3]), real(value), imag(value)))
            end
        else
            for (bond, value) in zip(observables.horizontal_bonds, observables.horizontal_bond_order)
                _write_csv_row(io, (bond[1], bond[2], "horizontal", real(value), imag(value)))
            end
            for (bond, value) in zip(observables.vertical_bonds, observables.vertical_bond_order)
                _write_csv_row(io, (bond[1], bond[2], "vertical", real(value), imag(value)))
            end
        end
    end
end

"""
    write_result(result, directory)

Write a non-overwriting, format-versioned result directory for a public
[`SolveResult`](@ref). Existing campaign writers are intentionally untouched.
The writer uses the established filenames `metadata.toml`, `observables.toml`,
`scf_history.csv`, `site_density.csv`, and `bond_order.csv`, adding
descriptive `input.toml` rather than attempting to serialize arbitrary Julia
functions as executable source.
"""
function write_result(result::SolveResult, directory::AbstractString)
    output = abspath(directory)
    ispath(output) && throw(ArgumentError("refusing to overwrite existing result directory: $output"))
    mkpath(output)
    _write_public_toml(joinpath(output, "input.toml"), _result_input(result))
    _write_public_toml(joinpath(output, "metadata.toml"), _result_metadata(result))
    _write_public_toml(joinpath(output, "observables.toml"), _result_observables_summary(result))
    _write_result_history(joinpath(output, "scf_history.csv"), result.diagnostics)
    _write_result_local_observables(output, result)
    return output
end

function _parse_csv_line(line::AbstractString)
    fields = String[]
    current = IOBuffer()
    quoted = false
    index = firstindex(line)
    while index <= lastindex(line)
        character = line[index]
        if character == '"'
            next_index = nextind(line, index)
            if quoted && next_index <= lastindex(line) && line[next_index] == '"'
                write(current, '"')
                index = nextind(line, next_index)
                continue
            end
            quoted = !quoted
        elseif character == ',' && !quoted
            push!(fields, String(take!(current)))
        else
            write(current, character)
        end
        index = nextind(line, index)
    end
    quoted && throw(ArgumentError("unterminated quoted CSV field: $line"))
    push!(fields, String(take!(current)))
    return fields
end

function _read_csv_rows(path::AbstractString)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("CSV file is empty: $path"))
    header = _parse_csv_line(only(first(lines, 1)))
    rows = Dict{String,String}[]
    for line in @view lines[2:end]
        isempty(line) && continue
        fields = _parse_csv_line(line)
        length(fields) == length(header) || throw(ArgumentError(
            "CSV row has $(length(fields)) fields but header has $(length(header)): $path",
        ))
        push!(rows, Dict(zip(header, fields)))
    end
    return rows
end

"""
    read_result(directory)

Read the normalized public result schema without reconstructing Julia
functions, a `System`, or MPO tensors. This is intentionally an analysis and
provenance reader, not a restart mechanism.
"""
function read_result(directory::AbstractString)
    input = abspath(directory)
    isdir(input) || throw(ArgumentError("result directory does not exist: $input"))
    required = ("input.toml", "metadata.toml", "observables.toml", "scf_history.csv", "site_density.csv", "bond_order.csv")
    for filename in required
        isfile(joinpath(input, filename)) || throw(ArgumentError(
            "result directory is missing required $filename: $input",
        ))
    end
    metadata = _read_public_toml(joinpath(input, "metadata.toml"))
    get(metadata, "format_version", nothing) == _PUBLIC_RESULT_FORMAT_VERSION || throw(ArgumentError(
        "unsupported public result format in $input: $(get(metadata, "format_version", "missing"))",
    ))
    history = _read_csv_rows(joinpath(input, "scf_history.csv"))
    density_rows = _read_csv_rows(joinpath(input, "site_density.csv"))
    bond_rows = _read_csv_rows(joinpath(input, "bond_order.csv"))
    density = [parse(Float64, row["density"]) for row in density_rows]
    bonds = [
        (
            site_left=parse(Int, row["site_left"]),
            site_right=parse(Int, row["site_right"]),
            orientation=Symbol(row["orientation"]),
            value=ComplexF64(parse(Float64, row["real"]), parse(Float64, row["imag"])),
        ) for row in bond_rows
    ]
    return StoredResult(
        input, _read_public_toml(joinpath(input, "input.toml")), metadata,
        _read_public_toml(joinpath(input, "observables.toml")), history, density, bonds,
    )
end

"""
    CaseSpec(; label, model, representation=QTTSettings(), solver=SCFSettings(), runtime=RuntimeSettings(),
               spectral_bounds=nothing)

One named public calculation in a [`CampaignSpec`](@ref). Function-valued
model terms remain Julia functions in the campaign source; this object never
attempts to serialize or infer their envelopes. `spectral_bounds` is relevant
only to `method=:mpo` and is deliberately omitted for `method=:dense`.
"""
struct CaseSpec
    label::String
    model::AbstractPhysicalModel
    representation::QTTSettings
    solver::SCFSettings
    kpm::KPMSettings
    runtime::RuntimeSettings
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}
end

function CaseSpec(
    ; label::AbstractString,
    model::AbstractPhysicalModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings(),
    kpm::KPMSettings=KPMSettings(),
    runtime::RuntimeSettings=RuntimeSettings(),
    spectral_bounds=nothing,
)
    isempty(strip(label)) && throw(ArgumentError("case label must not be empty"))
    bounds = isnothing(spectral_bounds) ? nothing : validate_spectral_bounds(spectral_bounds...)
    return CaseSpec(String(label), model, representation, solver, kpm, runtime, bounds)
end

"""
    CampaignSpec(; name, cases)

A reproducible list of public cases for the generic campaign CLI. It is
deliberately a Julia object, rather than a data-only format, so hopping,
potential, and seed functions keep their ordinary Julia definitions.
"""
struct CampaignSpec
    name::String
    cases::Vector{CaseSpec}
end

function CampaignSpec(; name::AbstractString, cases::AbstractVector{<:CaseSpec})
    isempty(strip(name)) && throw(ArgumentError("campaign name must not be empty"))
    isempty(cases) && throw(ArgumentError("campaign must contain at least one case"))
    labels = getfield.(cases, :label)
    length(unique(labels)) == length(labels) || throw(ArgumentError(
        "campaign case labels must be unique",
    ))
    return CampaignSpec(String(name), collect(cases))
end

function campaign_case(campaign::CampaignSpec, task::Integer)
    1 <= task <= length(campaign.cases) || throw(ArgumentError(
        "task=$task is outside the campaign range 1:$(length(campaign.cases))",
    ))
    return campaign.cases[Int(task)]
end

function preflight(
    campaign::CampaignSpec;
    task::Integer,
    method::Symbol=:mpo,
    runtime::RuntimeSettings=campaign_case(campaign, task).runtime,
)
    case = campaign_case(campaign, task)
    return preflight(
        case.model, case.representation, case.solver;
        runtime=runtime, method=method,
        spectral_bounds=method == :mpo ? case.spectral_bounds : nothing,
        kpm=case.kpm,
    )
end

_campaign_path_component(text::AbstractString) = begin
    component = replace(String(text), r"[^A-Za-z0-9._-]+" => "_")
    isempty(component) && throw(ArgumentError("campaign path component must contain a safe character"))
    component
end

"""
    run_campaign(campaign; task, method=:mpo, output_root, verbose=:nothing,
                 source_path=nothing)

Run one public campaign case and write its normalized, non-overwriting result
directory at `output_root/<campaign>/task_XXXX_<label>`. `source_path`, when
provided by the CLI, is copied as immutable `campaign.jl` provenance in
addition to the descriptive `input.toml` written by [`write_result`](@ref).

`method=:dense` is CPU-only. `method=:mpo` accepts `RuntimeSettings` and, on
CUDA, uses the currently supported hybrid path: purification and effective
Hamiltonian work execute on the device, while Hartree/Fock extraction and
observable measurement currently execute on the host.
`method=:kpm` accepts open `SquareModel` lattices and explicit `GraphModel`
graphs, using the per-case `KPMSettings`.
"""
function run_campaign(
    campaign::CampaignSpec;
    task::Integer,
    method::Symbol=:mpo,
    output_root::AbstractString,
    verbose::Symbol=:nothing,
    runtime::RuntimeSettings=campaign_case(campaign, task).runtime,
    source_path::Union{Nothing,AbstractString}=nothing,
    defer_observables::Bool=false,
    write_checkpoint::Bool=method == :mpo,
)
    case = campaign_case(campaign, task)
    report = preflight(campaign; task=task, method=method, runtime=runtime)
    report.runnable || throw(ArgumentError(
        "campaign case cannot run:\n" * _preflight_error_message(report),
    ))
    source = isnothing(source_path) ? nothing : abspath(source_path)
    !isnothing(source) && !isfile(source) && throw(ArgumentError(
        "campaign source does not exist: $source",
    ))
    root = abspath(output_root)
    isdir(root) || mkpath(root)
    directory = joinpath(
        root,
        _campaign_path_component(campaign.name),
        "task_" * lpad(string(task), 4, '0') * "_" * _campaign_path_component(case.label),
    )
    ispath(directory) && throw(ArgumentError(
        "refusing to overwrite existing result directory: $directory",
    ))
    defer_observables && method != :mpo && throw(ArgumentError(
        "defer_observables is currently supported only for method=:mpo",
    ))
    checkpoint_staging = write_checkpoint ? directory * ".checkpoint.h5" : nothing
    !isnothing(checkpoint_staging) && ispath(checkpoint_staging) && throw(ArgumentError(
        "refusing to overwrite existing checkpoint: $checkpoint_staging",
    ))
    result = if method == :mpo
        solve(
            case.model, case.representation, case.solver;
            runtime=runtime, method=method, spectral_bounds=case.spectral_bounds, verbose=verbose,
            checkpoint_path=checkpoint_staging,
            measure_observables=!defer_observables,
        )
    elseif method == :kpm
        solve(case.model, case.representation, case.solver;
            runtime=runtime, method=method, kpm=case.kpm, verbose=verbose)
    else
        solve(case.model, case.representation, case.solver; runtime=runtime, method=method)
    end
    write_result(result, directory)
    if !isnothing(checkpoint_staging) && isfile(checkpoint_staging)
        mv(checkpoint_staging, joinpath(directory, "converged_state.h5"); force=false)
    end
    if !isnothing(source)
        cp(source, joinpath(directory, "campaign.jl"); force=false)
    end
    return directory, result
end

function SCFSettings(
    ; purification::Symbol=:sp2,
    mixing::Real=0.5,
    tolerance::Real=0.1,
    maxiter::Integer=30,
    purification_maxiter::Integer=50,
    square_fock_method::Symbol=:binary_carry,
    sp2_idempotency_tolerance::Real=1e-3,
    sp2_relative_trace_tolerance::Real=1e-6,
    record_energy::Bool=true,
    stable_iterations::Integer=2,
    require_stationarity::Bool=true,
    measure_stationarity::Bool=true,
    detect_two_cycles::Bool=true,
    mixing_method::Symbol=:linear,
    pulay_history::Integer=4,
    pulay_warmup::Integer=4,
    pulay_regularization::Real=1e-12,
    pulay_coefficient_limit::Real=8.0,
    pulay_step_limit::Real=20.0,
)
    purification in (:sp2, :palser_manolopoulos, :adaptive_pm_mcweeny, :mcweeny_mu) ||
        throw(ArgumentError("unsupported purification method: $purification"))
    isfinite(mixing) && 0.0 <= mixing <= 1.0 || throw(ArgumentError(
        "mixing must be finite and lie in [0, 1], got $mixing",
    ))
    isfinite(tolerance) && tolerance > 0 || throw(ArgumentError(
        "tolerance must be finite and positive, got $tolerance",
    ))
    0 < maxiter <= typemax(Int) || throw(ArgumentError(
        "maxiter must be a positive Int-compatible integer, got $maxiter",
    ))
    0 < purification_maxiter <= typemax(Int) || throw(ArgumentError(
        "purification_maxiter must be a positive Int-compatible integer, got $purification_maxiter",
    ))
    square_fock_method in (:binary_carry, :tci) || throw(ArgumentError(
        "square_fock_method must be :binary_carry or :tci, got $square_fock_method",
    ))
    isfinite(sp2_idempotency_tolerance) && sp2_idempotency_tolerance > 0 || throw(ArgumentError(
        "sp2_idempotency_tolerance must be finite and positive, got $sp2_idempotency_tolerance",
    ))
    isfinite(sp2_relative_trace_tolerance) && sp2_relative_trace_tolerance > 0 || throw(ArgumentError(
        "sp2_relative_trace_tolerance must be finite and positive, got $sp2_relative_trace_tolerance",
    ))
    stable_iterations > 0 || throw(ArgumentError(
        "stable_iterations must be positive, got $stable_iterations",
    ))
    require_stationarity && !measure_stationarity && throw(ArgumentError(
        "measure_stationarity=false requires require_stationarity=false",
    ))
    mixing_method in (:linear, :pulay) || throw(ArgumentError(
        "mixing_method must be :linear or :pulay, got $mixing_method",
    ))
    pulay_history >= 2 || throw(ArgumentError("pulay_history must be at least 2"))
    2 <= pulay_warmup <= pulay_history || throw(ArgumentError(
        "pulay_warmup must lie in 2:pulay_history, got $pulay_warmup",
    ))
    isfinite(pulay_regularization) && pulay_regularization >= 0 || throw(ArgumentError(
        "pulay_regularization must be finite and nonnegative, got $pulay_regularization",
    ))
    isfinite(pulay_coefficient_limit) && pulay_coefficient_limit > 0 || throw(ArgumentError(
        "pulay_coefficient_limit must be finite and positive, got $pulay_coefficient_limit",
    ))
    isfinite(pulay_step_limit) && pulay_step_limit > 0 || throw(ArgumentError(
        "pulay_step_limit must be finite and positive, got $pulay_step_limit",
    ))
    return SCFSettings(
        purification, Float64(mixing), Float64(tolerance), Int(maxiter), Int(purification_maxiter),
        square_fock_method, Float64(sp2_idempotency_tolerance),
        Float64(sp2_relative_trace_tolerance), record_energy, Int(stable_iterations),
        require_stationarity, measure_stationarity, detect_two_cycles,
        mixing_method, Int(pulay_history), Int(pulay_warmup),
        Float64(pulay_regularization), Float64(pulay_coefficient_limit),
        Float64(pulay_step_limit),
    )
end

"""Return the number of binary/QTT bits required by the physical chain size."""
qtt_bits(model::ChainModel) = trailing_zeros(model.size)

"""Return the `(Lx, Ly)` QTT-bit counts for a two-dimensional physical model."""
qtt_bits(model::SquareModel) = (trailing_zeros(model.size[1]), trailing_zeros(model.size[2]))

"""Return the binary storage bit count for an explicit power-of-two graph."""
qtt_bits(model::GraphModel) = trailing_zeros(size(model.hopping, 1))

"""Return the total number of QTT sites used by the current implementation."""
qtt_levels(model::ChainModel) = qtt_bits(model)
qtt_levels(model::SquareModel) = sum(qtt_bits(model))
qtt_levels(model::GraphModel) = qtt_bits(model)

function _resolved_encoding(model::ChainModel, representation::QTTSettings)
    representation.encoding in (:auto, :binary) || throw(ArgumentError(
        "ChainModel supports only encoding=:binary, got $(representation.encoding)",
    ))
    return :binary
end

function _resolved_encoding(model::SquareModel, representation::QTTSettings)
    representation.encoding in (:auto, :interleaved) || throw(ArgumentError(
        "SquareModel supports only encoding=:interleaved, got $(representation.encoding)",
    ))
    return :interleaved
end

function _resolved_encoding(model::GraphModel, representation::QTTSettings)
    representation.encoding in (:auto, :binary) || throw(ArgumentError(
        "GraphModel supports only encoding=:binary for its storage order, got $(representation.encoding)",
    ))
    return :binary
end

"""
    legacy_parameters(model, representation, solver)

Convert a public configuration into the existing parameter type without
changing its physical or numerical values. Existing `System` and `run_scf!`
continue to use this conversion during the migration.
"""
function legacy_parameters(
    model::ChainModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings(),
)
    _resolved_encoding(model, representation)
    return Parameters1D(
        L=qtt_levels(model), t=model.hopping, U=model.interaction,
        W=model.potential, S=model.seed, tci_tol=representation.tci_tol,
        itensors_tol=representation.cutoff, itensors_maxdim=representation.maxdim,
        density=model.filling, purification_steps=solver.purification_maxiter,
        scf_mixing=solver.mixing, scf_tol=solver.tolerance,
        scf_max_iterations=solver.maxiter,
    )
end

function legacy_parameters(
    model::SquareModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings(),
)
    _resolved_encoding(model, representation)
    nx, ny = model.size
    nx == ny || throw(ArgumentError(
        "the current legacy ParametersSquare MPO/dense core supports only equal squares; " *
        "received SquareModel(size=($nx, $ny)). The public geometry is valid, but no " *
        "end-to-end rectangular MPO/dense backend is implemented yet.",
    ))
    return ParametersSquare(
        L=qtt_levels(model), t=model.hopping, U=model.interaction,
        W=model.potential, S=model.seed, tci_tol=representation.tci_tol,
        itensors_tol=representation.cutoff, itensors_maxdim=representation.maxdim,
        density=model.filling, purification_steps=solver.purification_maxiter,
        scf_mixing=solver.mixing, scf_tol=solver.tolerance,
        scf_max_iterations=solver.maxiter,
    )
end

function legacy_parameters(
    model::GraphModel,
    representation::QTTSettings=QTTSettings(),
    solver::SCFSettings=SCFSettings(),
)
    _resolved_encoding(model, representation)
    throw(ArgumentError(
        "GraphModel has no legacy Parameters1D/ParametersSquare conversion; use method=:kpm.",
    ))
end

"""Convert legacy chain parameters into public model/representation/solver settings."""
function public_configuration(params::Parameters1D; purification::Symbol=:sp2)
    size = Base.Checked.checked_pow(2, params.L)
    model = ChainModel(
        size=size, hopping=params.t, interaction=params.U, potential=params.W,
        seed=params.S, filling=params.density,
    )
    representation = QTTSettings(
        encoding=:binary, tci_tol=params.tci_tol, cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    solver = SCFSettings(
        purification=purification, mixing=params.scf_mixing, tolerance=params.scf_tol,
        maxiter=params.scf_max_iterations, purification_maxiter=params.purification_steps,
    )
    return (model=model, representation=representation, solver=solver)
end

"""Convert legacy square parameters into public model/representation/solver settings."""
function public_configuration(params::ParametersSquare; purification::Symbol=:sp2)
    side = Base.Checked.checked_pow(2, div(params.L, 2))
    model = SquareModel(
        size=(side, side), hopping=params.t, interaction=params.U, potential=params.W,
        seed=params.S, filling=params.density,
    )
    representation = QTTSettings(
        encoding=:interleaved, tci_tol=params.tci_tol, cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    solver = SCFSettings(
        purification=purification, mixing=params.scf_mixing, tolerance=params.scf_tol,
        maxiter=params.scf_max_iterations, purification_maxiter=params.purification_steps,
    )
    return (model=model, representation=representation, solver=solver)
end
