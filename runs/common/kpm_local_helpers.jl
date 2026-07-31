using Random
using LinearAlgebra

"""Jackson factors g_0,...,g_M for a degree-M Chebyshev expansion."""
function jackson_factors(degree::Int)
    angle = π / (degree + 1)
    return [
        ((degree - order + 1) * cos(order * angle) +
         sin(order * angle) * cot(angle)) / (degree + 1)
        for order in 0:degree
    ]
end

"""Coefficients for Theta(mu-H) in the convention c0/2 + sum(c_m T_m)."""
function projector_coefficients(degree::Int, scaled_mu::Float64)
    -1 < scaled_mu < 1 ||
        error("scaled chemical potential must lie strictly inside (-1,1)")
    theta = acos(scaled_mu)
    coefficients = Vector{Float64}(undef, degree + 1)
    coefficients[1] = 2 * (π - theta) / π
    for order in 1:degree
        coefficients[order + 1] = -2sin(order * theta) / (π * order)
    end
    return coefficients .* jackson_factors(degree)
end

function probing_matrix(
    N::Int,
    probes::Int,
    method::Symbol,
    seed::Int,
)
    probes <= N || error("PROBES=$probes exceeds matrix dimension N=$N")
    rng = Xoshiro(seed)
    if method == :rademacher
        return ifelse.(rand(rng, Bool, N, probes), 1.0, -1.0)
    end
    method == :hadamard || error("unknown probing method: $method")
    ispow2(N) || error("Hadamard probes require a power-of-two dimension")
    signs = ifelse.(rand(rng, Bool, N), 1.0, -1.0)
    matrix = Matrix{Float64}(undef, N, probes)
    for column in 0:(probes - 1)
        # Gray ordering retains the nested-prefix property when PROBES grows.
        gray = xor(column, column >> 1)
        for row in 0:(N - 1)
            matrix[row + 1, column + 1] =
                isodd(count_ones(row & gray)) ? -signs[row + 1] : signs[row + 1]
        end
    end
    return matrix
end

"""Hadamard probes whose hierarchy is defined by an explicit row code.

`codes` must be a zero-based permutation of `0:N-1`. The returned rows remain
in the original storage order; only the Walsh/Hadamard color associated with
each row changes. This permits physically contiguous storage with a separate
coordinate-aware probing hierarchy.
"""
function coded_hadamard_matrix(
    codes::AbstractVector{<:Integer},
    probes::Int,
    seed::Int,
)
    N = length(codes)
    probes <= N || error("PROBES=$probes exceeds matrix dimension N=$N")
    ispow2(N) || error("Hadamard probes require a power-of-two dimension")
    sort!(collect(Int, codes)) == collect(0:(N - 1)) ||
        error("Hadamard row codes must be a permutation of 0:N-1")
    rng = Xoshiro(seed)
    signs = ifelse.(rand(rng, Bool, N), 1.0, -1.0)
    matrix = Matrix{Float64}(undef, N, probes)
    for column in 0:(probes - 1)
        gray = xor(column, column >> 1)
        for row in eachindex(codes)
            matrix[row, column + 1] =
                isodd(count_ones(codes[row] & gray)) ?
                -signs[row] : signs[row]
        end
    end
    return matrix
end

"""A contiguous column block of a nested Hadamard probing matrix."""
function coded_hadamard_block(
    codes::AbstractVector{<:Integer},
    first_column::Int,
    columns::Int,
    seed::Int,
)
    N = length(codes)
    first_column >= 1 || error("first_column must be positive")
    columns > 0 || error("columns must be positive")
    first_column + columns - 1 <= N ||
        error("requested Hadamard block exceeds N=$N")
    ispow2(N) || error("Hadamard probes require a power-of-two dimension")
    sort!(collect(Int, codes)) == collect(0:(N - 1)) ||
        error("Hadamard row codes must be a permutation of 0:N-1")
    rng = Xoshiro(seed)
    signs = ifelse.(rand(rng, Bool, N), 1.0, -1.0)
    matrix = Matrix{Float64}(undef, N, columns)
    for local_column in 1:columns
        column = first_column + local_column - 2
        gray = xor(column, column >> 1)
        for row in eachindex(codes)
            matrix[row, local_column] =
                isodd(count_ones(codes[row] & gray)) ?
                -signs[row] : signs[row]
        end
    end
    return matrix
end

hadamard_block(N::Int, first_column::Int, columns::Int, seed::Int) =
    coded_hadamard_block(collect(0:(N - 1)), first_column, columns, seed)

"""A coordinate-interleaved zero-based Hadamard row code.

The least-significant `min(x_bits, y_bits)` coordinate bits are interleaved
as `y_0, x_0, y_1, x_1, ...`; any remaining high bits are appended. The map
is a permutation of `0:(2^(x_bits+y_bits)-1)` for power-of-two rectangles.
"""
function coordinate_interleaved_code(
    x::Integer,
    y::Integer,
    x_bits::Integer,
    y_bits::Integer,
)
    x_bits >= 0 && y_bits >= 0 || error("coordinate bit counts must be nonnegative")
    0 <= x < (1 << x_bits) || error("x=$x does not fit in $x_bits bits")
    0 <= y < (1 << y_bits) || error("y=$y does not fit in $y_bits bits")
    shared_bits = min(x_bits, y_bits)
    code = 0
    for bit in 0:(shared_bits - 1)
        code |= ((y >> bit) & 1) << (2bit)
        code |= ((x >> bit) & 1) << (2bit + 1)
    end
    position = 2shared_bits
    for bit in shared_bits:(x_bits - 1)
        code |= ((x >> bit) & 1) << position
        position += 1
    end
    for bit in shared_bits:(y_bits - 1)
        code |= ((y >> bit) & 1) << position
        position += 1
    end
    return code
end

"""Safeguarded Pulay/DIIS mixer for a deterministic fixed-point map.

`pulay_update!(mixer, input, output)` treats `output = F(input)` and stores
the residual `F(input) - input`. It uses linear mixing during `warmup` calls,
then forms an affine DIIS combination of recent outputs. Ill-conditioned,
non-finite, or excessively large extrapolations fall back to linear mixing.
"""
mutable struct PulayMixer
    outputs::Vector{Vector{Float64}}
    residuals::Vector{Vector{Float64}}
    history::Int
    warmup::Int
    regularization::Float64
    coefficient_limit::Float64
    step_limit::Float64
end

function PulayMixer(; history::Int=6, warmup::Int=4,
    regularization::Real=1e-12, coefficient_limit::Real=8.0,
    step_limit::Real=20.0)
    history >= 2 || error("Pulay history must be at least 2")
    warmup >= 2 || error("Pulay warmup must be at least 2")
    regularization >= 0 || error("Pulay regularization must be nonnegative")
    coefficient_limit > 0 || error("Pulay coefficient limit must be positive")
    step_limit > 0 || error("Pulay step limit must be positive")
    return PulayMixer(
        Vector{Vector{Float64}}(), Vector{Vector{Float64}}(), history,
        warmup, Float64(regularization), Float64(coefficient_limit),
        Float64(step_limit),
    )
end

function _linear_mixed(input, output, damping)
    return (1 - damping) .* input .+ damping .* output
end

function pulay_update!(mixer::PulayMixer, input::Vector{Float64},
    output::Vector{Float64}; damping::Real=0.5)
    0 < damping <= 1 || error("Pulay damping must lie in (0, 1]")
    length(input) == length(output) || error("Pulay input/output size mismatch")
    residual = output - input
    push!(mixer.outputs, copy(output))
    push!(mixer.residuals, residual)
    while length(mixer.outputs) > mixer.history
        popfirst!(mixer.outputs)
        popfirst!(mixer.residuals)
    end

    linear = _linear_mixed(input, output, damping)
    count = length(mixer.residuals)
    count < mixer.warmup && return linear, :linear

    gram = Matrix{Float64}(undef, count, count)
    for row in 1:count, column in 1:count
        gram[row, column] = dot(mixer.residuals[row], mixer.residuals[column])
    end
    scale = max(maximum(abs, gram), 1.0)
    @inbounds for index in 1:count
        gram[index, index] += mixer.regularization * scale
    end
    system = zeros(Float64, count + 1, count + 1)
    system[1:count, 1:count] .= gram
    system[1:count, count + 1] .= 1.0
    system[count + 1, 1:count] .= 1.0
    coefficients = try
        system \ vcat(zeros(Float64, count), 1.0)
    catch
        return linear, :linear_solve_fallback
    end
    weights = coefficients[1:count]
    if !all(isfinite, weights) || maximum(abs, weights) > mixer.coefficient_limit
        return linear, :linear_coefficient_fallback
    end

    candidate = (1 - damping) .* input
    for index in 1:count
        candidate .+= damping * weights[index] .* mixer.outputs[index]
    end
    if !all(isfinite, candidate) ||
       norm(candidate - input) > mixer.step_limit * max(norm(residual), eps(Float64))
        return linear, :linear_step_fallback
    end
    return candidate, :pulay
end

function kpm_apply(
    scaled_hamiltonian,
    probes,
    coefficients;
    synchronize=() -> nothing,
)
    degree = length(coefficients) - 1
    degree >= 1 || error("KPM application degree must be at least one")
    previous = copy(probes)
    current = similar(probes)
    following = similar(probes)
    result = similar(probes)
    mul!(current, scaled_hamiltonian, probes)
    @. result =
        (coefficients[1] / 2) * previous + coefficients[2] * current
    synchronize()
    for order in 2:degree
        mul!(following, scaled_hamiltonian, current)
        @. following = 2 * following - previous
        @. result += coefficients[order + 1] * following
        previous, current, following = current, following, previous
    end
    synchronize()
    return result
end

"""Estimated traces Tr[T_m(H)] for m=0,...,degree using fixed probes."""
function kpm_trace_moments(
    scaled_hamiltonian,
    probes,
    degree::Int;
    synchronize=() -> nothing,
)
    degree >= 0 || error("degree must be nonnegative")
    normalization = inv(size(probes, 2))
    moments = Vector{Float64}(undef, degree + 1)
    previous = copy(probes)
    current = similar(probes)
    following = similar(probes)
    moments[1] = real(dot(vec(probes), vec(previous))) * normalization
    degree == 0 && return moments

    mul!(current, scaled_hamiltonian, probes)
    moments[2] = real(dot(vec(probes), vec(current))) * normalization
    for order in 2:degree
        mul!(following, scaled_hamiltonian, current)
        @. following = 2 * following - previous
        moments[order + 1] =
            real(dot(vec(probes), vec(following))) * normalization
        previous, current, following = current, following, previous
    end
    synchronize()
    return moments
end

"""Trace of the degree-M Jackson projector from Chebyshev trace moments."""
function projector_trace_from_moments(
    trace_moments::AbstractVector{<:Real},
    scaled_mu::Float64,
)
    coefficients = projector_coefficients(length(trace_moments) - 1, scaled_mu)
    return coefficients[1] * trace_moments[1] / 2 +
           dot(@view(coefficients[2:end]), @view(trace_moments[2:end]))
end

"""Find a scaled chemical potential whose estimated projector trace is target."""
function find_scaled_chemical_potential(
    trace_moments::AbstractVector{<:Real},
    target::Real;
    tolerance::Real=1e-6,
    max_iterations::Int=80,
)
    isfinite(target) || error("target trace must be finite")
    isfinite(tolerance) && tolerance > 0 ||
        error("trace tolerance must be finite and positive")
    max_iterations > 0 || error("max_iterations must be positive")
    lower = -1.0 + 1e-10
    upper = 1.0 - 1e-10
    lower_trace = projector_trace_from_moments(trace_moments, lower)
    upper_trace = projector_trace_from_moments(trace_moments, upper)
    lower_trace <= target <= upper_trace || error(
        "target trace $target is outside KPM bracket " *
        "[$lower_trace, $upper_trace]",
    )
    midpoint = 0.0
    midpoint_trace = projector_trace_from_moments(trace_moments, midpoint)
    for iteration in 1:max_iterations
        midpoint = (lower + upper) / 2
        midpoint_trace =
            projector_trace_from_moments(trace_moments, midpoint)
        abs(midpoint_trace - target) <= tolerance &&
            return (; scaled_mu=midpoint, trace=midpoint_trace, iterations=iteration)
        if midpoint_trace < target
            lower = midpoint
        else
            upper = midpoint
        end
    end
    return (;
        scaled_mu=midpoint,
        trace=midpoint_trace,
        iterations=max_iterations,
    )
end
