"""QTT-native charge-density diagnostics for square MPO states.

These routines analyze the diagonal field `n(i) = rho[i,i]`. They do not
compute the density-matrix (orbital) inverse participation ratio. Instead,
their IPR is the participation ratio of the centered charge modulation
`delta_n(i) = n(i) - Tr(rho)/N`.
"""

"""
    density_diagonal_mps(sys)

Project the equal bra/ket physical bits of `sys.rho` locally and return the
charge density as an MPS/QTT. No lattice sites are enumerated.
"""
function density_diagonal_mps(sys::System)
    return MPS(_density_diagonal_qtt_tensors(sys.ρ, sys.sites))
end

function _ones_mps(sites::Vector{<:Index}, ::Type{T}=Float64) where {T<:Number}
    return MPS([ITensor(T[one(T), one(T)], site) for site in sites])
end

"""
    centered_density_mps(density, sites; cutoff, maxdim)

Return `(delta_density, mean_density, trace)`, where
`delta_density(i) = density(i) - trace/N`. The trace is obtained by a single
MPS contraction with the all-ones product state.
"""
function centered_density_mps(
    density::MPS,
    sites::Vector{<:Index};
    cutoff::Real,
    maxdim::Integer,
)
    length(density) == length(sites) || throw(ArgumentError(
        "density and site counts must agree",
    ))
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    ones_state = _ones_mps(sites, eltype(density[1]))
    trace_value = inner(ones_state, density)
    mean_density = trace_value / (1 << length(sites))
    centered = +(
        density,
        -mean_density * ones_state;
        cutoff=cutoff,
        maxdim=maxdim,
    )
    return centered, mean_density, trace_value
end

"""Return the (possibly complex) amplitude of an MPS at zero-based QTT index `i`."""
function qtt_mps_amplitude(
    state::MPS,
    sites::Vector{<:Index},
    index::Integer,
)
    return inner(_qtt_basis_mps(sites, index), state)
end

"""
    qtt_charge_ipr(delta_density; cutoff, maxdim)

Compute

`IPR = sum_i |delta_n(i)|^4 / (sum_i |delta_n(i)|^2)^2`

from QTT contractions. `participation = 1/IPR` is the effective number of
sites carrying the charge modulation. The elementwise square is explicitly
compressed with the supplied policy; `relative_fourth_moment_error` should be
estimated independently by repeating the analysis at a stricter policy.
"""
function qtt_charge_ipr(
    delta_density::MPS;
    cutoff::Real,
    maxdim::Integer,
)
    cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    norm2 = real(inner(delta_density, delta_density))
    norm2 > 0 || throw(ArgumentError("centered density has zero norm"))
    power = Quantics.automul(
        conj(delta_density), delta_density;
        cutoff=cutoff,
        maxdim=maxdim,
    )
    fourth_moment = real(inner(power, power))
    ipr = fourth_moment / norm2^2
    return (
        norm2=norm2,
        fourth_moment=fourth_moment,
        ipr=ipr,
        participation=inv(ipr),
        power_max_chi=maxlinkdim(power),
    )
end

function _contract_physical_site!(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    position::Int,
)
    site = sites[position]
    contracted = tensors[position] * ITensor(
        fill(one(eltype(tensors[position])), dim(site)), site,
    )
    if length(tensors) == 1
        tensors[position] = contracted
    elseif position == 1
        tensors[2] = contracted * tensors[2]
        deleteat!(tensors, 1)
        deleteat!(sites, 1)
    else
        tensors[position - 1] *= contracted
        deleteat!(tensors, position)
        deleteat!(sites, position)
    end
    return tensors
end

"""Return the box-mass MPS after summing fine physical bits."""
function _qtt_box_masses(
    measure::MPS,
    kept_sites::Integer;
    keep::Symbol,
)
    0 <= kept_sites <= length(measure) || throw(ArgumentError(
        "kept_sites must lie between zero and the MPS length",
    ))
    keep in (:prefix, :suffix) || throw(ArgumentError(
        "keep must be :prefix or :suffix",
    ))
    tensors = [copy(tensor) for tensor in measure]
    sites = collect(siteinds(measure))
    if keep == :prefix
        for position in length(tensors):-1:(kept_sites + 1)
            _contract_physical_site!(tensors, sites, position)
        end
    else
        for _ in 1:(length(tensors) - kept_sites)
            _contract_physical_site!(tensors, sites, 1)
        end
    end
    if kept_sites == 0
        scalar = only(Array(only(tensors)))
        return scalar
    end
    return MPS(tensors)
end

function _linear_fit(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    length(x) == length(y) >= 2 || throw(ArgumentError(
        "a linear fit requires at least two paired points",
    ))
    xmean = sum(x) / length(x)
    ymean = sum(y) / length(y)
    denominator = sum((value - xmean)^2 for value in x)
    denominator > 0 || throw(ArgumentError("fit coordinates have zero variance"))
    slope = sum((x[i] - xmean) * (y[i] - ymean) for i in eachindex(x)) / denominator
    intercept = ymean - slope * xmean
    residual = sum((y[i] - (intercept + slope * x[i]))^2 for i in eachindex(x))
    total = sum((value - ymean)^2 for value in y)
    r_squared = total == 0 ? 1.0 : 1.0 - residual / total
    return (; slope, intercept, r_squared)
end

"""
    qtt_multiscale_d2(amplitudes; cutoff, maxdim, keep=:prefix, fit_scales)

Treat `|amplitudes|^2` as a positive normalized measure and evaluate the
dyadic correlation sum

`Z2(s) = sum_b mu_b(s)^2`

for every scale `s=0,...,log2(side)`. The returned fit has
`D2 = tau(2)`, obtained from `log(Z2)` versus `log(epsilon)` with
`epsilon=2^-s`.

Use `keep=:prefix` when coarse coordinate bits occur first in the MPS (real
space). Quantics' Fourier transform emits the coordinate bits in reverse
order, so use `keep=:suffix` for its momentum-space output.
"""
function qtt_multiscale_d2(
    amplitudes::MPS;
    cutoff::Real,
    maxdim::Integer,
    keep::Symbol=:prefix,
    fit_scales=nothing,
)
    levels = length(amplitudes)
    iseven(levels) || throw(ArgumentError(
        "square multiscale analysis requires an even number of QTT sites",
    ))
    coordinate_bits = levels ÷ 2
    power = Quantics.automul(
        conj(amplitudes), amplitudes;
        cutoff=cutoff,
        maxdim=maxdim,
    )
    total_mass = real(inner(_ones_mps(siteinds(power), eltype(power[1])), power))
    total_mass > 0 || throw(ArgumentError("amplitude measure has zero total mass"))
    measure = inv(total_mass) * power

    scales = NamedTuple[]
    for scale in 0:coordinate_bits
        boxes = _qtt_box_masses(measure, 2scale; keep=keep)
        z2 = scale == 0 ? abs2(boxes) : real(inner(boxes, boxes))
        push!(scales, (
            scale=scale,
            box_linear_size=1 << (coordinate_bits - scale),
            epsilon=2.0^(-scale),
            boxes=1 << (2scale),
            z2=z2,
            log_epsilon=-scale * log(2.0),
            log_z2=log(z2),
        ))
    end

    selected = isnothing(fit_scales) ? collect(1:coordinate_bits) : collect(fit_scales)
    length(selected) >= 2 || throw(ArgumentError("fit_scales must contain at least two scales"))
    all(scale -> 1 <= scale <= coordinate_bits, selected) || throw(ArgumentError(
        "fit scales must lie in 1:$coordinate_bits",
    ))
    fit = _linear_fit(
        [scales[scale + 1].log_epsilon for scale in selected],
        [scales[scale + 1].log_z2 for scale in selected],
    )
    local_slopes = [(
        left_scale=scale - 1,
        right_scale=scale,
        d2=(scales[scale + 1].log_z2 - scales[scale].log_z2) /
           (scales[scale + 1].log_epsilon - scales[scale].log_epsilon),
    ) for scale in 1:coordinate_bits]
    return (
        d2=fit.slope,
        intercept=fit.intercept,
        r_squared=fit.r_squared,
        fit_scales=selected,
        total_mass=total_mass,
        finest_ipr=scales[end].z2,
        participation=inv(scales[end].z2),
        power_max_chi=maxlinkdim(power),
        measure_max_chi=maxlinkdim(measure),
        scales=scales,
        local_slopes=local_slopes,
    )
end

"""
    qtt_fourier_square(state, sites; sign=-1, cutoff_MPO, cutoff, maxdim)

Apply the normalized two-dimensional QTT Fourier transform to an interleaved
square state. Odd MPS positions are the `x` coordinate bits and even positions
are the `y` bits. Parseval therefore uses the same norm before and after the
transform.

Quantics stores each transformed coordinate with its output bits reversed
along the MPS chain. [`qtt_fourier_amplitude`](@ref) hides this representation
detail when querying ordinary momentum indices.
"""
function qtt_fourier_square(
    state::MPS,
    sites::Vector{<:Index};
    sign::Int=-1,
    cutoff_MPO::Real=1e-12,
    cutoff::Real=1e-10,
    maxdim::Integer=512,
)
    length(state) == length(sites) || throw(ArgumentError(
        "state and site counts must agree",
    ))
    iseven(length(sites)) || throw(ArgumentError(
        "interleaved square Fourier transform requires an even number of QTT sites",
    ))
    transformed_x = Quantics.fouriertransform(
        state;
        sign=sign,
        sitessrc=sites[1:2:end],
        cutoff_MPO=cutoff_MPO,
        cutoff=cutoff,
        maxdim=maxdim,
    )
    return Quantics.fouriertransform(
        transformed_x;
        sign=sign,
        sitessrc=sites[2:2:end],
        cutoff_MPO=cutoff_MPO,
        cutoff=cutoff,
        maxdim=maxdim,
    )
end

function _reverse_low_bits(value::Integer, bits::Integer)
    result = 0
    for bit in 0:(bits - 1)
        result |= ((Int(value) >> bit) & 1) << (bits - bit - 1)
    end
    return result
end

"""
    qtt_fourier_amplitude(transformed, kx, ky)

Query the normalized Fourier coefficient at integer momentum `(kx, ky)`, with
`kx,ky = 0,...,side-1`. Its physical momentum is `2pi*k/side`, in `[0,2pi)`.
"""
function qtt_fourier_amplitude(
    transformed::MPS,
    kx::Integer,
    ky::Integer,
)
    levels = length(transformed)
    iseven(levels) || throw(ArgumentError("square Fourier MPS must have even length"))
    bits = levels ÷ 2
    side = 1 << bits
    0 <= kx < side || throw(BoundsError("kx=$kx lies outside 0:$(side - 1)"))
    0 <= ky < side || throw(BoundsError("ky=$ky lies outside 0:$(side - 1)"))
    # Quantics' QFT emits least-significant output bits first for each selected
    # coordinate. Convert ordinary integer momenta to that positional order.
    stored_x = _reverse_low_bits(kx, bits)
    stored_y = _reverse_low_bits(ky, bits)
    index = square_lattice_index(stored_x, stored_y, levels) - 1
    return qtt_mps_amplitude(transformed, siteinds(transformed), index)
end
