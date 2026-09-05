"""Internal helpers for the fixed-H QTT KPM MPS-probe feasibility diagnostic.

These routines deliberately do *not* implement an SCF solver. They provide a
rank-one Walsh--Hadamard MPS probe and the same Chebyshev recurrence used by
the sparse-vector KPM path. The diagnostic can therefore separate the error
introduced by MPO--MPS compression from the independent finite-probe error.
"""

function _qtt_hadamard_probe_mps(sites::Vector{<:Index}, row::Integer)
    levels = length(sites)
    0 <= row < (1 << levels) || throw(ArgumentError(
        "Hadamard probe row $row is outside 0:$((1 << levels) - 1)",
    ))
    gray_row = xor(Int(row), Int(row) >> 1)
    tensors = Vector{ITensor}(undef, levels)
    for position in 1:levels
        bit = levels - position
        sign = isodd(count_ones(gray_row & (1 << bit))) ? -1.0 : 1.0
        # This product state is intentionally unnormalised: its amplitudes
        # are ±1, exactly matching the vector KPM Hadamard convention.
        tensors[position] = ITensor([1.0, sign], sites[position])
    end
    return MPS(tensors)
end

function _qtt_hadamard_probe_vector(levels::Integer, row::Integer)
    N = 1 << levels
    0 <= row < N || throw(ArgumentError("Hadamard probe row is outside the QTT basis"))
    gray_row = xor(Int(row), Int(row) >> 1)
    return [isodd(count_ones(index & gray_row)) ? -1.0 : 1.0 for index in 0:(N - 1)]
end

function _qtt_mps_mean_linkdim(psi::MPS)
    length(psi) <= 1 && return 1.0
    links = filter(!isnothing, [
        commonind(psi[position], psi[position + 1]) for position in 1:(length(psi) - 1)
    ])
    isempty(links) && return 1.0
    return sum(dim, links) / length(links)
end

function _qtt_basis_mps(sites::Vector{<:Index}, index::Integer)
    levels = length(sites)
    0 <= index < (1 << levels) || throw(BoundsError("QTT basis index $index"))
    states = Vector{String}(undef, levels)
    for position in 1:levels
        bit = levels - position
        states[position] = ((Int(index) >> bit) & 1) == 0 ? "0" : "1"
    end
    return MPS(sites, states)
end

function _qtt_mps_amplitudes(psi::MPS, sites::Vector{<:Index})
    N = 1 << length(sites)
    return [_qtt_mps_amplitude(psi, sites, index) for index in 0:(N - 1)]
end

_qtt_mps_amplitude(psi::MPS, sites::Vector{<:Index}, index::Integer) =
    real(inner(_qtt_basis_mps(sites, index), psi))

"""
    _qtt_hadamard_weighted_mps(psi, sites, row)

Return the MPS whose amplitude at QTT index `i` is `z_row(i) * psi(i)`, where
`z_row` is the unnormalised Walsh--Hadamard probe. Since `z_row` factorises
over quantics bits, the multiplication is an exact on-site diagonal operation
and does not enlarge any MPS link dimension.
"""
function _qtt_hadamard_weighted_mps(
    psi::MPS,
    sites::Vector{<:Index},
    row::Integer,
)
    length(psi) == length(sites) || throw(ArgumentError("MPS and site counts must agree"))
    levels = length(sites)
    0 <= row < (1 << levels) || throw(ArgumentError("Hadamard row is outside the QTT basis"))
    gray_row = xor(Int(row), Int(row) >> 1)
    weighted = copy(psi)
    for position in eachindex(sites)
        bit = levels - position
        sign = isodd(count_ones(gray_row & (1 << bit))) ? -1.0 : 1.0
        sign == 1.0 && continue
        site = sites[position]
        output_site = sim(site)
        diagonal = ITensor(output_site, site)
        diagonal[output_site => 1, site => 1] = 1.0
        diagonal[output_site => 2, site => 2] = sign
        weighted[position] = replaceind(weighted[position] * diagonal, output_site => site)
    end
    return weighted
end

"""
    _qtt_density_mps_from_hadamard_probes(states, sites, rows; cutoff, maxdim)

Construct the finite-probe KPM density estimator directly as an MPS:
`n(i) = mean(z_row(i) * states[row](i))`. This does not use TCI or evaluate
individual sites; `cutoff` and `maxdim` apply only to summing the MPS terms.
"""
function _qtt_density_mps_from_hadamard_probes(
    states::AbstractVector{<:MPS},
    sites::Vector{<:Index},
    rows::AbstractVector{<:Integer};
    cutoff::Real,
    maxdim::Integer,
)
    length(states) == length(rows) || throw(ArgumentError("states and rows must have equal length"))
    isempty(states) && throw(ArgumentError("at least one propagated state is required"))
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    scale = 1.0 / length(states)
    density = scale * _qtt_hadamard_weighted_mps(states[1], sites, rows[1])
    for position in 2:length(states)
        term = scale * _qtt_hadamard_weighted_mps(states[position], sites, rows[position])
        density = +(density, term; cutoff=cutoff, maxdim=maxdim)
    end
    return density
end

"""
    _qtt_hadamard_probe_block_mps(sites, rows; cutoff, maxdim)

Stack several rank-one Walsh--Hadamard probes into one block MPS with an
explicit probe-label index. The returned state has amplitudes
`Z(i, slot) = z_rows[slot](i)`. An MPO acting only on QTT site indices can
therefore propagate every block column through the same Chebyshev recurrence.
"""
function _qtt_hadamard_probe_block_mps(
    sites::Vector{<:Index},
    rows::AbstractVector{<:Integer};
    cutoff::Real,
    maxdim::Integer,
)
    isempty(rows) && throw(ArgumentError("at least one probe row is required"))
    length(unique(rows)) == length(rows) || throw(ArgumentError("probe rows must be distinct"))
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    cutoff > 0 || throw(ArgumentError("cutoff must be positive"))
    probe_index = Index(length(rows), "QTTProbe,Block")
    block = nothing
    for (slot, row) in enumerate(rows)
        probe = _qtt_hadamard_probe_mps(sites, row)
        probe[1] *= onehot(probe_index => slot)
        block = isnothing(block) ? probe : +(block, probe; cutoff=cutoff, maxdim=maxdim)
    end
    return block, probe_index
end

"""
    _qtt_hadamard_probe_register_mps(sites, probe_bits)

Construct a binary QTT register for `R = 2^probe_bits` Hadamard probes. The
MPS site order leaves the leading spatial bits unchanged and interleaves the
least-significant bits as `(…,i_{L-q+1},a₁,…,i_L,a_q)`, where `i` is the
spatial QTT coordinate and `a` is the Hadamard-code register. Its amplitudes
are `(-1)^popcount(i & a)`. This is the same *set* of columns as the existing
Gray-ordered first-`R` probes, merely in a different order.

The pairwise form has initial link dimension at most two; it avoids creating a
dense `N × R` array solely to initialise the probe block.
"""
function _qtt_hadamard_probe_register_mps(
    sites::Vector{<:Index},
    probe_bits::Integer,
)
    1 <= probe_bits <= length(sites) || throw(ArgumentError(
        "probe_bits must lie in 1:$(length(sites))",
    ))
    probe_sites = [Index(2, "Qubit,QTTProbe,n=$position") for position in 1:probe_bits]
    tensors = ITensor[]
    previous_bridge = nothing
    first_paired = length(sites) - probe_bits + 1
    for position in 1:(first_paired - 1)
        tensor = ITensor([1.0, 1.0], sites[position])
        !isnothing(previous_bridge) && (tensor *= onehot(previous_bridge => 1))
        next_bridge = Index(1, "Link,QTTProbe,n=$position")
        tensor *= onehot(next_bridge => 1)
        previous_bridge = next_bridge
        push!(tensors, tensor)
    end
    for probe_position in 1:probe_bits
        spatial_position = first_paired + probe_position - 1
        spatial_site, probe_site = sites[spatial_position], probe_sites[probe_position]
        pair = ITensor(spatial_site, probe_site)
        pair[spatial_site => 1, probe_site => 1] = 1.0
        pair[spatial_site => 1, probe_site => 2] = 1.0
        pair[spatial_site => 2, probe_site => 1] = 1.0
        pair[spatial_site => 2, probe_site => 2] = -1.0
        left, singular_values, right = svd(pair, [spatial_site])
        !isnothing(previous_bridge) && (left *= onehot(previous_bridge => 1))
        push!(tensors, left)
        pair_right = singular_values * right
        if spatial_position < length(sites)
            next_bridge = Index(1, "Link,QTTProbe,n=$spatial_position")
            pair_right *= onehot(next_bridge => 1)
            previous_bridge = next_bridge
        else
            previous_bridge = nothing
        end
        push!(tensors, pair_right)
    end
    ordered_sites = Index[]
    for position in 1:length(sites)
        push!(ordered_sites, sites[position])
        probe_position = position - first_paired + 1
        1 <= probe_position <= probe_bits && push!(ordered_sites, probe_sites[probe_position])
    end
    return MPS(tensors), probe_sites, ordered_sites
end

"""Insert spectator probe-register sites after the least-significant spatial QTT sites."""
function _qtt_extend_mpo_with_probe_register(
    H::MPO,
    spatial_sites::Vector{<:Index},
    probe_sites::Vector{<:Index},
)
    length(H) == length(spatial_sites) || throw(ArgumentError("MPO/site counts must agree"))
    length(probe_sites) <= length(spatial_sites) || throw(ArgumentError(
        "probe register cannot have more bits than the spatial QTT coordinate",
    ))
    source = [copy(H[position]) for position in 1:length(H)]
    tensors = ITensor[]
    first_paired = length(spatial_sites) - length(probe_sites) + 1
    for position in 1:length(source)
        probe_position = position - first_paired + 1
        if !(1 <= probe_position <= length(probe_sites))
            push!(tensors, source[position])
            continue
        end
        if position == length(source)
            final_link = Index(1, "Link,QTTProbe,final")
            source[position] *= onehot(final_link => 1)
            push!(tensors, source[position])
            push!(tensors, delta(prime(probe_sites[probe_position]), probe_sites[probe_position]) *
                onehot(final_link => 1))
            continue
        end
        push!(tensors, source[position])
        old_link = commonind(source[position], source[position + 1])
        isnothing(old_link) && throw(ArgumentError("neighbouring MPO tensors lack a common link"))
        new_link = sim(old_link)
        source[position + 1] = replaceind(source[position + 1], old_link => new_link)
        identity_probe = delta(prime(probe_sites[probe_position]), probe_sites[probe_position]) *
            delta(old_link, new_link)
        push!(tensors, identity_probe)
    end
    return MPO(tensors)
end

"""Evaluate the probe-register state at spatial QTT index `i` and code `a`."""
function _qtt_probe_register_amplitude(
    psi::MPS,
    spatial_sites::Vector{<:Index},
    probe_sites::Vector{<:Index},
    i::Integer,
    a::Integer,
)
    levels, probe_bits = length(spatial_sites), length(probe_sites)
    0 <= i < (1 << levels) || throw(BoundsError("spatial index $i is outside the QTT basis"))
    0 <= a < (1 << probe_bits) || throw(BoundsError("probe code $a is outside the register basis"))
    accumulator = ITensor(1.0)
    tensor_position = 1
    first_paired = levels - probe_bits + 1
    for position in 1:levels
        bit = levels - position
        accumulator *= psi[tensor_position] *
            onehot(spatial_sites[position] => (((Int(i) >> bit) & 1) + 1))
        tensor_position += 1
        probe_position = position - first_paired + 1
        if 1 <= probe_position <= probe_bits
            probe_bit = probe_bits - probe_position
            accumulator *= psi[tensor_position] *
                onehot(probe_sites[probe_position] => (((Int(a) >> probe_bit) & 1) + 1))
            tensor_position += 1
        end
    end
    return real(scalar(accumulator))
end

"""Return the three-index copy tensor `δ_{coefficient,input,output}`."""
function _qtt_copy_tensor(coefficient::Index, input::Index, output::Index)
    dim(coefficient) == dim(input) == dim(output) || throw(ArgumentError(
        "copy-tensor indices must have equal dimension",
    ))
    tensor = ITensor(coefficient, input, output)
    for value in 1:dim(coefficient)
        tensor[coefficient => value, input => value, output => value] = 1.0
    end
    return tensor
end

"""
    _qtt_diagonal_mpo_from_mps(function_mps, sites)

Promote a scalar function represented as an MPS to a diagonal MPO. The MPO
acts pointwise: `apply(D_f, ψ)` represents `f(x) ψ(x)`. Keeping this
conversion explicit avoids evaluating the function on all `2^L` coordinates.
"""
function _qtt_diagonal_mpo_from_mps(function_mps::MPS, sites::Vector{<:Index})
    length(function_mps) == length(sites) || throw(ArgumentError(
        "function-MPS/site counts must agree",
    ))
    tensors = ITensor[]
    for position in eachindex(sites)
        site = sites[position]
        coefficient_site = sim(site)
        coefficient_tensor = replaceind(function_mps[position], site => coefficient_site)
        push!(tensors, coefficient_tensor * _qtt_copy_tensor(
            coefficient_site, site, prime(site),
        ))
    end
    return MPO(tensors)
end

"""Contract a physical MPS site against weights and remove that site exactly."""
function _qtt_contract_mps_site(
    tensors::Vector{ITensor},
    site::Index,
    weights::AbstractVector{<:Number},
)
    length(weights) == dim(site) || throw(ArgumentError("weight/site dimensions must agree"))
    position = findfirst(tensor -> hasind(tensor, site), tensors)
    isnothing(position) && throw(ArgumentError("MPS does not contain the requested site"))
    weight_tensor = ITensor(collect(weights), site)
    collapsed = tensors[position] * weight_tensor
    if position == length(tensors)
        position == 1 && return ITensor[collapsed]
        tensors[position - 1] *= collapsed
    else
        tensors[position + 1] = collapsed * tensors[position + 1]
    end
    deleteat!(tensors, position)
    return tensors
end

"""
    _qtt_register_weighted_spatial_mps(state, weights, spatial_sites, probe_sites; cutoff, maxdim, scale)

For a register MPS `state(i,a)` and a register-function MPS `weights(i,a)`,
form the spatial MPS

That is, `f(i) = scale * sum_a weights(i,a) * state(i,a)`.

This is an exact tensor-network contraction apart from the explicit `apply`
compression. It is the primitive used for density and directed bond fields.
"""
function _qtt_register_weighted_spatial_mps(
    state::MPS,
    weights::MPS,
    spatial_sites::Vector{<:Index},
    probe_sites::Vector{<:Index};
    cutoff::Real,
    maxdim::Integer,
    scale::Real=1.0,
)
    length(state) == length(weights) || throw(ArgumentError("register MPS lengths must agree"))
    ordered_sites = Index[]
    first_paired = length(spatial_sites) - length(probe_sites) + 1
    for position in eachindex(spatial_sites)
        push!(ordered_sites, spatial_sites[position])
        probe_position = position - first_paired + 1
        1 <= probe_position <= length(probe_sites) && push!(ordered_sites, probe_sites[probe_position])
    end
    weighted = apply(
        _qtt_diagonal_mpo_from_mps(weights, ordered_sites), state;
        cutoff=cutoff, maxdim=maxdim,
    )
    tensors = [copy(weighted[position]) for position in eachindex(weighted)]
    for probe_site in reverse(probe_sites)
        _qtt_contract_mps_site(tensors, probe_site, ones(Float64, dim(probe_site)))
    end
    return scale * MPS(tensors)
end

"""Construct `n(i)=R⁻¹Σₐ zₐ(i)[P(H)zₐ](i)` directly as a spatial QTT MPS."""
function _qtt_density_mps_from_probe_register(
    propagated::MPS,
    initial_register::MPS,
    spatial_sites::Vector{<:Index},
    probe_sites::Vector{<:Index};
    cutoff::Real,
    maxdim::Integer,
)
    probes = 1 << length(probe_sites)
    return _qtt_register_weighted_spatial_mps(
        propagated, initial_register, spatial_sites, probe_sites;
        cutoff=cutoff, maxdim=maxdim, scale=inv(probes),
    )
end

"""
    _qtt_directed_bond_mps_from_probe_register(propagated, initial_register, translation, spatial_sites, probe_sites; cutoff, maxdim)

Return the directed finite-probe estimator
`rho_delta(i) = R^-1 sum_a [P z_a](i) z_a(i+delta)`. `translation` must use
the project's row convention, i.e. `(translation * f)[i] = f[i+delta]`; open
boundary rows are therefore zero automatically. The translation MPO is lifted
with identity action on the probe register before it is applied.
"""
function _qtt_directed_bond_mps_from_probe_register(
    propagated::MPS,
    initial_register::MPS,
    translation::MPO,
    spatial_sites::Vector{<:Index},
    probe_sites::Vector{<:Index};
    cutoff::Real,
    maxdim::Integer,
)
    extended_translation = _qtt_extend_mpo_with_probe_register(
        translation, spatial_sites, probe_sites,
    )
    shifted_register = apply(
        extended_translation, initial_register; cutoff=cutoff, maxdim=maxdim,
    )
    probes = 1 << length(probe_sites)
    return _qtt_register_weighted_spatial_mps(
        propagated, shifted_register, spatial_sites, probe_sites;
        cutoff=cutoff, maxdim=maxdim, scale=inv(probes),
    )
end

"""
    _qtt_square_fock_mpo_from_bond_mps(horizontal, vertical, sites, translations, U; cutoff, maxdim)

Assemble the real nearest-neighbour exchange field from spatial QTT MPS bond
coefficients. This preserves the repository's established convention
`(V_F)_{i,i+delta} = -U * real(rho_{i,i+delta})` and its open-boundary
translation orientation. No density-matrix MPO is required.
"""
function _qtt_square_fock_mpo_from_bond_mps(
    horizontal::MPS,
    vertical::MPS,
    sites::Vector{<:Index},
    translations,
    U::Real;
    cutoff::Real,
    maxdim::Integer,
)
    length(horizontal) == length(sites) == length(vertical) || throw(ArgumentError(
        "bond-MPS/site counts must agree",
    ))
    T_R, T_L, T_U, T_D = translations
    horizontal_coefficients = -Float64(U) * _qtt_diagonal_mpo_from_mps(horizontal, sites)
    vertical_coefficients = -Float64(U) * _qtt_diagonal_mpo_from_mps(vertical, sites)
    horizontal_field = +(
        apply(horizontal_coefficients, T_R; cutoff=cutoff, maxdim=maxdim),
        apply(T_L, ITensors.dag(horizontal_coefficients); cutoff=cutoff, maxdim=maxdim);
        cutoff=cutoff, maxdim=maxdim,
    )
    vertical_field = +(
        apply(vertical_coefficients, T_U; cutoff=cutoff, maxdim=maxdim),
        apply(T_D, ITensors.dag(vertical_coefficients); cutoff=cutoff, maxdim=maxdim);
        cutoff=cutoff, maxdim=maxdim,
    )
    return +(
        horizontal_field, vertical_field; cutoff=cutoff, maxdim=maxdim,
    ), (; horizontal=horizontal_field, vertical=vertical_field)
end

"""
    _qtt_probe_register_moments(H, register, degree; cutoff, maxdim)

Estimate Chebyshev moments using the complete binary probe register:
`mu_m = R^-1 sum_a <z_a|T_m(H)|z_a>`. Only the current and previous
recurrence MPS are retained. This is the register analogue of sparse-vector
KPM moments and supplies the trace-based chemical-potential solve in SCF.
"""
function _qtt_probe_register_moments(
    H::MPO,
    register::MPS,
    degree::Integer;
    cutoff::Real,
    maxdim::Integer,
)
    degree >= 1 || throw(ArgumentError("Chebyshev moment degree must be positive"))
    probes = 1
    for tensor in register
        for index in inds(tensor)
            hastags(index, "QTTProbe") && (probes *= dim(index))
        end
    end
    # Every probe site occurs in two adjacent tensors only if it is a bond;
    # physical register indices occur once, so the product above is exactly R.
    # The explicit tag makes this independent of spatial QTT bond dimensions.
    moments = Vector{Float64}(undef, degree + 1)
    previous = copy(register)
    moments[1] = real(inner(register, previous)) / probes
    current = apply(H, previous; cutoff=cutoff, maxdim=maxdim)
    moments[2] = real(inner(register, current)) / probes
    for order in 2:degree
        following = apply(H, current; cutoff=cutoff, maxdim=maxdim)
        following = +(2.0 * following, -previous; cutoff=cutoff, maxdim=maxdim)
        moments[order + 1] = real(inner(register, following)) / probes
        previous, current = current, following
    end
    return moments
end

"""Evaluate one `slot` of a block-MPS at a QTT basis `index`."""
function _qtt_block_mps_amplitude(
    psi::MPS,
    sites::Vector{<:Index},
    probe_index::Index,
    slot::Integer,
    index::Integer,
)
    1 <= slot <= dim(probe_index) || throw(BoundsError("block slot $slot is outside probe index"))
    position = findfirst(tensor -> hasind(tensor, probe_index), psi)
    isnothing(position) && throw(ArgumentError("block MPS does not contain the probe index"))
    basis = _qtt_basis_mps(sites, index)
    basis[position] *= onehot(probe_index => slot)
    return real(inner(basis, psi))
end

function _qtt_mps_chebyshev_apply(
    H::MPO,
    probe::MPS,
    coefficients::AbstractVector{<:Real};
    cutoff::Real,
    maxdim::Integer,
)
    length(coefficients) >= 2 || throw(ArgumentError("KPM requires at least degree one"))
    previous = copy(probe)
    current = apply(H, probe; cutoff=cutoff, maxdim=maxdim)
    result = +(
        (coefficients[1] / 2) * previous,
        coefficients[2] * current;
        cutoff=cutoff,
        maxdim=maxdim,
    )
    trajectory = NamedTuple[
        (order=0, state_max_chi=maxlinkdim(previous), state_mean_chi=_qtt_mps_mean_linkdim(previous),
         result_max_chi=maxlinkdim(result), result_mean_chi=_qtt_mps_mean_linkdim(result)),
        (order=1, state_max_chi=maxlinkdim(current), state_mean_chi=_qtt_mps_mean_linkdim(current),
         result_max_chi=maxlinkdim(result), result_mean_chi=_qtt_mps_mean_linkdim(result)),
    ]
    for order in 2:(length(coefficients) - 1)
        following = apply(H, current; cutoff=cutoff, maxdim=maxdim)
        following = +(2.0 * following, -previous; cutoff=cutoff, maxdim=maxdim)
        result = +(result, coefficients[order + 1] * following; cutoff=cutoff, maxdim=maxdim)
        push!(trajectory, (
            order=order, state_max_chi=maxlinkdim(following), state_mean_chi=_qtt_mps_mean_linkdim(following),
            result_max_chi=maxlinkdim(result), result_mean_chi=_qtt_mps_mean_linkdim(result),
        ))
        previous, current = current, following
    end
    return (; state=result, trajectory)
end
