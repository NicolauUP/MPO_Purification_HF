using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField

function checkerboard_field(amplitude)
    (x, y) -> iseven(Int(x) + Int(y)) ? amplitude : -amplitude
end

function fixed_profile_system(
    side_level::Int;
    cutoff::Float64,
    maxdim::Int,
    steps::Int,
)
    params = ParametersSquare(
        L=2side_level,
        t=(-0.6, -0.35),
        U=0.3,
        W=checkerboard_field(0.6),
        S=checkerboard_field(0.05),
        tci_tol=1e-10,
        itensors_tol=cutoff,
        itensors_maxdim=maxdim,
        density=0.5,
        purification_steps=steps,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=1,
    )
    sys = System(params)
    bounds = square_scf_spectral_bounds(
        params;
        potential_bounds=(-0.6, 0.6),
        margin=0.5,
    )
    rho0 = construct_rho_0(sys, params, bounds...; method=:sp2)
    purification = perform_purification(
        rho0,
        params;
        method=:sp2,
        spectral_bounds=bounds,
        spectral_bounds_validation=:supplied_unverified,
        sp2_idempotency_tolerance=2e-4,
        sp2_trace_tolerance=1e-6 * round(Int, params.density * 2^params.L),
        verbose=1,
        overwrite_progress=false,
    )
    purification.converged || error(
        "fixed-density SP2 preparation failed: " *
        "$(purification.termination_reason), trace_error=$(purification.trace_error), " *
        "idempotency=$(purification.idempotency_residual)",
    )
    sys.ρ = purification.rho
    return sys, purification, bounds
end

function bond_dimensions(mpo::MPO)
    length(mpo) <= 1 && return Int[]
    [dim(commonind(mpo[i], mpo[i + 1])) for i in 1:(length(mpo) - 1)]
end

function mean_bond_dimension(mpo::MPO)
    dimensions = bond_dimensions(mpo)
    isempty(dimensions) ? 1.0 : sum(dimensions) / length(dimensions)
end

function relative_mpo_error(candidate::MPO, reference::MPO)
    candidate_norm = max(0.0, real(inner(candidate, candidate)))
    reference_norm = max(0.0, real(inner(reference, reference)))
    overlap = real(inner(candidate, reference))
    difference_norm = sqrt(max(0.0, candidate_norm + reference_norm - 2overlap))
    return difference_norm / max(sqrt(reference_norm), sqrt(eps(Float64)))
end

function direct_hartree_probe_error(field::MPO, sys::System)
    side = 2^div(sys.params.L, 2)
    coordinates = unique((
        (0, 0),
        (side - 1, 0),
        (0, side - 1),
        (side - 1, side - 1),
        (div(side, 2), div(side, 2)),
    ))
    errors = Float64[]
    for (x, y) in coordinates
        site = square_lattice_index(x, y, sys.params.L)
        expected = sys.params.U * sum(
            real(MatrixChecker(
                sys.ρ,
                sys.sites,
                neighbour,
                neighbour,
                sys.bra_states,
                sys.ket_states,
            ))
            for neighbour in values(square_neighbours(site, sys.params.L))
            if !isnothing(neighbour)
        )
        observed = real(MatrixChecker(
            field,
            sys.sites,
            site,
            site,
            sys.bra_states,
            sys.ket_states,
        ))
        push!(errors, abs(observed - expected))
    end
    return maximum(errors), sum(errors) / length(errors)
end

function direct_fock_probe_rows(
    sys::System,
    horizontal_tci::MPO,
    vertical_tci::MPO,
    horizontal_carry::MPO,
    vertical_carry::MPO,
)
    side = 2^div(sys.params.L, 2)
    middle = div(side, 2)
    rows = NamedTuple[]
    valid_starts = unique((
        (0, 0),
        (side - 2, 0),
        (0, side - 1),
        (side - 2, side - 1),
        (min(middle, side - 2), middle),
    ))
    for orientation in (:horizontal, :vertical)
        tci = orientation == :horizontal ? horizontal_tci : vertical_tci
        carry = orientation == :horizontal ? horizontal_carry : vertical_carry
        for (x, y) in valid_starts
            source_x, source_y = orientation == :horizontal ? (x, y) : (y, x)
            source = square_lattice_index(source_x, source_y, sys.params.L)
            neighbour = orientation == :horizontal ?
                square_lattice_index(source_x + 1, source_y, sys.params.L) :
                square_lattice_index(source_x, source_y + 1, sys.params.L)
            expected = -sys.params.U * real(MatrixChecker(
                sys.ρ, sys.sites, source, neighbour,
                sys.bra_states, sys.ket_states,
            ))
            tci_value = real(MatrixChecker(
                tci, sys.sites, source, neighbour,
                sys.bra_states, sys.ket_states,
            ))
            carry_value = real(MatrixChecker(
                carry, sys.sites, source, neighbour,
                sys.bra_states, sys.ket_states,
            ))
            push!(rows, (
                orientation=orientation,
                probe="valid",
                source_x=source_x,
                source_y=source_y,
                neighbour_x=orientation == :horizontal ? source_x + 1 : source_x,
                neighbour_y=orientation == :horizontal ? source_y : source_y + 1,
                expected=expected,
                tci=tci_value,
                carry=carry_value,
                tci_error=abs(tci_value - expected),
                carry_error=abs(carry_value - expected),
            ))
        end

        # A carry leaving the coordinate's most-significant bit must be
        # discarded. Probe the corresponding wrapped coordinate explicitly:
        # it is not an open-boundary nearest-neighbour bond.
        for transverse in unique((0, side - 1, middle))
            source_x, source_y = orientation == :horizontal ?
                (side - 1, transverse) : (transverse, side - 1)
            neighbour_x, neighbour_y = orientation == :horizontal ?
                (0, transverse) : (transverse, 0)
            source = square_lattice_index(source_x, source_y, sys.params.L)
            neighbour = square_lattice_index(neighbour_x, neighbour_y, sys.params.L)
            tci_value = real(MatrixChecker(
                tci, sys.sites, source, neighbour,
                sys.bra_states, sys.ket_states,
            ))
            carry_value = real(MatrixChecker(
                carry, sys.sites, source, neighbour,
                sys.bra_states, sys.ket_states,
            ))
            push!(rows, (
                orientation=orientation,
                probe="open_boundary",
                source_x=source_x,
                source_y=source_y,
                neighbour_x=neighbour_x,
                neighbour_y=neighbour_y,
                expected=0.0,
                tci=tci_value,
                carry=carry_value,
                tci_error=abs(tci_value),
                carry_error=abs(carry_value),
            ))
        end
    end
    return rows
end

function truncated_adjacency_hartree(sys::System; cutoff::Float64, maxdim::Int)
    density_tensors = MPO_MeanField._density_diagonal_qtt_tensors(sys.ρ, sys.sites)
    density = MPS(density_tensors)
    adjacency = square_neighbour_adjacency_mpo(sys.sites)
    field_state = apply(adjacency, density; cutoff=cutoff, maxdim=maxdim)
    field = MPO_MeanField._diagonal_mpo_from_qtt_tensors(
        copy(field_state.data),
        sys.sites,
        sys.params;
        symmetrize=false,
    )
    return sys.params.U * field
end
