

abstract type AbstractModelParameters end #This is amazing?!


Base.@kwdef struct Parameters1D{Tt, Tu, Tw, Ts} <: AbstractModelParameters
    L::Int
    t::Tt #Type of the hopping, maybe either a number or a function to be used with TCI
    U::Tu #Type of the interaction
    W::Tw #Type of the potential.
    S::Ts #Type of the seed for TCI
    tci_tol::Float64
    itensors_tol::Float64
    itensors_maxdim::Int
    density::Float64
    purification_steps::Int
    scf_mixing::Float64
    scf_tol::Float64
    scf_max_iterations::Int
end
    
Base.@kwdef struct ParametersSquare{Tt<:Tuple, Tu, Tw, Ts} <: AbstractModelParameters
    L::Int
    t::Tt #Type of the hopping, maybe either a number or a function to be used with TCI
    U::Tu #Type of the interaction
    W::Tw #Type of the potential.
    S::Ts #Type of the seed for TCI
    tci_tol::Float64
    itensors_tol::Float64
    itensors_maxdim::Int
    density::Float64
    purification_steps::Int
    scf_mixing::Float64
    scf_tol::Float64
    scf_max_iterations::Int
end

"""A compact scalar record from one SCF iteration; it deliberately stores no MPOs.

`hartree_bond_dimension` and `fock_bond_dimension` refer to the fields freshly
extracted from the recorded density. `effective_hamiltonian_bond_dimension`
refers to the Hamiltonian that was purified to obtain that density.
"""
struct SCFIterationRecord
    iteration::Int
    trace::Float64
    vh_residual::Float64
    vf_residual::Float64
    rho_residual::Float64
    commutator_residual::Float64
    two_cycle_residual::Float64
    purification_converged::Bool
    purification_termination_reason::Symbol
    purification_iterations::Int
    rho_bond_dimension::Union{Nothing,Int}
    hartree_bond_dimension::Union{Nothing,Int}
    fock_bond_dimension::Union{Nothing,Int}
    effective_hamiltonian_bond_dimension::Union{Nothing,Int}
    energy_total::Union{Nothing,Float64}
end

"""Final status and compact scalar history from the most recent `run_scf!` call."""
struct SCFDiagnostics
    history::Vector{SCFIterationRecord}
    converged::Bool
    termination_reason::Symbol
end


mutable struct System{P}
    params::P
    sites::Vector{Index{Int64}}
    H0::MPO # This is just H0 + W
    VH::MPO # Dynamic: Hartree Potential
    VF::MPO # Dynamic: Fock Potential
    ρ::MPO # Dynamic: Density Matrix

    translations::Any # Immutable translation MPOs for this fixed basis
    bra_states::Any
    ket_states::Any
    scf_diagnostics::SCFDiagnostics
end

function _validate_common_parameters(params::AbstractModelParameters)
    params.L > 0 || throw(ArgumentError("L must be positive, got $(params.L)"))
    params.U isa Real || throw(ArgumentError("U must be a real scalar; spatially dependent and complex interactions are not implemented"))
    isfinite(params.U) || throw(ArgumentError("U must be finite, got $(params.U)"))
    isfinite(params.tci_tol) && params.tci_tol > 0 || throw(ArgumentError("tci_tol must be finite and positive, got $(params.tci_tol)"))
    isfinite(params.itensors_tol) && params.itensors_tol > 0 || throw(ArgumentError("itensors_tol must be finite and positive, got $(params.itensors_tol)"))
    params.itensors_maxdim > 0 || throw(ArgumentError("itensors_maxdim must be positive, got $(params.itensors_maxdim)"))
    params.purification_steps > 0 || throw(ArgumentError("purification_steps must be positive, got $(params.purification_steps)"))
    isfinite(params.scf_mixing) && 0.0 <= params.scf_mixing <= 1.0 || throw(ArgumentError("scf_mixing must be finite and lie in [0, 1], got $(params.scf_mixing)"))
    isfinite(params.scf_tol) && params.scf_tol > 0 || throw(ArgumentError("scf_tol must be finite and positive, got $(params.scf_tol)"))
    params.scf_max_iterations > 0 || throw(ArgumentError("scf_max_iterations must be positive, got $(params.scf_max_iterations)"))
    isfinite(params.density) && 0.0 < params.density < 1.0 || throw(ArgumentError("density must be finite and lie strictly between 0 and 1; empty and full fillings are not implemented"))

    N = try
        Base.Checked.checked_pow(2, params.L)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("L=$(params.L) is too large to represent N=2^L with Int"))
    end
    Ne = round(Int, N * params.density)
    0 < Ne < N || throw(ArgumentError("density=$(params.density) rounds to unsupported occupation Ne=$Ne for N=$N"))
    return nothing
end

function validate_parameters(params::Parameters1D)
    _validate_common_parameters(params)
    (params.t isa Number || params.t isa Function) || throw(ArgumentError("1D hopping t must be a number or function"))
    params.t isa Number && !isfinite(params.t) && throw(ArgumentError("numeric hopping t must be finite"))
    (isnothing(params.W) || params.W isa Function) || throw(ArgumentError("W must be nothing or a function"))
    (isnothing(params.S) || params.S isa Function) || throw(ArgumentError("S must be nothing or a function"))
    return nothing
end

function validate_parameters(params::ParametersSquare)
    _validate_common_parameters(params)
    iseven(params.L) || throw(ArgumentError("ParametersSquare requires even L, got $(params.L)"))
    length(params.t) == 2 || throw(ArgumentError("square hopping t must contain exactly (t_x, t_y)"))
    all(component -> component isa Number || component isa Function, params.t) ||
        throw(ArgumentError("each square hopping component must be a number or function"))
    all(component -> !(component isa Number) || isfinite(component), params.t) ||
        throw(ArgumentError("numeric square hopping components must be finite"))
    (isnothing(params.W) || params.W isa Function) || throw(ArgumentError("W must be nothing or a function"))
    (isnothing(params.S) || params.S isa Function) || throw(ArgumentError("S must be nothing or a function"))
    return nothing
end


function System(params::AbstractModelParameters)
    validate_parameters(params)
    sites = ITensors.siteinds("Qubit", params.L)
    

    translations = params isa Parameters1D ?
        build_translation_chain(sites) : build_translation_square(sites)
    H_static = build_H0(sites, params; translations=translations)
    VH_init = build_seed(sites, params)
    VF_init = Identity_MPO(sites) * 0.0 # We start with no Fock potential, but we could also build a seed for it.
    rho_init = Identity_MPO(sites) * 0.0 #nothing
    bra, ket = precompute_qtt_states(sites)
    diagnostics = SCFDiagnostics(SCFIterationRecord[], false, :not_run)
    
    return System{typeof(params)}(
        params,
        sites,
        H_static,
        VH_init,
        VF_init,
        rho_init,
        translations,
        bra,
        ket,
        diagnostics,
    )
end

"""Return diagnostics from the most recent call to `run_scf!` on `sys`."""
scf_diagnostics(sys::System) = sys.scf_diagnostics

function Base.show(io::IO, sys::System)
    println(io, "System (L=$(sys.params.L))")
    println(io, "  ├─ Hopping (t): $(typeof(sys.params.t))")
    println(io, "  ├─ Interaction (U): $(typeof(sys.params.U))")
    println(io, "  ├─ Potential (W): $(isnothing(sys.params.W) ? "None" : typeof(sys.params.W))")
    println(io, "  ├─ TCI Precision: $(sys.params.tci_tol)")
    println(io, "  ├─ ITensors Precision: $(sys.params.itensors_tol)")
    println(io, "  └─ ITensors MaxDim: $(sys.params.itensors_maxdim)")
    println(io, "  [Dynamic State]")
    println(io, "  ├─ VH MaxLinkDim: $(maxlinkdim(sys.VH))")
    println(io, "  ├─ VF MaxLinkDim: $(maxlinkdim(sys.VF))")
    println(io, "  └─ ρ MaxLinkDim:  $(maxlinkdim(sys.ρ))")
end
