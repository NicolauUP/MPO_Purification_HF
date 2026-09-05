const _MPO_CHECKPOINT_FORMAT_VERSION = 1

function _checkpoint_scalar(file, name::AbstractString)
    haskey(file, name) || throw(ArgumentError("checkpoint is missing $name"))
    return read(file, name)
end

function _checkpoint_sites(mpo::MPO)
    return [
        only(filter(index -> plev(index) == 0, collect(siteinds(mpo, site))))
        for site in 1:length(mpo)
    ]
end

"""
    write_mpo_checkpoint(sys, path)

Atomically write a host-side converged MPO state. The checkpoint contains
`H0`, `VH`, `VF`, and `rho`; arbitrary model functions remain in the campaign
source and are deliberately not serialized. Existing paths are never
overwritten.
"""
function write_mpo_checkpoint(sys::System, path::AbstractString)
    output = abspath(path)
    ispath(output) && throw(ArgumentError("refusing to overwrite checkpoint: $output"))
    mkpath(dirname(output))
    temporary = output * ".tmp.$(getpid()).$(rand(UInt))"
    diagnostics = scf_diagnostics(sys)
    try
        h5open(temporary, "w") do file
            write(file, "format_version", _MPO_CHECKPOINT_FORMAT_VERSION)
            write(file, "julia_version", string(VERSION))
            write(file, "parameter_type", string(nameof(typeof(sys.params))))
            write(file, "qtt_levels", sys.params.L)
            write(file, "matrix_dimension", 2^sys.params.L)
            write(file, "interaction", Float64(sys.params.U))
            write(file, "density", sys.params.density)
            write(file, "itensors_cutoff", sys.params.itensors_tol)
            write(file, "itensors_maxdim", sys.params.itensors_maxdim)
            write(file, "scf_converged", diagnostics.converged)
            write(file, "scf_termination_reason", string(diagnostics.termination_reason))
            write(file, "scf_iterations", length(diagnostics.history))
            write(file, "H0", sys.H0)
            write(file, "VH", sys.VH)
            write(file, "VF", sys.VF)
            write(file, "rho", sys.ρ)
        end
        mv(temporary, output; force=false)
    catch
        ispath(temporary) && rm(temporary; force=true)
        rethrow()
    end
    return output
end

"""
    read_mpo_checkpoint(path, params)

Load a checkpoint using parameters reconstructed from its original campaign.
The stored scalar conventions are checked before a host-side `System` is
returned. Translation MPOs are not rebuilt because final observables require
only the stored Hamiltonian, fields, and density matrix.
"""
function read_mpo_checkpoint(path::AbstractString, params::AbstractModelParameters)
    input = abspath(path)
    isfile(input) || throw(ArgumentError("checkpoint does not exist: $input"))
    validate_parameters(params)
    loaded = h5open(input, "r") do file
        format_version = Int(_checkpoint_scalar(file, "format_version"))
        format_version == _MPO_CHECKPOINT_FORMAT_VERSION || throw(ArgumentError(
            "unsupported MPO checkpoint format $format_version in $input",
        ))
        qtt_levels = Int(_checkpoint_scalar(file, "qtt_levels"))
        qtt_levels == params.L || throw(ArgumentError(
            "checkpoint L=$qtt_levels does not match campaign L=$(params.L)",
        ))
        interaction = Float64(_checkpoint_scalar(file, "interaction"))
        interaction == Float64(params.U) || throw(ArgumentError(
            "checkpoint interaction=$interaction does not match campaign U=$(params.U)",
        ))
        density = Float64(_checkpoint_scalar(file, "density"))
        density == params.density || throw(ArgumentError(
            "checkpoint density=$density does not match campaign density=$(params.density)",
        ))
        cutoff = Float64(_checkpoint_scalar(file, "itensors_cutoff"))
        cutoff == params.itensors_tol || throw(ArgumentError(
            "checkpoint cutoff=$cutoff does not match campaign cutoff=$(params.itensors_tol)",
        ))
        maxdim = Int(_checkpoint_scalar(file, "itensors_maxdim"))
        maxdim == params.itensors_maxdim || throw(ArgumentError(
            "checkpoint maxdim=$maxdim does not match campaign maxdim=$(params.itensors_maxdim)",
        ))
        matrix_dimension = Int(_checkpoint_scalar(file, "matrix_dimension"))
        matrix_dimension == 2^params.L || throw(ArgumentError(
            "checkpoint dimension=$matrix_dimension does not match campaign dimension=$(2^params.L)",
        ))
        parameter_type = String(_checkpoint_scalar(file, "parameter_type"))
        parameter_type == string(nameof(typeof(params))) || throw(ArgumentError(
            "checkpoint parameter type $parameter_type does not match $(nameof(typeof(params)))",
        ))
        (
            H0=read(file, "H0", MPO),
            VH=read(file, "VH", MPO),
            VF=read(file, "VF", MPO),
            rho=read(file, "rho", MPO),
            converged=Bool(_checkpoint_scalar(file, "scf_converged")),
            termination_reason=Symbol(String(_checkpoint_scalar(file, "scf_termination_reason"))),
            iterations=Int(_checkpoint_scalar(file, "scf_iterations")),
            format_version=format_version,
        )
    end
    lengths = length.((loaded.H0, loaded.VH, loaded.VF, loaded.rho))
    all(==(params.L), lengths) || throw(ArgumentError(
        "checkpoint MPO lengths $lengths do not match L=$(params.L)",
    ))
    sites = _checkpoint_sites(loaded.H0)
    for (name, mpo) in (("VH", loaded.VH), ("VF", loaded.VF), ("rho", loaded.rho))
        _checkpoint_sites(mpo) == sites || throw(ArgumentError(
            "checkpoint $name physical indices do not match H0",
        ))
    end
    bra_states, ket_states = precompute_qtt_states(sites)
    diagnostics = SCFDiagnostics(
        SCFIterationRecord[], loaded.converged, loaded.termination_reason,
    )
    system = System{typeof(params)}(
        params, sites, loaded.H0, loaded.VH, loaded.VF, loaded.rho,
        nothing, bra_states, ket_states, diagnostics,
    )
    return (
        system=system,
        format_version=loaded.format_version,
        converged=loaded.converged,
        termination_reason=loaded.termination_reason,
        scf_iterations=loaded.iterations,
    )
end
