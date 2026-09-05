#!/usr/bin/env julia

"""Compare square binary-carry mean-field extraction on CPU and CUDA.

The benchmark prepares one fixed-H SP2 trajectory state on CUDA, copies that
same MPO once to the host, and then evaluates identical Hartree and horizontal/
vertical Fock kernels on both backends. Timings are synchronized and exclude
JIT warmup. Device-to-host result transfers are reported separately. Numerical
agreement is checked after both results are on the CPU.
"""

using CUDA
using Dates
using ITensors
using ITensorMPS
using LinearAlgebra
using MPO_MeanField
using TOML

CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

length(ARGS) == 4 || error(
    "usage: component_profile_cpu_cuda.jl CAMPAIGN TASK SP2_STEPS OUTPUT",
)
campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
sp2_steps = parse(Int, ARGS[3])
output = abspath(ARGS[4])
ispath(output) && error("refusing to overwrite existing output: $output")
Base.include(Main, campaign_file)
1 <= task_index <= length(campaign.cases) || error("task index is outside campaign")
case = campaign.cases[task_index]
params = legacy_parameters(case.model, case.representation, case.solver)
params isa ParametersSquare || error("benchmark requires a square QTT case")
isnothing(case.spectral_bounds) && error("campaign case requires spectral bounds")

to_device(value) = ITensors.adapt(CUDA.CuArray, value)
to_host(value) = ITensors.cpu(value)

function synchronized_timing(f, synchronize=() -> nothing)
    synchronize()
    timing = @timed begin
        value = f()
        synchronize()
        value
    end
    return timing
end

function relative_mpo_error(reference::MPO, candidate::MPO)
    nr = max(0.0, real(inner(reference, reference)))
    nc = max(0.0, real(inner(candidate, candidate)))
    overlap = real(inner(reference, candidate))
    return sqrt(max(0.0, nr + nc - 2overlap) / max(nr, eps(Float64)))
end

mean_chi(mpo::MPO) = length(mpo) <= 1 ? 1.0 :
    sum(dim(linkind(mpo, bond)) for bond in 1:(length(mpo) - 1)) / (length(mpo) - 1)

function csv(value)
    return '"' * replace(string(value), '"' => "\"\"") * '"'
end
write_row(io, values) = println(io, join(csv.(values), ','))

sys = System(params)
translations_host = sys.translations
translations_device_timing = synchronized_timing(
    () -> to_device(translations_host), CUDA.synchronize,
)
translations_device = translations_device_timing.value

# Construct exactly the affine SP2 start used by run_scf!, on the GPU.
sys.H0 = to_device(sys.H0)
sys.VH = to_device(sys.VH)
sys.VF = to_device(sys.VF)
rho_device = construct_rho_0(
    sys, params, case.spectral_bounds...; method=:sp2, to_gpu=to_device,
)
Ne = round(Int, params.density * 2^params.L)
trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)
for step in 1:sp2_steps
    rho_squared = apply(
        rho_device, rho_device;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    CUDA.synchronize()
    branch = MPO_MeanField._sp2_branch(
        real(tr(rho_device)), real(tr(rho_squared)), Ne, trace_tolerance,
    )
    global rho_device = branch == :square ? rho_squared : +(
        2.0 * rho_device, -rho_squared;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    CUDA.synchronize()
    println("prepared SP2 step $step/$sp2_steps chi=$(maxlinkdim(rho_device))")
    flush(stdout)
end

rho_to_host = synchronized_timing(() -> to_host(rho_device), CUDA.synchronize)
rho_host = rho_to_host.value

function extract_components(sys, backend::Symbol)
    if backend == :cpu
        hartree = extract_hartree_mpo_binary_carry_square_adjacency(sys)
        horizontal = extract_fock_mpo_binary_carry_square_horizontal(sys)
        vertical = extract_fock_mpo_binary_carry_square_vertical(sys)
    else
        hartree = extract_hartree_mpo_binary_carry_square_adjacency(
            sys; to_backend=to_device,
        )
        horizontal = extract_fock_mpo_binary_carry_square_horizontal(
            sys; to_backend=to_device, translations=translations_device,
        )
        vertical = extract_fock_mpo_binary_carry_square_vertical(
            sys; to_backend=to_device, translations=translations_device,
        )
    end
    fock = +(
        horizontal, vertical;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    return (hartree=hartree, horizontal=horizontal, vertical=vertical, fock=fock)
end

# Warm compilation independently. GPU warmup is synchronized before timing.
sys.ρ = rho_host
warm_cpu = extract_components(sys, :cpu)
warm_cpu = nothing
GC.gc(true)
sys.ρ = rho_device
warm_gpu = extract_components(sys, :cuda)
CUDA.synchronize()
warm_gpu = nothing
GC.gc(true)
CUDA.reclaim()

sys.ρ = rho_host
cpu = synchronized_timing(() -> extract_components(sys, :cpu))
sys.ρ = rho_device
gpu = synchronized_timing(() -> extract_components(sys, :cuda), CUDA.synchronize)
gpu_to_host = synchronized_timing(
    () -> map(to_host, gpu.value), CUDA.synchronize,
)

function audit_fock_orientation(
    orientation::Symbol,
    rho::MPO,
    cpu_field::MPO,
    gpu_field::MPO,
    sys::System,
    max_bonds::Int,
)
    side = 2^div(params.L, 2)
    total_bonds = side * (side - 1)
    sampling_stride = max(1, cld(total_bonds, max_bonds))
    expected_errors_cpu = Float64[]
    expected_errors_gpu = Float64[]
    cpu_gpu_differences = Float64[]
    hermiticity_errors_cpu = Float64[]
    hermiticity_errors_gpu = Float64[]
    max_record = (error=-Inf, x=-1, y=-1, expected=NaN, cpu=NaN, gpu=NaN)
    bond_index = 0
    for y in 0:(side - 1), x in 0:(side - 1)
        neighbour_x, neighbour_y = orientation == :horizontal ? (x + 1, y) : (x, y + 1)
        (neighbour_x < side && neighbour_y < side) || continue
        bond_index += 1
        (bond_index - 1) % sampling_stride == 0 || continue
        source = square_lattice_index(x, y, params.L)
        neighbour = square_lattice_index(neighbour_x, neighbour_y, params.L)
        expected = -params.U * real(MatrixChecker(
            rho, sys.sites, source, neighbour, sys.bra_states, sys.ket_states,
        ))
        cpu_value = real(MatrixChecker(
            cpu_field, sys.sites, source, neighbour, sys.bra_states, sys.ket_states,
        ))
        gpu_value = real(MatrixChecker(
            gpu_field, sys.sites, source, neighbour, sys.bra_states, sys.ket_states,
        ))
        cpu_reverse = real(MatrixChecker(
            cpu_field, sys.sites, neighbour, source, sys.bra_states, sys.ket_states,
        ))
        gpu_reverse = real(MatrixChecker(
            gpu_field, sys.sites, neighbour, source, sys.bra_states, sys.ket_states,
        ))
        push!(expected_errors_cpu, abs(cpu_value - expected))
        push!(expected_errors_gpu, abs(gpu_value - expected))
        difference = abs(cpu_value - gpu_value)
        push!(cpu_gpu_differences, difference)
        push!(hermiticity_errors_cpu, abs(cpu_value - cpu_reverse))
        push!(hermiticity_errors_gpu, abs(gpu_value - gpu_reverse))
        difference > max_record.error && (max_record = (
            error=difference, x=x, y=y, expected=expected,
            cpu=cpu_value, gpu=gpu_value,
        ))
    end
    rms(values) = sqrt(sum(abs2, values) / length(values))
    return (
        audited_bonds=length(cpu_gpu_differences),
        total_bonds=total_bonds,
        sampling_stride=sampling_stride,
        cpu_expected_max=maximum(expected_errors_cpu),
        cpu_expected_rms=rms(expected_errors_cpu),
        gpu_expected_max=maximum(expected_errors_gpu),
        gpu_expected_rms=rms(expected_errors_gpu),
        cpu_gpu_max=maximum(cpu_gpu_differences),
        cpu_gpu_rms=rms(cpu_gpu_differences),
        cpu_hermiticity_max=maximum(hermiticity_errors_cpu),
        gpu_hermiticity_max=maximum(hermiticity_errors_gpu),
        max_record=max_record,
    )
end

function audit_open_boundary_wraps(
    orientation::Symbol,
    cpu_field::MPO,
    gpu_field::MPO,
    sys::System,
)
    side = 2^div(params.L, 2)
    cpu_values = Float64[]
    gpu_values = Float64[]
    for transverse in 0:(side - 1)
        source_x, source_y = orientation == :horizontal ?
            (side - 1, transverse) : (transverse, side - 1)
        neighbour_x, neighbour_y = orientation == :horizontal ?
            (0, transverse) : (transverse, 0)
        source = square_lattice_index(source_x, source_y, params.L)
        neighbour = square_lattice_index(neighbour_x, neighbour_y, params.L)
        push!(cpu_values, abs(MatrixChecker(
            cpu_field, sys.sites, source, neighbour, sys.bra_states, sys.ket_states,
        )))
        push!(gpu_values, abs(MatrixChecker(
            gpu_field, sys.sites, source, neighbour, sys.bra_states, sys.ket_states,
        )))
    end
    return (cpu_max=maximum(cpu_values), gpu_max=maximum(gpu_values))
end

max_audit_bonds = parse(Int, get(ENV, "MPOHF_FIELD_AUDIT_BONDS", "4096"))
max_audit_bonds > 0 || error("MPOHF_FIELD_AUDIT_BONDS must be positive")
println("auditing up to $max_audit_bonds deterministic bonds per Fock orientation")
flush(stdout)
fock_audits = Dict{Symbol,Any}()
for (orientation, component) in ((:horizontal, :horizontal), (:vertical, :vertical))
    fock_audits[orientation] = audit_fock_orientation(
        orientation, rho_host, getproperty(cpu.value, component),
        getproperty(gpu_to_host.value, component), sys, max_audit_bonds,
    )
    fock_audits[Symbol(orientation, :_wrap)] = audit_open_boundary_wraps(
        orientation, getproperty(cpu.value, component),
        getproperty(gpu_to_host.value, component), sys,
    )
end

mkpath(output)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(Dates.now()),
        "campaign" => campaign_file,
        "case_label" => case.label,
        "matrix_dimension" => 2^params.L,
        "qtt_length" => params.L,
        "sp2_steps" => sp2_steps,
        "rho_max_chi" => maxlinkdim(rho_host),
        "rho_mean_chi" => mean_chi(rho_host),
        "cutoff" => params.itensors_tol,
        "maxdim" => params.itensors_maxdim,
        "device_name" => CUDA.name(CUDA.device()),
        "rho_to_host_time_s" => rho_to_host.time,
        "translations_to_device_time_s" => translations_device_timing.time,
        "gpu_results_to_host_time_s" => gpu_to_host.time,
    ))
end

open(joinpath(output, "summary.csv"), "w") do io
    write_row(io, (
        "component", "cpu_time_s", "cuda_time_s", "speedup",
        "relative_mpo_error", "cpu_max_chi", "cuda_max_chi",
    ))
    # The aggregate synchronized timing is the trustworthy end-to-end number.
    for component in (:all_fields,)
        write_row(io, (
            component, cpu.time, gpu.time, cpu.time / gpu.time,
            maximum(relative_mpo_error(
                getproperty(cpu.value, name), getproperty(gpu_to_host.value, name),
            ) for name in (:hartree, :horizontal, :vertical, :fock)),
            maximum(maxlinkdim(getproperty(cpu.value, name)) for name in (:hartree, :fock)),
            maximum(maxlinkdim(getproperty(gpu_to_host.value, name)) for name in (:hartree, :fock)),
        ))
    end
    for name in (:hartree, :horizontal, :vertical, :fock)
        write_row(io, (
            name, "included_in_all_fields", "included_in_all_fields", "not_separately_timed",
            relative_mpo_error(
                getproperty(cpu.value, name), getproperty(gpu_to_host.value, name),
            ),
            maxlinkdim(getproperty(cpu.value, name)),
            maxlinkdim(getproperty(gpu_to_host.value, name)),
        ))
    end
end

open(joinpath(output, "fock_local_audit.csv"), "w") do io
    write_row(io, (
        "orientation", "audited_bonds", "total_physical_bonds", "sampling_stride",
        "cpu_expected_max_abs_error",
        "cpu_expected_rms_error", "gpu_expected_max_abs_error",
        "gpu_expected_rms_error", "cpu_gpu_max_abs_difference",
        "cpu_gpu_rms_difference", "cpu_hermiticity_max_abs_error",
        "gpu_hermiticity_max_abs_error", "cpu_wrap_max_abs",
        "gpu_wrap_max_abs", "max_difference_x", "max_difference_y",
        "expected_at_max", "cpu_at_max", "gpu_at_max",
    ))
    for orientation in (:horizontal, :vertical)
        audit = fock_audits[orientation]
        wrap = fock_audits[Symbol(orientation, :_wrap)]
        record = audit.max_record
        write_row(io, (
            orientation, audit.audited_bonds, audit.total_bonds, audit.sampling_stride,
            audit.cpu_expected_max,
            audit.cpu_expected_rms, audit.gpu_expected_max,
            audit.gpu_expected_rms, audit.cpu_gpu_max, audit.cpu_gpu_rms,
            audit.cpu_hermiticity_max, audit.gpu_hermiticity_max,
            wrap.cpu_max, wrap.gpu_max, record.x, record.y,
            record.expected, record.cpu, record.gpu,
        ))
    end
end

println("CPU/CUDA mean-field benchmark complete: $output")
println("CPU=$(cpu.time)s CUDA=$(gpu.time)s speedup=$(cpu.time / gpu.time)x")
println("GPU result transfer=$(gpu_to_host.time)s rho host transfer=$(rho_to_host.time)s")
