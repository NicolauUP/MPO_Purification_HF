#!/usr/bin/env julia

"""Validate QTT charge FFT and D2 against a directly diagonalized projector.

The validation deliberately separates four layers:
1. exact diagonalization of a small open square Hamiltonian;
2. compression of its exact occupied projector as an MPO;
3. local extraction of the density-diagonal QTT;
4. QTT Fourier transformation and dyadic D2 analysis.

No Hartree--Fock self-consistency is performed here.
"""

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using TOML

length(ARGS) in (1, 2) || error(
    "usage: $(PROGRAM_FILE) OUTPUT_DIRECTORY [SIDE_BITS=4]",
)
output = abspath(ARGS[1])
side_bits = length(ARGS) == 2 ? parse(Int, ARGS[2]) : 4
side_bits >= 2 || error("SIDE_BITS must be at least 2")
ispath(output) && error("refusing to overwrite existing output: $output")

levels = 2side_bits
side = 1 << side_bits
N = side^2
Ne = N ÷ 2
tau = sqrt(2.0) - 5.0 / 6.0
hopping_amplitude = 0.1
seed_amplitude = 2.0
tx(x, y) = -1.0 - hopping_amplitude * cos(2pi * tau * (Float64(x) + 0.5))
ty(x, y) = -1.0 - hopping_amplitude * cos(2pi * tau * (Float64(y) + 0.5))
seed(x, y) = iseven(Int(x) + Int(y)) ? seed_amplitude : -seed_amplitude

params = ParametersSquare(
    L=levels,
    t=(tx, ty),
    U=0.0,
    W=nothing,
    S=seed,
    tci_tol=1e-12,
    itensors_tol=1e-12,
    itensors_maxdim=N,
    density=0.5,
    purification_steps=30,
    scf_mixing=0.5,
    scf_tol=1e-8,
    scf_max_iterations=10,
)

function dense_hamiltonian()
    H = zeros(Float64, N, N)
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_lattice_index(x, y, levels)
        H[site, site] = seed(x, y)
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, levels)
            H[site, neighbour] = H[neighbour, site] = tx(x, y)
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, levels)
            H[site, neighbour] = H[neighbour, site] = ty(x, y)
        end
    end
    return H
end

function projector_tensor(matrix::AbstractMatrix, sites)
    values = zeros(eltype(matrix), ntuple(_ -> 2, 2levels))
    coordinates = Vector{Int}(undef, 2levels)
    for row in 1:N, column in 1:N
        for position in 1:levels
            shift = levels - position
            coordinates[2position - 1] = ((row - 1) >> shift & 1) + 1
            coordinates[2position] = ((column - 1) >> shift & 1) + 1
        end
        values[coordinates...] = matrix[row, column]
    end
    indices = Index[]
    for site in sites
        push!(indices, prime(site), dag(site))
    end
    return ITensor(values, indices...)
end

function direct_fourier(values::AbstractVector)
    transformed = Matrix{ComplexF64}(undef, side, side)
    for kx in 0:(side - 1), ky in 0:(side - 1)
        transformed[kx + 1, ky + 1] = sum(
            values[square_lattice_index(x, y, levels)] *
            exp(-2im * pi * (kx * x + ky * y) / side)
            for x in 0:(side - 1), y in 0:(side - 1)
        ) / side
    end
    return transformed
end

function direct_z2(amplitudes; scales=0:side_bits)
    mass = abs2.(amplitudes)
    mass ./= sum(mass)
    rows = NamedTuple[]
    for scale in scales
        box_side = 1 << (side_bits - scale)
        box_masses = [sum(
            mass[square_lattice_index(x, y, levels)]
            for x in bx:(bx + box_side - 1), y in by:(by + box_side - 1)
        ) for bx in 0:box_side:(side - 1), by in 0:box_side:(side - 1)]
        push!(rows, (; scale, z2=sum(abs2, box_masses)))
    end
    return rows
end

function direct_fourier_z2(amplitudes; scales=0:side_bits)
    mass = abs2.(amplitudes)
    mass ./= sum(mass)
    rows = NamedTuple[]
    for scale in scales
        box_side = 1 << (side_bits - scale)
        box_masses = [sum(
            mass[kx + 1, ky + 1]
            for kx in bx:(bx + box_side - 1), ky in by:(by + box_side - 1)
        ) for bx in 0:box_side:(side - 1), by in 0:box_side:(side - 1)]
        push!(rows, (; scale, z2=sum(abs2, box_masses)))
    end
    return rows
end

function fit_d2(rows, fit_scales)
    x = [-scale * log(2.0) for scale in fit_scales]
    y = [log(rows[scale + 1].z2) for scale in fit_scales]
    xmean = sum(x) / length(x)
    ymean = sum(y) / length(y)
    slope = sum((x[i] - xmean) * (y[i] - ymean) for i in eachindex(x)) /
            sum((value - xmean)^2 for value in x)
    return slope
end

csv(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_row(io, values) = println(io, join(csv.(values), ','))
rms(values) = sqrt(sum(abs2, values) / length(values))

system = System(params)
dense_ed = @timed begin
    decomposition = eigen(Symmetric(dense_hamiltonian()))
    occupied = @view decomposition.vectors[:, 1:Ne]
    projector = Matrix(occupied * occupied')
    (; decomposition, projector)
end
exact_projector = dense_ed.value.projector
exact_density = real.(diag(exact_projector))
exact_centered = exact_density .- sum(exact_density) / N

compression = @timed MPO(
    projector_tensor(exact_projector, system.sites), system.sites;
    cutoff=1e-13, maxdim=N,
)
system.ρ = compression.value
density = density_diagonal_mps(system)
centered, mean_density, trace_value = centered_density_mps(
    density, system.sites; cutoff=1e-13, maxdim=N,
)
qtt_density = [real(qtt_mps_amplitude(density, system.sites, index)) for index in 0:(N - 1)]
density_error = qtt_density - exact_density

qtt_fft_timing = @timed qtt_fourier_square(
    centered, system.sites;
    sign=-1, cutoff_MPO=1e-13, cutoff=1e-12, maxdim=N,
)
qtt_fft = qtt_fft_timing.value
direct_fft_timing = @timed direct_fourier(exact_centered)
direct_fft = direct_fft_timing.value
qtt_fft_values = [qtt_fourier_amplitude(qtt_fft, kx, ky)
                  for kx in 0:(side - 1), ky in 0:(side - 1)]
fft_error = qtt_fft_values - direct_fft

fit_scales = collect(1:side_bits)
real_qtt = qtt_multiscale_d2(
    centered; cutoff=1e-13, maxdim=N, keep=:prefix, fit_scales=fit_scales,
)
fourier_qtt = qtt_multiscale_d2(
    qtt_fft; cutoff=1e-13, maxdim=N, keep=:suffix, fit_scales=fit_scales,
)
real_direct = direct_z2(exact_centered)
fourier_direct = direct_fourier_z2(direct_fft)
real_direct_d2 = fit_d2(real_direct, fit_scales)
fourier_direct_d2 = fit_d2(fourier_direct, fit_scales)

parseval_direct = sum(abs2, exact_centered)
parseval_qtt = real(inner(qtt_fft, qtt_fft))

temporary = output * ".tmp.$(getpid())"
mkpath(temporary)
try
    open(joinpath(temporary, "metadata.toml"), "w") do io
        TOML.print(io, Dict(
            "created_at" => string(now()),
            "diagnostic" => "qtt_fft_d2_vs_exact_diagonalization",
            "side" => side,
            "matrix_dimension" => N,
            "target_particles" => Ne,
            "hopping_amplitude" => hopping_amplitude,
            "seed_amplitude" => seed_amplitude,
            "exact_fermi_gap" => dense_ed.value.decomposition.values[Ne + 1] -
                                 dense_ed.value.decomposition.values[Ne],
            "fit_scales" => fit_scales,
            "fourier_normalization" => "1/sqrt(N)",
            "momentum_interval" => "[0, 2pi)",
        ))
    end
    open(joinpath(temporary, "summary.toml"), "w") do io
        TOML.print(io, Dict(
            "exact_diagonalization_time_s" => dense_ed.time,
            "projector_mpo_compression_time_s" => compression.time,
            "projector_mpo_max_chi" => maxlinkdim(system.ρ),
            "density_qtt_max_chi" => maxlinkdim(density),
            "density_trace" => real(trace_value),
            "density_mean" => real(mean_density),
            "density_max_abs_error" => maximum(abs, density_error),
            "density_rms_error" => rms(density_error),
            "direct_fourier_time_s" => direct_fft_timing.time,
            "qtt_fourier_time_s" => qtt_fft_timing.time,
            "qtt_fourier_max_chi" => maxlinkdim(qtt_fft),
            "fourier_max_abs_error" => maximum(abs, fft_error),
            "fourier_rms_error" => rms(fft_error),
            "parseval_relative_error" => abs(parseval_qtt - parseval_direct) / parseval_direct,
            "real_d2_exact" => real_direct_d2,
            "real_d2_qtt" => real_qtt.d2,
            "real_d2_abs_error" => abs(real_qtt.d2 - real_direct_d2),
            "real_d2_qtt_r_squared" => real_qtt.r_squared,
            "fourier_d2_exact" => fourier_direct_d2,
            "fourier_d2_qtt" => fourier_qtt.d2,
            "fourier_d2_abs_error" => abs(fourier_qtt.d2 - fourier_direct_d2),
            "fourier_d2_qtt_r_squared" => fourier_qtt.r_squared,
        ))
    end
    open(joinpath(temporary, "scales.csv"), "w") do io
        write_row(io, ("space", "scale", "epsilon", "box_linear_size", "z2_exact", "z2_qtt", "abs_error", "local_d2_qtt"))
        for (space, exact_rows, qtt_result) in (
            ("real", real_direct, real_qtt),
            ("fourier", fourier_direct, fourier_qtt),
        )
            for scale in 0:side_bits
                exact_z2 = exact_rows[scale + 1].z2
                qtt_row = qtt_result.scales[scale + 1]
                local_d2 = scale == 0 ? NaN : qtt_result.local_slopes[scale].d2
                write_row(io, (space, scale, qtt_row.epsilon, qtt_row.box_linear_size,
                    exact_z2, qtt_row.z2, abs(exact_z2 - qtt_row.z2), local_d2))
            end
        end
    end
    mv(temporary, output; force=false)
catch
    ispath(temporary) && rm(temporary; recursive=true, force=true)
    rethrow()
end

println("QTT FFT/D2 validation complete: $output")
println("density_max_abs_error=$(maximum(abs, density_error))")
println("fourier_max_abs_error=$(maximum(abs, fft_error))")
println("parseval_relative_error=$(abs(parseval_qtt - parseval_direct) / parseval_direct)")
println("real_D2 exact=$real_direct_d2 qtt=$(real_qtt.d2)")
println("fourier_D2 exact=$fourier_direct_d2 qtt=$(fourier_qtt.d2)")
