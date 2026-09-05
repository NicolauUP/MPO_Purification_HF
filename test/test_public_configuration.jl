@testset "Public model/representation/solver configuration compatibility" begin
    chain_potential = x -> 0.1 * Int(x)
    chain_seed = x -> iseven(Int(x)) ? 0.05 : -0.05
    legacy_chain = parameters_1d(
        L=4, t=x -> -1.0 - 0.1 * Int(x), U=0.7, W=chain_potential,
        S=chain_seed, tci_tol=1e-9, itensors_tol=1e-8, itensors_maxdim=48,
        density=0.375, purification_steps=23, scf_mixing=0.4,
        scf_tol=0.02, scf_max_iterations=17,
    )
    chain_configuration = public_configuration(legacy_chain; purification=:palser_manolopoulos)
    @test chain_configuration.model isa ChainModel
    @test chain_configuration.model.size == 16
    @test qtt_bits(chain_configuration.model) == 4
    @test qtt_levels(chain_configuration.model) == 4
    @test chain_configuration.representation.encoding == :binary
    @test chain_configuration.solver.purification == :palser_manolopoulos
    reconstructed_chain = legacy_parameters(
        chain_configuration.model, chain_configuration.representation, chain_configuration.solver,
    )
    @test reconstructed_chain.L == legacy_chain.L
    @test reconstructed_chain.t === legacy_chain.t
    @test reconstructed_chain.U == legacy_chain.U
    @test reconstructed_chain.W === legacy_chain.W
    @test reconstructed_chain.S === legacy_chain.S
    @test reconstructed_chain.tci_tol == legacy_chain.tci_tol
    @test reconstructed_chain.itensors_tol == legacy_chain.itensors_tol
    @test reconstructed_chain.itensors_maxdim == legacy_chain.itensors_maxdim
    @test reconstructed_chain.density == legacy_chain.density
    @test reconstructed_chain.purification_steps == legacy_chain.purification_steps
    @test reconstructed_chain.scf_mixing == legacy_chain.scf_mixing
    @test reconstructed_chain.scf_tol == legacy_chain.scf_tol
    @test reconstructed_chain.scf_max_iterations == legacy_chain.scf_max_iterations

    square_potential = (x, y) -> 0.1 * Int(x) - 0.05 * Int(y)
    square_seed = (x, y) -> iseven(Int(x) + Int(y)) ? 0.03 : -0.03
    legacy_square = parameters_square(
        L=6, t=((x, y) -> -1.0 - 0.01 * Int(x), (x, y) -> -0.8), U=0.2,
        W=square_potential, S=square_seed, tci_tol=1e-9, itensors_tol=1e-8,
        itensors_maxdim=96, density=0.5, purification_steps=31, scf_mixing=0.35,
        scf_tol=0.03, scf_max_iterations=19,
    )
    square_configuration = public_configuration(legacy_square)
    @test square_configuration.model isa SquareModel
    @test square_configuration.model.size == (8, 8)
    @test qtt_bits(square_configuration.model) == (3, 3)
    @test qtt_levels(square_configuration.model) == 6
    @test square_configuration.representation.encoding == :interleaved
    reconstructed_square = legacy_parameters(
        square_configuration.model, square_configuration.representation, square_configuration.solver,
    )
    @test reconstructed_square.L == legacy_square.L
    @test reconstructed_square.t === legacy_square.t
    @test reconstructed_square.U == legacy_square.U
    @test reconstructed_square.W === legacy_square.W
    @test reconstructed_square.S === legacy_square.S
    @test reconstructed_square.tci_tol == legacy_square.tci_tol
    @test reconstructed_square.itensors_tol == legacy_square.itensors_tol
    @test reconstructed_square.itensors_maxdim == legacy_square.itensors_maxdim
    @test reconstructed_square.density == legacy_square.density
    @test reconstructed_square.purification_steps == legacy_square.purification_steps
    @test reconstructed_square.scf_mixing == legacy_square.scf_mixing
    @test reconstructed_square.scf_tol == legacy_square.scf_tol
    @test reconstructed_square.scf_max_iterations == legacy_square.scf_max_iterations

    @test legacy_parameters(
        ChainModel(size=4, hopping=-1.0, interaction=0.0, filling=0.5),
        QTTSettings(encoding=:auto), SCFSettings(),
    ) isa Parameters1D
    @test legacy_parameters(
        SquareModel(size=(4, 4), hopping=(-1.0, -1.0), interaction=0.0, filling=0.5),
        QTTSettings(encoding=:auto), SCFSettings(),
    ) isa ParametersSquare

    @test_throws ArgumentError ChainModel(size=12, hopping=-1.0, interaction=0.0)
    rectangle = SquareModel(
        size=(8, 4), hopping=(-1.0, -1.0), interaction=0.0,
    )
    @test rectangle.size == (8, 4)
    @test qtt_bits(rectangle) == (3, 2)
    @test qtt_levels(rectangle) == 5
    @test_throws ArgumentError legacy_parameters(rectangle)
    @test_throws ArgumentError SquareModel(size=(4, 4), hopping=(-1.0, -1.0), interaction=0.0, boundary=:periodic)
    @test_throws ArgumentError legacy_parameters(
        SquareModel(size=(4, 4), hopping=(-1.0, -1.0), interaction=0.0),
        QTTSettings(encoding=:binary), SCFSettings(),
    )
    @test_throws ArgumentError QTTSettings(cutoff=0.0)
    @test_throws ArgumentError SCFSettings(mixing=1.1)
    operational_solver = SCFSettings(
        square_fock_method=:tci,
        sp2_idempotency_tolerance=2e-4,
        sp2_relative_trace_tolerance=1e-6,
        record_energy=false,
        stable_iterations=3,
        require_stationarity=false,
        measure_stationarity=false,
        detect_two_cycles=false,
    )
    @test operational_solver.square_fock_method == :tci
    @test operational_solver.sp2_idempotency_tolerance == 2e-4
    @test operational_solver.sp2_relative_trace_tolerance == 1e-6
    @test !operational_solver.record_energy
    @test operational_solver.stable_iterations == 3
    @test !operational_solver.require_stationarity
    @test !operational_solver.measure_stationarity
    @test !operational_solver.detect_two_cycles
    pulay_solver = SCFSettings(
        mixing_method=:pulay, pulay_history=4, pulay_warmup=3,
        pulay_regularization=1e-11, pulay_coefficient_limit=7.0,
        pulay_step_limit=15.0,
    )
    @test pulay_solver.mixing_method == :pulay
    @test pulay_solver.pulay_history == 4
    @test pulay_solver.pulay_warmup == 3
    @test pulay_solver.pulay_regularization == 1e-11
    @test pulay_solver.pulay_coefficient_limit == 7.0
    @test pulay_solver.pulay_step_limit == 15.0
    @test_throws ArgumentError SCFSettings(square_fock_method=:unknown)
    @test_throws ArgumentError SCFSettings(sp2_idempotency_tolerance=0.0)
    @test_throws ArgumentError SCFSettings(sp2_relative_trace_tolerance=0.0)
    @test_throws ArgumentError SCFSettings(stable_iterations=0)
    @test_throws ArgumentError SCFSettings(measure_stationarity=false)
    @test_throws ArgumentError SCFSettings(mixing_method=:unknown)
    @test_throws ArgumentError SCFSettings(pulay_history=1)
    @test_throws ArgumentError SCFSettings(pulay_history=3, pulay_warmup=4)
    @test_throws ArgumentError RuntimeSettings(backend=:unsupported)
    @test_throws ArgumentError RuntimeSettings(device_scalar_type=Float32)
    cpu_runtime = runtime_preflight(RuntimeSettings())
    @test cpu_runtime.runnable
    @test cpu_runtime.active_backend == :cpu
    @test cpu_runtime.execution_path == :cpu_end_to_end
    dense_cuda_runtime = runtime_preflight(RuntimeSettings(backend=:cuda); method=:dense)
    @test !dense_cuda_runtime.runnable
    @test occursin("CPU-only", dense_cuda_runtime.message)

    preflight_square = preflight(
        SquareModel(size=(4, 4), hopping=(-1.0, -0.5), interaction=0.3),
        QTTSettings(maxdim=128), SCFSettings();
        method=:mpo, spectral_bounds=(-4.0, 4.0),
    )
    @test preflight_square isa PreflightReport
    @test preflight_square.runnable
    @test preflight_square.geometry == :square
    @test preflight_square.physical_size == (4, 4)
    @test preflight_square.physical_sites == 16
    @test preflight_square.qtt_bit_counts == (2, 2)
    @test preflight_square.qtt_levels == 4
    @test preflight_square.encoding == :interleaved
    @test preflight_square.target_particles == 8
    @test preflight_square.spectral_bounds == (-4.0, 4.0)
    @test preflight_square.spectral_bounds_status == :supplied_unverified
    @test isempty(preflight_square.issues)
    @test occursin("physical lattice: 4 × 4", sprint(show, preflight_square))

    preflight_chain = preflight(
        ChainModel(size=8, hopping=-1.0, interaction=0.0),
        QTTSettings(), SCFSettings(), method=:mpo,
    )
    @test !preflight_chain.runnable
    @test preflight_chain.geometry == :chain
    @test preflight_chain.qtt_bit_counts == (3, 0)
    @test preflight_chain.spectral_bounds_status == :required_not_supplied
    @test any(issue -> issue.code == :spectral_bounds_required, preflight_chain.issues)

    preflight_rectangle = preflight(rectangle, QTTSettings(), SCFSettings(); method=:mpo)
    @test !preflight_rectangle.runnable
    @test preflight_rectangle.geometry == :rectangle
    @test any(issue -> issue.code == :rectangular_backend_unavailable, preflight_rectangle.issues)

    preflight_bad_encoding = preflight(
        SquareModel(size=(4, 4), hopping=(-1.0, -1.0), interaction=0.0),
        QTTSettings(encoding=:binary), SCFSettings(); method=:mpo,
    )
    @test !preflight_bad_encoding.runnable
    @test any(issue -> issue.code == :unsupported_encoding, preflight_bad_encoding.issues)
    preflight_dense_cuda = preflight(
        ChainModel(size=8, hopping=-1.0, interaction=0.0),
        QTTSettings(), SCFSettings(); method=:dense, runtime=RuntimeSettings(backend=:cuda),
    )
    @test !preflight_dense_cuda.runnable
    @test any(issue -> issue.code == :runtime_unavailable, preflight_dense_cuda.issues)
    @test_throws ArgumentError preflight(
        ChainModel(size=4, hopping=-1.0, interaction=0.0), method=:unknown,
    )
end
