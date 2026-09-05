@testset "Atomic MPO checkpoint and deferred observables" begin
    model = ChainModel(size=8, hopping=-0.7, interaction=0.0, filling=0.5)
    representation = QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64)
    solver = SCFSettings(
        purification=:sp2,
        mixing=0.5,
        tolerance=0.1,
        maxiter=5,
        purification_maxiter=35,
        record_energy=false,
    )
    params = legacy_parameters(model, representation, solver)

    mktempdir() do directory
        checkpoint = joinpath(directory, "converged_state.h5")
        result = solve(
            model, representation, solver;
            method=:mpo,
            spectral_bounds=(-3.0, 3.0),
            checkpoint_path=checkpoint,
            measure_observables=false,
        )
        @test result.converged
        @test result.observables === nothing
        @test isfile(checkpoint)

        loaded = read_mpo_checkpoint(checkpoint, params)
        @test loaded.converged
        @test loaded.termination_reason == :converged
        @test isapprox(real(tr(loaded.system.ρ)), 4.0; atol=3e-3)
        measured = observables(loaded.system)
        @test isapprox(measured.particle_number, 4.0; atol=3e-3)
        @test isfinite(measured.energy.total)

        incompatible = Parameters1D(
            L=params.L,
            t=params.t,
            U=0.1,
            W=params.W,
            S=params.S,
            tci_tol=params.tci_tol,
            itensors_tol=params.itensors_tol,
            itensors_maxdim=params.itensors_maxdim,
            density=params.density,
            purification_steps=params.purification_steps,
            scf_mixing=params.scf_mixing,
            scf_tol=params.scf_tol,
            scf_max_iterations=params.scf_max_iterations,
        )
        @test_throws ArgumentError read_mpo_checkpoint(checkpoint, incompatible)
        @test_throws ArgumentError write_mpo_checkpoint(loaded.system, checkpoint)
    end
end
