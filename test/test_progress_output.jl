@testset "terminal iteration progress" begin
    terminal = IOBuffer()
    overwriting = MPO_MeanField.print_iteration_progress(
        terminal, "Purification", 35, 60, "Tr=2"; overwrite=true,
    )
    MPO_MeanField.print_iteration_progress(
        terminal, "Purification", 36, 60, "Tr=2"; overwrite=true,
    )
    MPO_MeanField.finish_iteration_progress(terminal, overwriting)
    @test String(take!(terminal)) == "\r\e[2KPurification 35/60 | Tr=2\r\e[2KPurification 36/60 | Tr=2\n"

    log = IOBuffer()
    overwriting = MPO_MeanField.print_iteration_progress(
        log, "Purification", 35, 60, "Tr=2"; overwrite=false,
    )
    MPO_MeanField.finish_iteration_progress(log, overwriting)
    @test String(take!(log)) == "Purification 35/60 | Tr=2\n"
end

@testset "SP2 redirected progress is complete" begin
    sys, params, _, bounds, rho0 = sp2_static_system()
    log = IOBuffer()
    result = perform_purification(
        rho0, params;
        method=:sp2,
        verbose=1,
        io=log,
        overwrite_progress=false,
        spectral_bounds=bounds,
    )
    output = String(take!(log))
    @test result.converged
    @test occursin("SP2 purifying", output)
    @test occursin("SP2 1/", output)
    @test all(line -> !occursin('\r', line), split(chomp(output), '\n'))
end

@testset "SCF progress is terminal-aware" begin
    params = parameters_1d(
        L=2,
        t=-0.7,
        U=0.0,
        W=x -> (0.0, 0.13, -0.08, 0.21)[Int(x) + 1],
        S=nothing,
        purification_steps=35,
        scf_max_iterations=5,
    )
    redirected = IOBuffer()
    redirected_sys = System(params)
    @test run_scf!(redirected_sys, -5.0, 5.0;
        purification_method=:sp2,
        verbose=:all,
        io=redirected,
        overwrite_progress=false,
    )
    redirected_output = String(take!(redirected))
    @test occursin("SCF 1/", redirected_output)
    @test occursin("Tr=", redirected_output)
    @test all(line -> !occursin('\r', line), split(chomp(redirected_output), '\n'))

    terminal = IOBuffer()
    terminal_sys = System(params)
    @test run_scf!(terminal_sys, -5.0, 5.0;
        purification_method=:sp2,
        verbose=:all,
        io=terminal,
        overwrite_progress=true,
    )
    terminal_output = String(take!(terminal))
    @test count("\r\e[2KSCF", terminal_output) >= 2
    @test endswith(terminal_output, "SCF converged in 3 iterations.\n")
end
