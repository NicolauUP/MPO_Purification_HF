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
