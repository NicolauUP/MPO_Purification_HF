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
