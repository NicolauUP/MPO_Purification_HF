#!/usr/bin/env julia

"""Collect deterministic MPO-validation data for the separable 2D campaign.

Usage
-----
    julia --project=. runs/common/collect_separable_mpo_report.jl MANIFEST OUTPUT

The manifest only names *existing* result directories. This collector never
changes those results: it normalizes their small summary and trajectory files
into a presentation-friendly report directory. It deliberately has no plotting
dependency; `notebooks/analyze_separable_mpo_report.ipynb` consumes the output.
"""

using Dates
using Printf
using TOML

length(ARGS) == 2 || error(
    "usage: collect_separable_mpo_report.jl MANIFEST.toml OUTPUT_DIRECTORY",
)

manifest_path = abspath(ARGS[1])
output = abspath(ARGS[2])
isfile(manifest_path) || error("manifest does not exist: $manifest_path")
ispath(output) && error("refusing to overwrite existing output directory: $output")
manifest = TOML.parsefile(manifest_path)
inputs = get(manifest, "inputs", Dict{String,Any}())

directories(key) = begin
    value = get(inputs, key, String[])
    value isa Vector || error("inputs.$key must be an array of directories")
    paths = abspath.(String.(value))
    all(isdir, paths) || error("one or more inputs.$key directories do not exist")
    paths
end

compression_dirs = directories("compression_directories")
sp2_dirs = directories("sp2_directories")
dense_dirs = directories("dense_directories")
mpo_scf_dirs = directories("mpo_scf_directories")

function csv_fields(line::AbstractString)
    fields = String[]
    buffer = IOBuffer()
    quoted = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            next = nextind(line, i)
            if quoted && next <= lastindex(line) && line[next] == '"'
                print(buffer, '"')
                i = nextind(line, next)
                continue
            end
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(fields, String(take!(buffer)))
        elseif c != '\r' && c != '\n'
            print(buffer, c)
        end
        i = nextind(line, i)
    end
    quoted && error("unterminated quoted CSV field: $line")
    push!(fields, String(take!(buffer)))
    return fields
end

function read_csv(path)
    open(path) do io
        eof(io) && error("empty CSV file: $path")
        header = csv_fields(readline(io))
        rows = Vector{Dict{String,String}}()
        for line in eachline(io)
            isempty(strip(line)) && continue
            fields = csv_fields(line)
            length(fields) == length(header) || error("malformed row in $path")
            push!(rows, Dict(zip(header, fields)))
        end
        return rows
    end
end

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(csv_escape.(header), ','))
        for row in rows
            println(io, join(csv_escape.([get(row, column, "") for column in header]), ','))
        end
    end
end

number(value) = tryparse(Float64, value)
number_or_blank(row, key) = begin
    value = get(row, key, "")
    parsed = number(value)
    isnothing(parsed) ? "" : @sprintf("%.16g", parsed)
end
toml_or_blank(table, key) = haskey(table, key) ? string(table[key]) : ""
label_from(directory, table) = begin
    candidate = get(table, "label", basename(directory))
    string(candidate)
end

mkpath(output)

# Exact-projector compression ladder: one `summary.csv` per supplied directory.
compression_rows = Dict{String,String}[]
for directory in compression_dirs
    path = joinpath(directory, "summary.csv")
    isfile(path) || error("compression input lacks summary.csv: $directory")
    metadata_path = joinpath(directory, "metadata.toml")
    metadata = isfile(metadata_path) ? TOML.parsefile(metadata_path) : Dict{String,Any}()
    for row in read_csv(path)
        push!(compression_rows, Dict(
            "label" => label_from(directory, metadata),
            "directory" => directory,
            "cutoff" => get(row, "cutoff", toml_or_blank(metadata, "cutoff")),
            "requested_maxdim" => get(row, "requested_maxdim", ""),
            "max_chi" => get(row, "max_chi", ""),
            "mean_chi" => get(row, "mean_chi", ""),
            "trace_error" => number_or_blank(row, "trace_error"),
            "density_max_abs_error" => number_or_blank(row, "density_max_abs_error"),
            "density_rms_error" => number_or_blank(row, "density_rms_error"),
            "horizontal_bond_rms_error" => number_or_blank(row, "horizontal_bond_rms_error"),
            "vertical_bond_rms_error" => number_or_blank(row, "vertical_bond_rms_error"),
            "one_body_energy_error" => number_or_blank(row, "one_body_energy_error"),
            "source" => "exact_projector_compression",
        ))
    end
end
write_csv(joinpath(output, "compression_ladder.csv"), [
    "label", "directory", "cutoff", "requested_maxdim", "max_chi", "mean_chi",
    "trace_error", "density_max_abs_error", "density_rms_error",
    "horizontal_bond_rms_error", "vertical_bond_rms_error", "one_body_energy_error", "source",
], compression_rows)

# Fixed-H SP2 summaries and trajectories. Supports both the compact comparison
# result (`summary.toml` + sp2_progress.txt) and the fixed-H diagnostic
# (`summary.toml` + iterations.csv).
sp2_summary_rows = Dict{String,String}[]
sp2_trajectory_rows = Dict{String,String}[]
for directory in sp2_dirs
    summary_path = joinpath(directory, "summary.toml")
    isfile(summary_path) || error("SP2 input lacks summary.toml: $directory")
    summary = TOML.parsefile(summary_path)
    label = label_from(directory, summary)
    requested_maxdim = get(summary, "itensors_maxdim", get(summary, "maxdim", ""))
    push!(sp2_summary_rows, Dict(
        "label" => label, "directory" => directory,
        "requested_maxdim" => string(requested_maxdim),
        "final_max_chi" => toml_or_blank(summary, "sp2_final_max_chi"),
        "work_max_chi" => toml_or_blank(summary, "sp2_work_max_chi"),
        "iterations" => toml_or_blank(summary, "sp2_iterations"),
        "converged" => toml_or_blank(summary, "sp2_converged"),
        "termination_reason" => toml_or_blank(summary, "sp2_termination_reason"),
        "trace_error" => toml_or_blank(summary, "sp2_trace_error"),
        "idempotency_residual" => toml_or_blank(summary, "sp2_idempotency_residual"),
        "density_max_abs_error" => toml_or_blank(summary, "density_max_abs_error"),
        "density_rms_error" => toml_or_blank(summary, "density_rms_error"),
        "horizontal_bond_rms_error" => toml_or_blank(summary, "horizontal_bond_rms_error"),
        "vertical_bond_rms_error" => toml_or_blank(summary, "vertical_bond_rms_error"),
        "one_body_energy_error" => toml_or_blank(summary, "one_body_energy_error"),
    ))

    iteration_path = joinpath(directory, "iterations.csv")
    if isfile(iteration_path)
        for row in read_csv(iteration_path)
            push!(sp2_trajectory_rows, Dict(
                "label" => label, "directory" => directory,
                "requested_maxdim" => string(requested_maxdim),
                "iteration" => get(row, "iteration", ""),
                "branch" => get(row, "branch", ""),
                "trace" => number_or_blank(row, "trace"),
                "trace_error" => number_or_blank(row, "trace_error"),
                "idempotency_residual" => number_or_blank(row, "idempotency_residual"),
                "hermiticity_residual" => number_or_blank(row, "hermiticity_residual"),
                "rho_max_chi" => get(row, "rho_max_chi", ""),
                "rho_mean_chi" => get(row, "rho_mean_chi", ""),
                "cap_reached" => get(row, "cap_reached", ""),
                "step_time_s" => number_or_blank(row, "step_time_s"),
            ))
        end
    else
        # `compare_square_sp2_dense_projector.jl` predates the structured
        # fixed-H diagnostic and logs its trajectory only as progress text.
        # Preserve the available quantities rather than silently omitting that
        # controlled trajectory from the P0 figures.
        progress_path = joinpath(directory, "sp2_progress.txt")
        if isfile(progress_path)
            pattern = r"SP2\s+(\d+)/\d+\s+\|\s+Tr=([^|]+)\|\s+err=([^|]+)\|\s+idem=([^|]+)\|\s+herm=([^|]+)\|\s+χ=\((\d+),(\d+)\)"
            for line in eachline(progress_path)
                matched = match(pattern, line)
                isnothing(matched) && continue
                capture = matched.captures
                rho_chi = parse(Int, capture[6])
                squared_chi = parse(Int, capture[7])
                push!(sp2_trajectory_rows, Dict(
                    "label" => label, "directory" => directory,
                    "requested_maxdim" => string(requested_maxdim),
                    "iteration" => capture[1], "branch" => "",
                    "trace" => strip(capture[2]), "trace_error" => strip(capture[3]),
                    "idempotency_residual" => strip(capture[4]),
                    "hermiticity_residual" => strip(capture[5]),
                    "rho_max_chi" => string(rho_chi), "rho_mean_chi" => "",
                    "cap_reached" => string(max(rho_chi, squared_chi) >= requested_maxdim),
                    "step_time_s" => "",
                ))
            end
        end
    end
end
write_csv(joinpath(output, "sp2_summary.csv"), [
    "label", "directory", "requested_maxdim", "final_max_chi", "work_max_chi", "iterations",
    "converged", "termination_reason", "trace_error", "idempotency_residual",
    "density_max_abs_error", "density_rms_error", "horizontal_bond_rms_error",
    "vertical_bond_rms_error", "one_body_energy_error",
], sp2_summary_rows)
write_csv(joinpath(output, "sp2_trajectory.csv"), [
    "label", "directory", "requested_maxdim", "iteration", "branch", "trace", "trace_error",
    "idempotency_residual", "hermiticity_residual", "rho_max_chi", "rho_mean_chi",
    "cap_reached", "step_time_s",
], sp2_trajectory_rows)

function collect_scf(directory, solver)
    metadata_path = joinpath(directory, "metadata.toml")
    observables_path = joinpath(directory, "observables.toml")
    isfile(metadata_path) || error("SCF input lacks metadata.toml: $directory")
    isfile(observables_path) || error("SCF input lacks observables.toml: $directory")
    metadata = TOML.parsefile(metadata_path)
    observables = TOML.parsefile(observables_path)
    label = label_from(directory, metadata)
    history_path = joinpath(directory, "scf_history.csv")
    trajectory = Dict{String,String}[]
    if isfile(history_path)
        for row in read_csv(history_path)
            push!(trajectory, Dict(
                "label" => label, "directory" => directory, "solver" => solver,
                "iteration" => get(row, "iteration", ""),
                "trace" => number_or_blank(row, "trace"),
                "vh_residual" => number_or_blank(row, "vh_residual"),
                "vf_residual" => number_or_blank(row, "vf_residual"),
                "rho_residual" => number_or_blank(row, "rho_residual"),
                "commutator_residual" => number_or_blank(row, "commutator_residual"),
                "two_cycle_residual" => number_or_blank(row, "two_cycle_residual"),
                "energy_total" => number_or_blank(row, "energy_total"),
                "rho_bond_dimension" => get(row, "rho_bond_dimension", ""),
            ))
        end
    end
    summary = Dict(
        "label" => label, "directory" => directory, "solver" => solver,
        "converged" => string(get(observables, "scf_converged", get(metadata, "scf_converged", ""))),
        "termination_reason" => string(get(observables, "scf_termination_reason", get(metadata, "scf_termination_reason", ""))),
        "iterations" => string(get(metadata, "scf_iterations", "")),
        "particle_number" => toml_or_blank(observables, "particle_number"),
        "energy_total" => toml_or_blank(observables, "energy_total"),
        "checkerboard_order" => toml_or_blank(observables, "checkerboard_density_order"),
        "stationarity_residual" => toml_or_blank(observables, "stationarity_residual"),
        "idempotency_residual" => toml_or_blank(observables, "idempotency_residual"),
        "rho_bond_dimension" => toml_or_blank(observables, "rho_bond_dimension"),
    )
    return summary, trajectory
end

scf_summary_rows = Dict{String,String}[]
scf_trajectory_rows = Dict{String,String}[]
for (solver, dirs) in (("dense", dense_dirs), ("mpo_scf", mpo_scf_dirs))
    for directory in dirs
        summary, trajectory = collect_scf(directory, solver)
        push!(scf_summary_rows, summary)
        append!(scf_trajectory_rows, trajectory)
    end
end
write_csv(joinpath(output, "scf_summary.csv"), [
    "label", "directory", "solver", "converged", "termination_reason", "iterations",
    "particle_number", "energy_total", "checkerboard_order", "stationarity_residual",
    "idempotency_residual", "rho_bond_dimension",
], scf_summary_rows)
write_csv(joinpath(output, "scf_trajectory.csv"), [
    "label", "directory", "solver", "iteration", "trace", "vh_residual", "vf_residual",
    "rho_residual", "commutator_residual", "two_cycle_residual", "energy_total",
    "rho_bond_dimension",
], scf_trajectory_rows)

report_metadata = Dict(
    "created_at_utc" => string(now(UTC)),
    "manifest" => manifest_path,
    "label" => get(get(manifest, "project", Dict{String,Any}()), "label", "separable_mpo_report"),
    "compression_directories" => compression_dirs,
    "sp2_directories" => sp2_dirs,
    "dense_directories" => dense_dirs,
    "mpo_scf_directories" => mpo_scf_dirs,
    "density_comparison_directory" => get(inputs, "density_comparison_directory", ""),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, report_metadata)
end

println("P0 report written: $output")
println("  compression points: $(length(compression_rows))")
println("  fixed-H SP2 runs:   $(length(sp2_summary_rows))")
println("  SP2 trajectory rows: $(length(sp2_trajectory_rows))")
println("  SCF runs:           $(length(scf_summary_rows))")
