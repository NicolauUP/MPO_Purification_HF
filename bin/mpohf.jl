#!/usr/bin/env julia

"""Generic public campaign entry point for MPO_MeanField.

Usage:
    julia --project=. bin/mpohf.jl preflight CAMPAIGN.jl --task N --method mpo|dense|kpm [--backend cpu|cuda]
    julia --project=. bin/mpohf.jl run CAMPAIGN.jl --task N --method mpo|dense|kpm --output-root DIR [--backend cpu|cuda] [--verbose nothing|all] [--defer-observables true|false]

Campaign files are Julia source and must bind `campaign` to a `CampaignSpec`.
CUDA is loaded before `MPO_MeanField` when `--backend cuda` is requested, then
preflighted explicitly. The current MPO CUDA path is hybrid, with
host-side Hartree/Fock extraction, while the public KPM path is available for
open squares, rectangles, and explicit sparse graphs.
"""

# CUDA's ITensors/NDTensors extensions must exist before the solver package is
# loaded. Loading CUDA lazily from inside `solve` leaves already-compiled calls
# in an older Julia world and can silently select CPU-bound generic work for
# GPU-backed tensors. This lightweight raw-argument scan precedes the strict
# parser below; malformed options are still rejected by `parse_options`.
const _MPOHF_CUDA_REQUESTED = any(
    index -> ARGS[index] == "--backend" && ARGS[index + 1] == "cuda",
    1:max(0, length(ARGS) - 1),
)
if _MPOHF_CUDA_REQUESTED
    using CUDA
end
using MPO_MeanField

function usage(io::IO=stderr)
    println(io, "usage:")
    println(io, "  mpohf.jl preflight CAMPAIGN.jl --task N --method mpo|dense|kpm [--backend cpu|cuda]")
    println(io, "  mpohf.jl run CAMPAIGN.jl --task N --method mpo|dense|kpm --output-root DIR [--backend cpu|cuda] [--verbose nothing|all] [--defer-observables true|false]")
end

function parse_options(arguments::Vector{String})
    length(arguments) >= 2 || throw(ArgumentError("missing command or campaign source"))
    command = arguments[1]
    command in ("preflight", "run") || throw(ArgumentError("command must be `preflight` or `run`"))
    campaign_path = abspath(arguments[2])
    options = Dict{String,String}()
    index = 3
    while index <= length(arguments)
        option = arguments[index]
        startswith(option, "--") || throw(ArgumentError("expected option beginning with --, got $option"))
        index < length(arguments) || throw(ArgumentError("missing value for $option"))
        haskey(options, option) && throw(ArgumentError("option supplied more than once: $option"))
        options[option] = arguments[index + 1]
        index += 2
    end
    allowed = command == "run" ?
        Set(["--task", "--method", "--backend", "--output-root", "--verbose", "--defer-observables"]) :
        Set(["--task", "--method", "--backend"])
    all(key -> key in allowed, keys(options)) || throw(ArgumentError(
        "unsupported option(s): $(join(sort!(collect(setdiff(Set(keys(options)), allowed))), ", "))",
    ))
    haskey(options, "--task") || throw(ArgumentError("--task N is required"))
    haskey(options, "--method") || throw(ArgumentError("--method mpo|dense|kpm is required"))
    command == "run" && !haskey(options, "--output-root") && throw(ArgumentError(
        "--output-root DIR is required for run",
    ))
    task = tryparse(Int, options["--task"])
    isnothing(task) && throw(ArgumentError("--task must be an integer"))
    task > 0 || throw(ArgumentError("--task must be positive"))
    method = Symbol(options["--method"])
    method in (:mpo, :dense, :kpm) || throw(ArgumentError("--method must be mpo, dense, or kpm"))
    backend = Symbol(get(options, "--backend", "cpu"))
    backend in (:cpu, :cuda) || throw(ArgumentError("--backend must be cpu or cuda"))
    verbose = Symbol(get(options, "--verbose", "nothing"))
    verbose in (:nothing, :all) || throw(ArgumentError("--verbose must be nothing or all"))
    defer_text = get(options, "--defer-observables", "false")
    defer_text in ("true", "false") || throw(ArgumentError(
        "--defer-observables must be true or false",
    ))
    defer_observables = defer_text == "true"
    defer_observables && method != :mpo && throw(ArgumentError(
        "--defer-observables=true is currently supported only by --method mpo",
    ))
    return (; command, campaign_path, task, method, backend, verbose,
        defer_observables,
        output_root=get(options, "--output-root", nothing))
end

function load_campaign(path::AbstractString)
    isfile(path) || throw(ArgumentError("campaign source does not exist: $path"))
    namespace = Module(gensym(:MPOHFPublicCampaign))
    Core.eval(namespace, :(using MPO_MeanField))
    Base.include(namespace, path)
    Base.invokelatest(isdefined, namespace, :campaign) || throw(ArgumentError(
        "campaign source must bind `campaign` to a CampaignSpec: $path",
    ))
    # `include` creates bindings at a newer Julia world age than this helper.
    campaign = Base.invokelatest(getfield, namespace, :campaign)
    campaign isa CampaignSpec || throw(ArgumentError(
        "campaign must be a CampaignSpec, got $(typeof(campaign)); migrate the source instead of passing a legacy campaign",
    ))
    return campaign
end

function main(arguments::Vector{String})
    parsed = parse_options(arguments)
    campaign = load_campaign(parsed.campaign_path)
    runtime = RuntimeSettings(backend=parsed.backend)
    if parsed.command == "preflight"
        report = preflight(campaign; task=parsed.task, method=parsed.method, runtime=runtime)
        println(report)
        return report.runnable ? 0 : 2
    end
    # The campaign source is evaluated in a fresh module so it can contain
    # ordinary Julia hopping/seed functions. Dispatching through
    # `invokelatest` is required on Julia 1.12: those closures are newer than
    # this CLI method's world age.
    directory, result = Base.invokelatest(
        run_campaign, campaign;
        task=parsed.task,
        method=parsed.method,
        output_root=parsed.output_root,
        verbose=parsed.verbose,
        runtime=runtime,
        source_path=parsed.campaign_path,
        defer_observables=parsed.defer_observables,
    )
    println("Result directory: $directory")
    println("SCF: converged=$(result.converged) termination=$(result.termination_reason)")
    return result.converged ? 0 : 2
end

try
    exit(main(copy(ARGS)))
catch error
    if get(ENV, "MPO_DEBUG_BACKTRACE", "0") == "1"
        println(stderr, "mpohf.jl: ", sprint(showerror, error, catch_backtrace()))
    else
        println(stderr, "mpohf.jl: ", sprint(showerror, error))
    end
    usage(stderr)
    exit(1)
end
