#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 OUTPUT_DIRECTORY [extraction_scaling_2d.jl options]" >&2
  exit 2
fi

output_dir=$1
shift
mkdir -p "$output_dir"

# One Julia process makes process_time.txt a whole-run peak-RSS measurement.
if /usr/bin/time -v true >/dev/null 2>&1; then
  /usr/bin/time -v julia --startup-file=no --project=. \
    benchmark/extraction_scaling_2d/extraction_scaling_2d.jl --output "$output_dir" "$@" \
    2> >(tee "$output_dir/process_time.txt" >&2)
else
  /usr/bin/time -l julia --startup-file=no --project=. \
    benchmark/extraction_scaling_2d/extraction_scaling_2d.jl --output "$output_dir" "$@" \
    2> >(tee "$output_dir/process_time.txt" >&2)
fi
