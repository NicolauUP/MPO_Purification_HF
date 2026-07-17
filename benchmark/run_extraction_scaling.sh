#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 OUTPUT_DIRECTORY [extraction_scaling.jl options]" >&2
  exit 2
fi

output_dir=$1
shift
mkdir -p "$output_dir"

# This wrapper intentionally launches exactly one Julia process. The resulting
# process_time.txt is a process-wide peak-RSS measurement, not a per-method
# measurement; the per-method allocations/times are in samples.csv.
if /usr/bin/time -v true >/dev/null 2>&1; then
  /usr/bin/time -v julia --project=. benchmark/extraction_scaling.jl \
    --output "$output_dir" "$@" 2>"$output_dir/process_time.txt"
else
  /usr/bin/time -l julia --project=. benchmark/extraction_scaling.jl \
    --output "$output_dir" "$@" 2>"$output_dir/process_time.txt"
fi
