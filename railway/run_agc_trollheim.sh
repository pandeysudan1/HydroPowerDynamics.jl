#!/usr/bin/env bash
set -euo pipefail
cd /workspace/HydroPowerDynamics.jl
mkdir -p /workspace/results

echo "=== HydroPowerDynamics.jl Trollheim AGC comparison ==="
julia --version
julia --project=. quick_examples/trollheim_agc_surgetank_compare/trollheim_agc_surgetank_compare.jl \
  2>&1 | tee /workspace/results/agc_trollheim_hpd.log

echo "=== GENERATED PLOTS ==="
find quick_examples/trollheim_agc_surgetank_compare/plots -maxdepth 1 -type f -print || true
