#!/bin/bash

# Environment variables to tell Julia about our nonstandard install location
# Change these if you change users:
export JULIA_DIR="/dartfs-hpc/scratch/alex/julia-1.7.3"
export JULIA_DEPOT_PATH="/dartfs-hpc/scratch/alex/.julia"
# Probably don't have to change these:
export JULIA_PROJECT="$JULIA_DEPOT_PATH/environments/v1.7/Project.toml"
export JULIA_LOAD_PATH="$JULIA_DEPOT_PATH/environments/v1.7/Project.toml:$JULIA_DIR/share/julia/stdlib/v1.7"

# Run Julia!
$JULIA_DIR/bin/julia --project=$JULIA_PROJECT ./strongscaling.jl