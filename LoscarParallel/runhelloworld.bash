#!/bin/bash
# Environment variables to tell Julia about our nonstandard install location
export JULIA_BINDIR="/dartfs-hpc/scratch/alex/julia-1.5.3/bin"
export JULIA_DEPOT_PATH="/dartfs-hpc/scratch/alex/.julia"
export JULIA_PROJECT="/dartfs-hpc/scratch/alex/.julia/environments/v1.5/Project.toml"
export JULIA_LOAD_PATH="/dartfs-hpc/scratch/alex/.julia/environments/v1.5/Project.toml:/dartfs-hpc/scratch/alex/julia-1.5.3/share/julia/stdlib/v1.5"
# Run Julia!
$JULIA_BINDIR/julia --project=$JULIA_PROJECT ./helloworld.jl
