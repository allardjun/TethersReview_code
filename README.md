# 2D Brownian Bridge Simulation with Constrained Displacement

This project simulates 2D Brownian bridges with displacement projection constraints for molecular tether analysis. The simulation generates trajectories where each step has a fixed displacement magnitude (representing physical constraints like tether length) while maintaining the statistical properties of a Brownian bridge.

## Overview

This Julia project implements a specialized Brownian bridge simulation that:

1. **Constrains displacement magnitude**: Each step in the trajectory has exactly the prescribed displacement magnitude `delta` (representing physical tether length)
2. **Maintains bridge properties**: The trajectory still exhibits Brownian bridge behavior, attempting to reach a target endpoint
3. **Generates publication figures**: Batch processes parameter sets to create figure panels for research publications

## Core Algorithm

The simulation uses a discrete Brownian bridge approach with adaptive biasing:

- **Adaptive Drift**: Computes drift strength as `distance_to_target / (delta * steps_remaining)`, providing stronger bias when far from target with few steps remaining
- **Direction Combination**: Blends drift direction toward target with random direction: `drift_strength * drift_direction + sqrt(1 - drift_strength²) * random_direction`
- **Stability Control**: Clamps drift strength to [0.0, 0.9] to prevent numerical instabilities
- **Fixed Step Size**: Every step has exactly magnitude `delta` by construction
- **Error Correction**: Redistributes final positioning error over the last few steps to ensure exact target hitting

## Files

### `discrete_bridge3.jl`
Core simulation engine containing:
- `discrete_brownian_bridge_2d()`: Main simulation function implementing adaptive biased random walk
- `plot_trajectory()`: Visualization function with support for trajectory-only mode
- Commented example simulation with default parameters

### `batch_brownian_bridge.jl`
Batch processing system that:
- Reads simulation parameters from CSV files
- Processes multiple parameter sets automatically
- Generates individual PDF figures for each parameter set
- Uses Figure/Panel naming convention for output files
- Supports uniform axis limits across all figures
- Optional trajectory-only mode for clean publication figures

### `sizes.csv`
Parameter file containing simulation configurations:
- **Figure/Panel**: Output filename identifiers
- **x_end_x, x_end_y**: Target endpoint coordinates
- **delta**: Displacement magnitude (tether length)
- **N**: Number of trajectory segments
- **seed**: Random seed (optional)

## Usage

### Single Simulation
```julia
julia> include("discrete_bridge3.jl")
# Runs example simulation and saves "brownian_bridge_2d_projected.pdf"
```

### Batch Processing
```julia
# Standard mode with all markers and legend
julia> batch_generate_brownian_bridges("sizes.csv")

# Clean trajectory-only mode (no markers or legend)
julia> batch_generate_brownian_bridges("sizes.csv"; trajectory_only=true)

# Processes all parameter sets in sizes.csv
# Generates individual PDFs: fig_1_A.pdf, fig_2_A.pdf, etc.
```

## Dependencies

- `LinearAlgebra`: Vector operations and norms
- `CairoMakie`: High-quality plotting and PDF generation
- `Random`: Random number generation and seeding
- `Statistics`: Statistical calculations
- `CSV`: Parameter file reading
- `DataFrames`: Data manipulation

## Parameters

- **x_end**: Target endpoint `[x, y]` coordinates
- **delta**: Fixed displacement magnitude per step (tether length)
- **N**: Number of trajectory steps
- **seed**: Random seed for reproducible results

## Output

Each simulation generates:
1. **PDF plot** showing the trajectory:
   - **Standard mode**: Black trajectory line with green start marker, red target marker, purple actual end marker, and orange projection points
   - **Trajectory-only mode**: Clean black trajectory line only (no markers or legend)
   - **Uniform scaling**: All batch-generated figures use identical axis limits for easy comparison
2. **Console statistics** including displacement magnitude verification and final distance from target

## Applications

This simulation is designed for analyzing molecular tethers and constrained random walks in biological systems. The improved discrete Brownian bridge model provides more accurate bridge conditioning while maintaining the physical constraint that each step has exactly the prescribed magnitude (representing tether length). This makes it particularly suitable for studying polymer dynamics, protein conformational changes, and other biophysical processes with length constraints.