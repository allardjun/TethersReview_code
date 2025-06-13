# 2D Brownian Bridge Simulation with Constrained Displacement

This project simulates 2D Brownian bridges with displacement projection constraints for molecular tether analysis. The simulation generates trajectories where each step has a fixed displacement magnitude (representing physical constraints like tether length) while maintaining the statistical properties of a Brownian bridge.

## Overview

This Julia project implements a specialized Brownian bridge simulation that:

1. **Constrains displacement magnitude**: Each step in the trajectory has exactly the prescribed displacement magnitude `delta` (representing physical tether length)
2. **Maintains bridge properties**: The trajectory still exhibits Brownian bridge behavior, attempting to reach a target endpoint
3. **Generates publication figures**: Batch processes parameter sets to create figure panels for research publications

## Core Algorithm

The simulation uses a modified Brownian bridge approach:

- **Brownian Bridge Drift**: At each step, applies a drift term `(x_end - x_current) / (T_remaining)` to bias the trajectory toward the target
- **Displacement Projection**: Projects each step to have exactly magnitude `delta`, preserving the direction but enforcing the constraint
- **Diffusion Coefficient**: Uses `D = delta²/6` based on the relationship between step size and diffusion

## Files

### `discrete_bridge3.jl`
Core simulation engine containing:
- `brownian_bridge_2d_projected()`: Main simulation function
- `plot_trajectory()`: Visualization function
- Example simulation with default parameters

### `batch_brownian_bridge.jl`
Batch processing system that:
- Reads simulation parameters from CSV files
- Processes multiple parameter sets automatically
- Generates individual PDF figures for each parameter set
- Uses Figure/Panel naming convention for output files

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
julia> include("batch_brownian_bridge.jl")
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
- **N**: Number of trajectory segments
- **T_end**: Total simulation time (automatically set to `delta * N`)
- **D**: Diffusion coefficient (automatically calculated as `delta²/6`)

## Output

Each simulation generates:
1. **PDF plot** showing the complete trajectory with:
   - Blue trajectory line
   - Green start point marker
   - Red target endpoint marker
   - Purple actual endpoint marker
   - Orange projection points
2. **Console statistics** including displacement magnitude verification and final distance from target

## Applications

This simulation is designed for analyzing molecular tethers and constrained random walks in biological systems, where the physical constraint of tether length must be maintained while studying the statistical properties of the resulting trajectories.