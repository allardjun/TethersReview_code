using LinearAlgebra
using CairoMakie
using Random
using Statistics

function discrete_brownian_bridge_2d(x_end, delta, N; seed=123)
    """
    Simulate a 2D discrete Brownian Bridge with fixed step size.
    Uses a bias in the random walk direction based on bridge conditioning.
    
    Parameters:
    - x_end: target ending position [x, y]
    - delta: step size
    - N: number of steps
    - seed: random seed
    
    Returns:
    - times: time points
    - trajectory: 2×(N+1) array of positions
    """
    Random.seed!(seed)
    
    x_start = [0.0, 0.0]
    
    # Initialize arrays
    times = range(0, N, length=N+1)
    trajectory = zeros(2, N+1)
    trajectory[:, 1] = x_start
    
    for i in 2:N+1
        x_current = trajectory[:, i-1]
        steps_remaining = N - (i - 1)
        
        if steps_remaining == 0
            # Last step: go directly to target
            trajectory[:, i] = x_end
        else
            # Vector pointing toward the target
            to_target = x_end - x_current
            distance_to_target = norm(to_target)
            
            if distance_to_target < 1e-10
                # Already at target, pure random walk
                angle = 2π * rand()
                direction = [cos(angle), sin(angle)]
            else
                # Compute drift direction
                drift_direction = to_target / distance_to_target
                
                # Compute drift strength based on how many steps we have left
                # and how far we need to go
                drift_strength = distance_to_target / (delta * steps_remaining)
                
                # Clamp drift strength to prevent instability
                drift_strength = clamp(drift_strength, 0.0, 0.9)
                
                # Generate random direction
                angle = 2π * rand()
                random_direction = [cos(angle), sin(angle)]
                
                # Combine drift and random components
                direction = drift_strength * drift_direction + 
                           sqrt(1 - drift_strength^2) * random_direction
                direction = direction / norm(direction)
            end
            
            # Take step of size delta
            trajectory[:, i] = x_current + delta * direction
        end
    end
    
    # Adjust the last few steps to ensure we exactly hit the target
    # This is a simple correction that redistributes the error
    if norm(trajectory[:, N+1] - x_end) > 1e-10
        error_vec = x_end - trajectory[:, N+1]
        # Distribute error over last few steps
        correction_steps = min(5, N)
        for j in (N+1-correction_steps+1):N+1
            weight = (j - (N+1-correction_steps)) / correction_steps
            trajectory[:, j] += weight * error_vec
        end
        trajectory[:, N+1] = x_end
    end
    
    return times, trajectory
end

function plot_trajectory(trajectory, x_end, delta; axis_limits=nothing, trajectory_only=false)
    """
    Plot the 2D trajectory with line segments.
    
    Parameters:
    - trajectory: 2×(N+1) array of positions
    - x_end: target ending position [x, y]
    - delta: displacement magnitude
    - axis_limits: optional uniform axis limits (default: auto-scale)
    - trajectory_only: if true, plot only the trajectory line without markers (default: false)
    """
    fig = Figure(size = (900, 700))
    ax = Makie.Axis(fig[1, 1], 
              xlabel = "x", 
              ylabel = "y", 
              title = "2D Brownian Bridge with Delta Projection\n(l_K = $delta)",
              aspect = DataAspect())
    
    # Plot trajectory as connected line segments
    lines!(ax, trajectory[1, :], trajectory[2, :], 
           color = :black, linewidth = 2, label = "Trajectory")
    
    # Add markers only if not trajectory_only mode
    if !trajectory_only
        x_start = trajectory[:, 1]

        # Mark start and end points
        scatter!(ax, [x_start[1]], [x_start[2]], 
                 color = :green, markersize = 20, 
                 marker = :circle, label = "Start")
        scatter!(ax, [x_end[1]], [x_end[2]], 
                 color = :red, markersize = 20, 
                 marker = :star5, label = "Target End")
        
        # Mark actual end point
        scatter!(ax, [trajectory[1, end]], [trajectory[2, end]], 
                 color = :purple, markersize = 15, 
                 marker = :diamond, label = "Actual End")
        
        # Mark all projection points (every delta_t)
        scatter!(ax, trajectory[1, :], trajectory[2, :], 
                 color = :orange, markersize = 8, 
                 alpha = 0.7, label = "Projection Points")
    end
    
    # Set axis limits if provided
    if axis_limits !== nothing
        xlims!(ax, -axis_limits, axis_limits)
        ylims!(ax, -axis_limits, axis_limits)
    end
    
    # Add grid and legend (only show legend if not trajectory_only)
    # ax.grid = true
    if !trajectory_only
        axislegend(ax, position = :lt)
    end
    
    # Display displacement magnitudes for verification
    displacements = [norm(trajectory[:, i] - trajectory[:, i-1]) for i in 2:size(trajectory, 2)]
    println("Displacement magnitudes (should all be ≈ $delta):")
    println("Mean: $(round(mean(displacements), digits=4))")
    println("Std:  $(round(std(displacements), digits=6))")
    println("Min:  $(round(minimum(displacements), digits=6))")
    println("Max:  $(round(maximum(displacements), digits=6))")
    
    return fig
end

# # Example usage and simulation
# begin
#     # Parameters
#     x_start = [0.0, 0.0]      # Starting position
#     x_end = [1.2, 0.1]        # Target ending position  
#     N = 60                    # Number of time segments
#     T_end = N              # Total time
#     delta = 0.1               # Prescribed displacement magnitude
    
#     println("Simulation Parameters:")
#     # println("- Start: $x_start")
#     println("- Target End: $x_end") 
#     println("- Delta: $delta")
#     println("- Number of segments: $N")
#     println("- Diffusion coefficient D = l_K²/6 = $(delta^2/6)")
#     println()
    
#     # Run simulation
#     times, trajectory = brownian_bridge_2d_projected(x_end, delta, N)
    
#     # Create plot
#     fig = plot_trajectory(trajectory, x_end, delta)
    
#     # Display the figure
#     # display(fig)
    
#     # Save the figure to a file
#     save("brownian_bridge_2d_projected.pdf", fig)

#     # Calculate final distance from target
#     final_distance = norm(trajectory[:, end] - x_end)
#     println("\nFinal distance from target: $(round(final_distance, digits=4))")
# end