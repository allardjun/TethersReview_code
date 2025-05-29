using LinearAlgebra
using CairoMakie
using Random
using Statistics

function brownian_bridge_2d_projected(x_start, x_end, T_end, delta, N; seed=123)
    """
    Simulate a 2D Brownian Bridge with displacement projection.
    
    Parameters:
    - x_start: starting position [x, y]
    - x_end: target ending position [x, y] 
    - T_end: total time
    - delta: prescribed displacement magnitude for each segment
    - N: number of time segments
    - seed: random seed for reproducibility
    
    Returns:
    - times: time points
    - trajectory: 2×(N+1) array of positions
    """
    # Random.seed!(seed)
    
    # Parameters
    D = delta^2 / 6           # Diffusion coefficient
    delta_t = T_end / N       # Time step duration
    
    # Initialize arrays
    times = range(0, T_end, length=N+1)
    trajectory = zeros(2, N+1)
    trajectory[:, 1] = x_start
    
    # Track the last projection point for displacement calculation
    last_projection_point = copy(x_start)
    
    for i in 2:N+1
        t_prev = times[i-1]
        x_prev = trajectory[:, i-1]
        
        # Brownian bridge drift term: (x_end - x) / (T_end - t)
        remaining_time = T_end - t_prev
        if remaining_time > 1e-10
            drift = (x_end - x_prev) / remaining_time
        else
            drift = [0.0, 0.0]
        end
        
        # Generate 2D Wiener increment: sqrt(delta_t) * randn(2)
        dW = sqrt(delta_t) * randn(2)
        
        # Brownian bridge step: dx = drift*dt + sqrt(2D)*dW
        dx = drift * delta_t + sqrt(2*D) * dW
        x_candidate = x_prev + dx
        
        # Project displacement from last projection point to have magnitude exactly delta
        displacement = x_candidate - last_projection_point
        displacement_magnitude = norm(displacement)
        
        if displacement_magnitude > 1e-12
            # Scale displacement to have magnitude exactly delta
            x_new = last_projection_point + delta * (displacement / displacement_magnitude)
        else
            # If displacement is essentially zero, choose random direction
            random_direction = randn(2)
            random_direction = random_direction / norm(random_direction)
            x_new = last_projection_point + delta * random_direction
        end
        
        trajectory[:, i] = x_new
        last_projection_point = copy(x_new)
    end
    
    return times, trajectory
end

function plot_trajectory(times, trajectory, x_start, x_end, delta)
    """
    Plot the 2D trajectory with line segments.
    """
    fig = Figure(resolution = (900, 700))
    ax = Makie.Axis(fig[1, 1], 
              xlabel = "x", 
              ylabel = "y", 
              title = "2D Brownian Bridge with Delta Projection\n(l_K = $delta)",
              aspect = DataAspect())
    
    # Plot trajectory as connected line segments
    lines!(ax, trajectory[1, :], trajectory[2, :], 
           color = :blue, linewidth = 2, label = "Trajectory")
    
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
    
    # Add grid and legend
    # ax.grid = true
    axislegend(ax, position = :lt)
    
    # Display displacement magnitudes for verification
    displacements = [norm(trajectory[:, i] - trajectory[:, i-1]) for i in 2:size(trajectory, 2)]
    println("Displacement magnitudes (should all be ≈ $delta):")
    println("Mean: $(round(mean(displacements), digits=4))")
    println("Std:  $(round(std(displacements), digits=6))")
    println("Min:  $(round(minimum(displacements), digits=6))")
    println("Max:  $(round(maximum(displacements), digits=6))")
    
    return fig
end

# Example usage and simulation
begin
    # Parameters
    x_start = [0.0, 0.0]      # Starting position
    x_end = [5.0, 3.0]        # Target ending position  
    N = 100                    # Number of time segments
    T_end = N              # Total time
    delta = 0.5               # Prescribed displacement magnitude
    
    println("Simulation Parameters:")
    println("- Start: $x_start")
    println("- Target End: $x_end") 
    println("- Total Time: $T_end")
    println("- Delta: $delta")
    println("- Number of segments: $N")
    println("- Diffusion coefficient D = l_K²/6 = $(delta^2/6)")
    println()
    
    # Run simulation
    times, trajectory = brownian_bridge_2d_projected(x_start, x_end, T_end, delta, N)
    
    # Create plot
    fig = plot_trajectory(times, trajectory, x_start, x_end, delta)
    
    # Display the figure
    display(fig)
    
    # Calculate final distance from target
    final_distance = norm(trajectory[:, end] - x_end)
    println("\nFinal distance from target: $(round(final_distance, digits=4))")
end