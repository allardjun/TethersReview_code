using Random
using LinearAlgebra
using CairoMakie

"""
    brownian_bridge_2d(N, l_K, target; noise_strength=1.0, start=[0.0, 0.0])

Simulate a 2D random walk that ends at a specified target location.

# Arguments
- `N`: number of steps
- `l_K`: Kuhn length (step size)
- `target`: target endpoint as [x, y]
- `noise_strength`: controls randomness (0=deterministic, 1=balanced)
- `start`: starting position (default: origin)

# Returns
- `path`: matrix of size (N+1, 2) containing the full path
"""
function brownian_bridge_2d(N::Int, l_K::Real, target::Vector{<:Real}; 
                           noise_strength::Real=1.0, start::Vector{<:Real}=[0.0, 0.0])
    
    # Initialize path array
    path = zeros(N+1, 2)
    current_pos = copy(start)
    path[1, :] = current_pos
    
    for i in 1:N
        # Calculate remaining displacement needed
        remaining_steps = N - i + 1
        remaining_displacement = target - current_pos
        
        # Ideal step direction (deterministic component)
        ideal_step = remaining_displacement / remaining_steps
        
        # Add random component
        random_angle = 2*pi * rand()
        random_step = [cos(random_angle), sin(random_angle)]
        
        # Combine ideal and random components
        combined_direction = ideal_step + noise_strength * random_step
        
        # Normalize to unit vector and scale by l_K
        if norm(combined_direction) > 0
            step_direction = combined_direction / norm(combined_direction)
        else
            step_direction = random_step
        end
        
        # Take the step
        current_pos += l_K * step_direction
        path[i+1, :] = current_pos
    end
    
    return path
end

"""
    plot_brownian_bridge(path; title="2D Brownian Bridge")

Plot a 2D Brownian bridge path.
"""
function plot_brownian_bridge(path::Matrix{<:Real}; title::String="2D Brownian Bridge")
    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1, 1], 
              title=title,
              xlabel="X",
              ylabel="Y",
              aspect=DataAspect())
    
    # Plot the path
    lines!(ax, path[:, 1], path[:, 2], 
           linewidth=2, 
           alpha=0.7, 
           color=:blue,
           label="Path")
    
    # Mark start and end points
    scatter!(ax, [path[1, 1]], [path[1, 2]], 
             color=:green, 
             markersize=15, 
             label="Start")
    
    scatter!(ax, [path[end, 1]], [path[end, 2]], 
             color=:red, 
             markersize=15, 
             label="Target")
    
    # Add legend
    # axislegend(ax)
    
    return fig
end

# Example usage
function example_usage()
    println("Simulating 2D Brownian Bridge...")
    
    # Parameters
    N = 100          # number of steps
    l_K = 0.1          # step size
    target = [5.0, 3.0]  # target location
    noise_strength = 0.8
    
    # Generate the path
    # Random.seed!(42)  # for reproducibility
    path = brownian_bridge_2d(N, l_K, target; noise_strength=noise_strength)
    
    # Print some statistics
    println("Number of steps: $N")
    println("Step size: $l_K")
    println("Target: $target")
    println("Final position: $(path[end, :])")
    println("Distance from target: $(norm(path[end, :] - target))")
    
    # Plot the result
    fig = plot_brownian_bridge(path)
    display(fig)
    
    return path
end


example_usage()
