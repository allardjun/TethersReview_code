using Random
using LinearAlgebra
using CairoMakie


"""
    true_brownian_bridge_2d(N, δ, target; start=[0.0, 0.0])

Generate a true discrete Brownian Bridge using the conditional distribution approach.
This has the correct covariance structure without free parameters.
"""
function true_brownian_bridge_2d(N::Int, δ::Real, target::Vector{<:Real}; 
                                start::Vector{<:Real}=[0.0, 0.0])
    
    path = zeros(N+1, 2)
    path[1, :] = start
    
    # Total distance and time
    total_displacement = target - start
    total_time = N * δ^2  # "time" in Brownian motion scales as δ²
    
    for i in 1:N
        # Current time and remaining time
        current_time = i * δ^2
        remaining_time = total_time - current_time
        
        # Mean of conditional distribution (linear interpolation)
        mean_pos = start + (current_time / total_time) * total_displacement
        
        # Variance of conditional distribution
        # For Brownian Bridge: Var = current_time * remaining_time / total_time
        variance = current_time * remaining_time / total_time
        
        # Sample from conditional distribution
        if variance > 0
            noise = sqrt(variance) * randn(2)
            path[i+1, :] = mean_pos + noise
        else
            path[i+1, :] = mean_pos
        end
    end
    
    # Ensure exact endpoint (due to numerical precision)
    path[end, :] = target
    
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
    axislegend(ax)
    
    return fig
end

# Example usage
function example_usage()
    println("Simulating 2D Brownian Bridge...")
    
    # Parameters
    N = 100          # number of steps
    δ = 0.1          # step size
    target = [5.0, 3.0]  # target location
    noise_strength = 0.8
    
    # Generate the path
    Random.seed!(42)  # for reproducibility
    path = true_brownian_bridge_2d(N, δ, target)
    
    # Print some statistics
    println("Number of steps: $N")
    println("Step size: $δ")
    println("Target: $target")
    println("Final position: $(path[end, :])")
    println("Distance from target: $(norm(path[end, :] - target))")
    
    # Plot the result
    fig = plot_brownian_bridge(path)
    display(fig)
    
    return path
end


# Run the example
example_usage()
