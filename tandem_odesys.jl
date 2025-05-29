using OrdinaryDiffEq
using CairoMakie
using Colors

# Define the ODE system (Extended Lotka-Volterra with additional dynamics)
function ode_system!(du, u, p, t)
    # Parameters
    r1, r2, K1, K2, a12, a21, b1, b2 = p
    
    # Variables: u[1] = prey1, u[2] = predator1, u[3] = prey2, u[4] = predator2
    x1, y1, x2, y2 = u
    
    # System of 4 coupled ODEs
    du[1] = r1 * x1 * (1 - (x1 + a12*x2)/K1) - b1 * x1 * y1  # prey1
    du[2] = -0.5 * y1 + 0.3 * b1 * x1 * y1                    # predator1
    du[3] = r2 * x2 * (1 - (x2 + a21*x1)/K2) - b2 * x2 * y2  # prey2
    du[4] = -0.4 * y2 + 0.25 * b2 * x2 * y2                   # predator2
end

# Parameters and initial conditions
tspan = (0.0, 50.0)
toff = 25.0  # Time when parameter change occurs
u0 = [10.0, 5.0, 8.0, 3.0]  # Initial conditions

# Two parameter sets
param_sets = [
    # Set 1: [r1, r2, K1, K2, a12, a21, b1, b2]
    [1.2, 1.0, 50.0, 40.0, 0.8, 0.6, 0.1, 0.12],
    # Set 2: Different competition coefficients
    [1.2, 1.0, 50.0, 40.0, 1.2, 1.0, 0.1, 0.12]
]

# Parameter that changes at t=toff (let's change r1 - growth rate of prey1)
function create_time_varying_params(base_params, change_factor, toff)
    function params_func(t)
        p = copy(base_params)
        if t >= toff
            p[1] = p[1] * change_factor  # Change r1
        end
        return p
    end
    return params_func
end

# Storage for solutions
solutions = []
labels = ["Parameter Set 1", "Parameter Set 2"]
colors = [colorant"#1f77b4", colorant"#ff7f0e"]

# Solve for both parameter sets
for (i, base_params) in enumerate(param_sets)
    println("Solving for parameter set $i...")
    
    # Create time-varying parameter function
    # At t=toff, r1 increases by 50%
    params_func = create_time_varying_params(base_params, 1.5, toff)
    
    # Split integration at toff to handle parameter change
    # Phase 1: t = 0 to toff
    prob1 = ODEProblem(ode_system!, u0, (0.0, toff), params_func(0.0))
    sol1 = solve(prob1, Tsit5(), saveat=0.1)
    
    # Phase 2: t = toff to end (using final state from phase 1)
    u_mid = sol1.u[end]
    prob2 = ODEProblem(ode_system!, u_mid, (toff, tspan[2]), params_func(toff + 0.1))
    sol2 = solve(prob2, Tsit5(), saveat=0.1)
    
    # Combine solutions
    t_combined = vcat(sol1.t, sol2.t[2:end])  # Avoid duplicate at toff
    u_combined = vcat(sol1.u, sol2.u[2:end])
    
    push!(solutions, (t=t_combined, u=u_combined))
end

# Create the plot
fig = Figure(resolution=(1200, 800))

# Define y-axis limits
y_limits = [(0, 60), (0, 15), (0, 50), (0, 12)]
titles = ["Prey Species 1", "Predator Species 1", "Prey Species 2", "Predator Species 2"]

for i in 1:4
    ax = Axis(fig[div(i-1, 2)+1, ((i-1) % 2)+1], 
             title=titles[i],
             xlabel="Time",
             ylabel="Population",
             limits=(nothing, y_limits[i]))
    
    # Plot solutions for both parameter sets
    for (j, sol) in enumerate(solutions)
        u_values = [u[i] for u in sol.u]
        lines!(ax, sol.t, u_values, 
               color=colors[j], 
               linewidth=2.5,
               label=labels[j])
    end
    
    # Add vertical line at parameter change time
    vlines!(ax, [toff], color=:red, linestyle=:dash, alpha=0.7, linewidth=2)
    
    # Add legend to first subplot
    if i == 1
        axislegend(ax, position=:rt)
    end
end

# Add overall title and parameter change annotation
Label(fig[0, :], "4-Species Population Dynamics with Parameter Step Change", 
      fontsize=18, font="bold")

# Add text annotation about parameter change
text_str = "Parameter change at t = $toff\n(r₁ increases by 50%)"
Label(fig[3, :], text_str, fontsize=12, color=:red)

# Display the figure
display(fig)

# Export to vectorized PDF
save("ode_simulation_results.pdf", fig, pt_per_unit=1)
println("Plot exported to 'ode_simulation_results.pdf'")

# Print some summary statistics
println("\nSummary:")
println("========")
println("Integration time span: $(tspan)")
println("Parameter change at t = $toff")
println("Number of parameter sets: $(length(param_sets))")
for (i, sol) in enumerate(solutions)
    final_state = sol.u[end]
    println("Parameter set $i - Final populations: [$(round.(final_state, digits=2))]")
end