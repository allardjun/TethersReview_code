using OrdinaryDiffEq
using CairoMakie
using Colors
using ComponentArrays

function ode_system!(du, u, p, t)
    # Parameters
    kon, koff, kon_tethered = p

    transition_matrix = [
        -2*kon  koff  koff  0.0;
         kon   -koff-kon_tethered  0.0  +koff;
         kon    kon_tethered  -koff-kon_tethered  +koff;
         0.0    0.0  kon_tethered  -2*koff
    ]

    # System of 4 coupled ODEs
    du[:] = transition_matrix * u 
end

# Parameters and initial conditions
tspan = (0.0, 120.0)
toff = 60.0  # Time when parameter change occurs
u0 = [1.0, 0.0, 0.0, 0.0]  # Initial conditions

D = 40.0 # um^2/s
c_in_microMolar = 10.0 # uM
uMum3 = 602.2 # Conversion factor for uM to particles per um^3
c = c_in_microMolar * uMum3  # Convert to particles per um^3
l = 0.01 # um
R = 0.005 # um

alpha=1e-3


# Two parameter sets
param_sets = [
    # ComponentArray(
    #     kon = 0.1,
    #     koff = 0.1,
    #     kon_tethered = 0.1,
    # ),
    # ComponentArray(
    #     kon = 4*pi*R*D*c/uMum3,
    #     koff = 0.1,
    #     kon_tethered = 4*pi*R*D*1/(2*pi*l^2)^(3/2),
    # ),
    ComponentArray(
        kon = alpha*4*pi*R*D*c,
        koff = 100.0,
        kon_tethered = 4*pi*R*D*1/(2*pi*l^2)^(3/2),
    ),
]

# Parameter that changes at t=toff 
function create_time_varying_params(base_params, toff)
    function params_func(t)
        p = copy(base_params)
        if t >= toff
            p.kon = 0.0
        end
        return p
    end
    return params_func
end

# Storage for solutions
solutions = []
labels = ["Parameter Set 1", "Parameter Set 2", "Parameter Set 3"]
colors = [colorant"#1f77b4", colorant"#ff7f0e", colorant"#ff7f0e"]

# Solve for both parameter sets
for (i, base_params) in enumerate(param_sets)
    println("Solving for parameter set $i...")
    
    @show base_params

    # Create time-varying parameter function
    params_func = create_time_varying_params(base_params, toff)
    
    @show params_func(0.0)  # Initial parameters

    # Split integration at toff to handle parameter change
    # Phase 1: t = 0 to toff
    prob1 = ODEProblem(ode_system!, u0, (0.0, toff), params_func(0.0))
    sol1 = solve(prob1, Rosenbrock23(), saveat=0.1)
    
    # Phase 2: t = toff to end (using final state from phase 1)
    u_mid = sol1.u[end]
    prob2 = ODEProblem(ode_system!, u_mid, (toff, tspan[2]), params_func(toff + 0.1))
    sol2 = solve(prob2, Rosenbrock23(), saveat=0.1)
    
    # Combine solutions
    t_combined = vcat(sol1.t, sol2.t[2:end])  # Avoid duplicate at toff
    u_combined = vcat(sol1.u, sol2.u[2:end])
    
    push!(solutions, (t=t_combined, u=u_combined))
end

# Create the plot
fig = Figure(size=(600, 400))

# Define y-axis limits
y_limits = (0, 1)
title = "Tandem binding"

ax = Makie.Axis(fig[1,1], 
        title=title,
        xlabel="Time (seconds)",
        ylabel="Bound fraction",
        limits=(nothing, y_limits),
        )
    
# Plot solutions for both parameter sets
for (j, sol) in enumerate(solutions)
    u_values = [u[2].+u[3].+u[4] for u in sol.u]
    lines!(ax, sol.t, u_values, 
            color=colors[j],
            linewidth=2.5,
            label=labels[j])
end

# Add vertical line at parameter change time
vlines!(ax, [toff], color=:red, linestyle=:dash, alpha=0.7, linewidth=2)

# Add legend to first subplot
# if i == 1
#     axislegend(ax, position=:rt)
# end

# # Add overall title and parameter change annotation
# Label(fig[0, :], "4-Species Population Dynamics with Parameter Step Change", 
#       fontsize=18, font="bold")

# # Add text annotation about parameter change
# text_str = "Parameter change at t = $toff\n(r₁ increases by 50%)"
# Label(fig[3, :], text_str, fontsize=12, color=:red)

# Display the figure
display(fig)

# Export to vectorized PDF
save("fig_monovalent.pdf", fig, pt_per_unit=1)
