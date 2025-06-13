using LinearAlgebra
using CairoMakie
using Random
using Statistics
using XLSX
using CSV
using DataFrames

include("discrete_bridge3.jl")

function batch_generate_brownian_bridges(input_file::String)
    """
    Read parameters from xlsx or csv file and generate PDF for each row.
    
    Expected format:
    - Column 1: filename (base name for PDF)
    - Column 2: x_end_x (x-coordinate of target end)
    - Column 3: x_end_y (y-coordinate of target end)
    - Column 4: delta (displacement magnitude)
    - Column 5: N (number of segments)
    - Column 6: seed (random seed, optional)
    """
    
    println("Reading parameters from: $input_file")
    
    # Read the file based on extension
    # if endswith(input_file, ".csv")
    df = CSV.read(input_file, DataFrame)
    headers = names(df)
    # elseif endswith(input_file, ".xlsx")
    #     xf = XLSX.readtable(input_file, "Sheet1", header=true)
    #     df = DataFrame(xf[1])  # Convert to DataFrame
    #     headers = xf[2]  # Column names
    # else
    #     error("Unsupported file format. Use .csv or .xlsx")
    # end
    
    println("Found $(nrow(df)) parameter sets to process")
    println("Headers: $(join(headers, ", "))")
    
    # Compute maximum axis limit for consistent scaling across all figures
    valid_rows = findall(i -> !ismissing(df.Figure[i]) && string(df.Figure[i]) != "", 1:nrow(df))
    if !isempty(valid_rows)
        max_axis_limit = maximum(Float64(df.delta[i]) * parse(Int, string(df.N[i])) for i in valid_rows)
        println("Using uniform axis limits: ±$(round(max_axis_limit, digits=2))")
    else
        max_axis_limit = nothing
        println("No valid rows found, using auto-scaling")
    end
    println()
    
    # Process each row
    for i in 1:nrow(df)
        # Skip blank rows (check if Figure column is missing or empty)
        if ismissing(df.Figure[i]) || string(df.Figure[i]) == ""
            println("Skipping blank row $i")
            continue
        end
        
        println("Processing row $i: $(df.Figure[i])_$(df.Panel[i])")
        
        try

            @show df[i,:]

            # Extract parameters
            figure_num = string(df.Figure[i])
            panel_name = string(df.Panel[i])
            x_end = [Float64(df.x_end_x[i]), Float64(df.x_end_y[i])]
            delta = Float64(df.delta[i])
            @show df.N[i]
            @show typeof(df.N[i])
            N = parse(Int, string(df.N[i]))

            
            # Handle optional seed
            if "seed" in names(df) && !ismissing(df.seed[i])
                Random.seed!(parse(Int, string(df.seed[i])))
                println("  Using seed: $(df.seed[i])")
            end
            
            println("  Parameters: x_end=$x_end, delta=$delta, N=$N")
            
            # Run simulation
            times, trajectory = brownian_bridge_2d_projected(x_end, delta, N)
            
            # Create plot with uniform axis limits
            fig = plot_trajectory(trajectory, x_end, delta; axis_limits=max_axis_limit)
            
            # Save with filename combining Figure and Panel columns
            pdf_filename = "fig_$(figure_num)_$(panel_name).pdf"
            save(pdf_filename, fig)
            println("  Saved: $pdf_filename")
            
            # Calculate and display final distance
            final_distance = norm(trajectory[:, end] - x_end)
            println("  Final distance from target: $(round(final_distance, digits=4))")
            println()
            
        catch e
            println("  ERROR processing row $i: $e")
            println()
            continue
        end
    end
    
    println("Batch processing complete!")
end



# For REPL debugging - call this function directly
batch_generate_brownian_bridges("sizes.csv")