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
    if endswith(input_file, ".csv")
        df = CSV.read(input_file, DataFrame)
        headers = names(df)
    elseif endswith(input_file, ".xlsx")
        xf = XLSX.readtable(input_file, "Sheet1", header=true)
        df = DataFrame(xf[1])  # Convert to DataFrame
        headers = xf[2]  # Column names
    else
        error("Unsupported file format. Use .csv or .xlsx")
    end
    
    println("Found $(nrow(df)) parameter sets to process")
    println("Headers: $(join(headers, ", "))")
    println()
    
    # Process each row
    for i in 1:nrow(df)
        println("Processing row $i: $(df.filename[i])")
        
        try
            # Extract parameters
            filename = string(df.filename[i])
            x_end = [Float64(df.x_end_x[i]), Float64(df.x_end_y[i])]
            delta = Float64(df.delta[i])
            N = Int(df.N[i])
            
            # Handle optional seed
            if "seed" in names(df) && !ismissing(df.seed[i])
                Random.seed!(Int(df.seed[i]))
                println("  Using seed: $(df.seed[i])")
            end
            
            println("  Parameters: x_end=$x_end, delta=$delta, N=$N")
            
            # Run simulation
            times, trajectory = brownian_bridge_2d_projected(x_end, delta, N)
            
            # Create plot
            fig = plot_trajectory(trajectory, x_end, delta)
            
            # Save with filename from spreadsheet (prefixed with "fig_")
            pdf_filename = "fig_$(filename).pdf"
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

# Example usage
if abspath(PROGRAM_FILE) == @__FILE__
    # Try CSV first, then xlsx
    input_files = ["sizes.csv", "sizes.xlsx"]
    
    for input_file in input_files
        if isfile(input_file)
            println("Found $input_file - using this file")
            batch_generate_brownian_bridges(input_file)
            break
        end
    end
    
    if !any(isfile.(input_files))
        println("Error: No input file found!")
        println("Please create either sizes.csv or sizes.xlsx with the expected format:")
        println("  Column 1: filename")
        println("  Column 2: x_end_x") 
        println("  Column 3: x_end_y")
        println("  Column 4: delta")
        println("  Column 5: N")
        println("  Column 6: seed (optional)")
    end
end