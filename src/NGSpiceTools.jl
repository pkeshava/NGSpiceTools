module NGSpiceTools

export parse_ngspice_log, parse_ngspice_raw
export extract_simulation_data, plot_iv_curves, parse_simulation_type
export NGSpiceResult, read_iv_data, simulate_iv_from_log
export parse_ngspice_iv_curves, plot_iv_family_curves, parse_and_plot_iv_curves

using Plots

"""
    NGSpiceResult

Structure containing parsed ngspice simulation results.

# Fields
- `simulation_type::String`: Type of simulation (dc, ac, tran, etc.)
- `variables::Vector{String}`: Names of variables in the simulation
- `data::Dict{String, Vector{Float64}}`: Data values for each variable
- `metadata::Dict{String, String}`: Additional simulation metadata
"""
struct NGSpiceResult
    simulation_type::String
    variables::Vector{String}
    data::Dict{String, Vector{Float64}}
    metadata::Dict{String, String}
end

"""
    parse_ngspice_log(filename::String) -> NGSpiceResult

Parse an ngspice log file and extract simulation data.

# Arguments
- `filename::String`: Path to the ngspice log file

# Returns
- `NGSpiceResult`: Structured data containing simulation results
"""
function parse_ngspice_log(filename::String)
    content = read(filename, String)
    lines = split(content, r"\r?\n")
    
    # Extract simulation type
    sim_type = parse_simulation_type(lines)
    
    # Extract data rows info
    data_rows = 0
    for line in lines
        m = match(r"No\. of Data Rows\s*:\s*(\d+)", line)
        if !isnothing(m)
            data_rows = parse(Int, m[1])
            break
        end
    end
    
    # Create empty result structure
    metadata = Dict{String, String}()
    metadata["data_rows"] = string(data_rows)
    
    return NGSpiceResult(
        sim_type,
        String[],  # We don't have variable names from log files
        Dict{String, Vector{Float64}}(),  # No data available directly in log
        metadata
    )
end

"""
    parse_ngspice_raw(filename::String) -> NGSpiceResult

Parse an ngspice raw file and extract simulation data.

# Arguments
- `filename::String`: Path to the ngspice raw file

# Returns
- `NGSpiceResult`: Structured data containing simulation results
"""
function parse_ngspice_raw(filename::String)
    content = read(filename, String)
    lines = split(content, r"\r?\n")
    
    # Initialize variables
    variables = String[]
    data = Dict{String, Vector{Float64}}()
    metadata = Dict{String, String}()
    simulation_type = "unknown"
    in_variables = false
    in_values = false
    num_points = 0
    data_format = "unknown"
    
    i = 1
    while i <= length(lines)
        line = lines[i]
        
        # Parse metadata section
        if startswith(line, "Title: ")
            metadata["title"] = strip(line[8:end])
        elseif startswith(line, "Date: ")
            metadata["date"] = strip(line[7:end])
        elseif startswith(line, "Plotname: ")
            simulation_type = lowercase(strip(line[11:end]))
        elseif startswith(line, "Flags: ")
            metadata["flags"] = strip(line[8:end])
        elseif startswith(line, "No. Variables: ")
            num_vars = parse(Int, strip(line[15:end]))
            metadata["num_variables"] = string(num_vars)
        elseif startswith(line, "No. Points: ")
            num_points = parse(Int, strip(line[12:end]))
            metadata["num_points"] = string(num_points)
        elseif startswith(line, "Variables:")
            in_variables = true
            i += 1
            # Parse variable names
            for j in 1:num_vars
                if i > length(lines) 
                    break
                end
                var_parts = split(lines[i])
                if length(var_parts) >= 3
                    push!(variables, var_parts[2])
                    data[var_parts[2]] = Float64[]
                end
                i += 1
            end
            in_variables = false
            i -= 1  # Adjust for the outer loop increment
        elseif startswith(line, "Values:")
            in_values = true
            i += 1
            # Parse data values
            point_idx = 0
            while point_idx < num_points && i <= length(lines)
                if match(r"^\s*\d+\s*:", lines[i]) !== nothing
                    # This is an index line
                    point_idx += 1
                    i += 1
                    
                    # Read values for each variable
                    for j in 1:length(variables)
                        if i > length(lines)
                            break
                        end
                        if !isempty(lines[i])
                            val = parse(Float64, strip(lines[i]))
                            push!(data[variables[j]], val)
                        end
                        i += 1
                    end
                    i -= 1  # Adjust for outer loop increment
                else
                    i += 1
                end
            end
            in_values = false
        elseif startswith(line, "Binary:")
            metadata["format"] = "binary"
            # Binary parsing would go here
            break
        end
        
        i += 1
    end
    
    return NGSpiceResult(
        simulation_type,
        variables,
        data,
        metadata
    )
end

"""
    parse_simulation_type(lines::Vector{<:AbstractString}) -> String

Extract the simulation type from ngspice output lines.

# Arguments
- `lines::Vector{<:AbstractString}`: Lines from the ngspice output file

# Returns
- `String`: The identified simulation type or "unknown"
"""
function parse_simulation_type(lines::Vector{<:AbstractString})
    for line in lines
        if occursin(r"\.dc\s+", line)
            return "dc"
        elseif occursin(r"\.ac\s+", line)
            return "ac"
        elseif occursin(r"\.tran\s+", line)
            return "tran"
        elseif occursin(r"IV Curves", line)
            return "iv"
        end
    end
    return "unknown"
end

"""
    extract_simulation_data(spfile::String, raw_file::String, var_names::Vector{String})

Extract specific variable data from a raw file based on the simulation setup in spfile.

# Arguments
- `spfile::String`: Path to the ngspice circuit file
- `raw_file::String`: Path to the raw output file
- `var_names::Vector{String}`: Variable names to extract

# Returns
- `Dict{String, Vector{Float64}}`: Dictionary of variable data
"""
function extract_simulation_data(spfile::String, raw_file::String, var_names::Vector{String})
    # Read circuit file to understand simulation context
    sp_content = read(spfile, String)
    sp_lines = split(sp_content, r"\r?\n")
    
    # Parse raw file
    result = parse_ngspice_raw(raw_file)
    
    # Extract and return specific variables
    extracted_data = Dict{String, Vector{Float64}}()
    for var in var_names
        if haskey(result.data, var)
            extracted_data[var] = result.data[var]
        end
    end
    
    return extracted_data
end

"""
    plot_iv_curves(ngspice_result::NGSpiceResult; 
                  x_var::String="v-sweep", 
                  y_var::String="i-sweep",
                  title::String="IV Curves",
                  xlabel::String="Voltage (V)",
                  ylabel::String="Current (A)")

Generate IV curve plots from ngspice simulation results.

# Arguments
- `ngspice_result::NGSpiceResult`: Parsed ngspice results
- `x_var::String`: Name of the voltage variable
- `y_var::String`: Name of the current variable
- `title::String`: Plot title
- `xlabel::String`: X-axis label
- `ylabel::String`: Y-axis label

# Returns
- Plots.Plot: The generated plot
"""
function plot_iv_curves(ngspice_result::NGSpiceResult; 
                        x_var::String="v-sweep", 
                        y_var::String="i-sweep",
                        title::String="IV Curves",
                        xlabel::String="Voltage (V)",
                        ylabel::String="Current (A)")
    
    if !haskey(ngspice_result.data, x_var) || !haskey(ngspice_result.data, y_var)
        error("Variables $x_var or $y_var not found in the simulation data")
    end
    
    x_data = ngspice_result.data[x_var]
    y_data = ngspice_result.data[y_var]
    
    plot(x_data, y_data, 
         title=title, 
         xlabel=xlabel, 
         ylabel=ylabel, 
         label="IV Curve",
         linewidth=2)
end

"""
    read_iv_data(datafile::String; is_nmos::Bool=true)

Read IV curve data from a text file saved by NGSpice.
Applies correct sign conventions for NMOS/PMOS devices.

# Arguments
- `datafile::String`: Path to NGSpice data file (text format)
- `is_nmos::Bool`: Whether the device is an NMOS (true) or PMOS (false)

# Returns
- `Dict{String, Vector{Float64}}`: Dictionary with voltage and current data
"""
function read_iv_data(datafile::String; is_nmos::Bool=true)
    if !isfile(datafile)
        error("Data file not found: $datafile")
    end
    
    data = Dict{String, Vector{Float64}}()
    
    # Read text format file
    lines = readlines(datafile)
    
    # First line might be headers
    headers = split(lines[1])
    if all(x -> !occursin(r"^[-+]?\d", x), headers)
        # This is a header line
        data_start = 2
        var_names = headers
    else
        # No headers, just data
        data_start = 1
        var_names = ["v", "i"]  # Default names
    end
    
    # Initialize data arrays
    for name in var_names
        data[name] = Float64[]
    end
    
    # Read data values
    for i in data_start:length(lines)
        values = split(lines[i])
        if length(values) == length(var_names)
            for (j, name) in enumerate(var_names)
                if !isempty(values[j])
                    push!(data[name], parse(Float64, values[j]))
                end
            end
        end
    end
    
    # Apply correct sign conventions based on device type
    # We assume the first variable is voltage and the second is current
    if length(var_names) >= 2
        v_name = var_names[1]
        i_name = var_names[2]
        
        if is_nmos
            # For NMOS: VDS should be positive, IDS should be positive
            data[v_name] = abs.(data[v_name])
            data[i_name] = abs.(data[i_name])
        else
            # For PMOS: VDS should be negative, IDS should be positive
            data[v_name] = -abs.(data[v_name])
            data[i_name] = abs.(data[i_name])
        end
    end
    
    return data
end

"""
    simulate_iv_from_log(logfile::String; datafile::String="", is_nmos::Bool=true)

Parse log file and look for accompanying data file to plot IV curves.
If no data file is provided, it will look for default data filenames.
Applies correct sign conventions for NMOS/PMOS devices.

# Arguments
- `logfile::String`: Path to NGSpice log file
- `datafile::String`: (Optional) Path to data file
- `is_nmos::Bool`: Whether the device is an NMOS (true) or PMOS (false)

# Returns
- `Plots.Plot`: A plot of the IV curves from the data
"""
function simulate_iv_from_log(logfile::String; datafile::String="", is_nmos::Bool=true)
    log_result = parse_ngspice_log(logfile)
    
    # Try to auto-detect device type from logfile if not specified
    if !is_nmos && !occursin(r"pmos|pfet|pchan", lowercase(logfile))
        # Check if the user explicitly set is_nmos=false
        # If is_nmos was not explicitly set, try to detect from log content
        if occursin(r"nmos|nfet|nchan", lowercase(logfile))
            is_nmos = true
            @info "Detected NMOS device from filename"
        else
            # Check log content to determine device type
            content = read(logfile, String)
            if occursin(r"nmos ", lowercase(content))
                is_nmos = true
                @info "Detected NMOS device from log content"
            elseif occursin(r"pmos ", lowercase(content))
                is_nmos = false
                @info "Detected PMOS device from log content"
            end
            # Otherwise keep the is_nmos value as provided
        end
    elseif occursin(r"pmos|pfet|pchan", lowercase(logfile))
        is_nmos = false
        @info "Detected PMOS device from filename"
    end
    
    # Device type for display purposes
    device_type = is_nmos ? "NMOS" : "PMOS"
    
    # Generate possible data filenames to look for
    possible_files = String[]
    
    # If datafile is specified directly, use that
    if !isempty(datafile)
        push!(possible_files, datafile)
    else
        # Try to infer data filenames based on the log filename
        base_dir = dirname(logfile)
        base_name = replace(basename(logfile), r"_results\.log$|\.log$" => "")
        
        # Common data file patterns - only look for text files
        push!(possible_files, joinpath(base_dir, "$(base_name)_iv_data.txt"))
        push!(possible_files, joinpath(base_dir, "$(base_name)_data.txt"))
        push!(possible_files, joinpath(base_dir, "$(base_name).txt"))
        # Look for the raw file we received
        push!(possible_files, joinpath(base_dir, "$(base_name)_iv_data copy.txt"))
        # Also look for .dat files which are often used for data export
        push!(possible_files, joinpath(base_dir, "$(base_name)_iv_data.dat"))
        push!(possible_files, joinpath(base_dir, "$(base_name)_data.dat"))
        push!(possible_files, joinpath(base_dir, "$(base_name).dat"))
    end
    
    # Try to find and read a data file
    data = Dict{String, Vector{Float64}}()
    found_file = ""
    
    for file in possible_files
        if isfile(file)
            try
                data = read_iv_data(file, is_nmos=is_nmos)
                found_file = file
                break
            catch e
                @warn "Failed to read data from $file: $e"
            end
        end
    end
    
    if isempty(found_file)
        error("""
        No usable data file found. Please run NGSpice with the proper save commands.
        
        Add these lines to your .control section:
            set wr_vecnames
            set wr_singlescale
            wrdata $(basename(logfile, ".log"))_iv_data.txt v(d) i(VDS)
            
        Then run ngspice again.
        """)
    end
    
    # Get variable names for plotting
    var_names = collect(keys(data))
    if length(var_names) < 2
        error("Insufficient data columns in $found_file")
    end
    
    # Determine which variables to plot (voltage and current)
    # If we can't find specific names, just use the first two columns
    v_var = ""
    i_var = ""
    
    # Try to find voltage and current variables by name
    for name in var_names
        if occursin(r"v\(|voltage|vds", lowercase(name))
            v_var = name
        elseif occursin(r"i\(|current|id", lowercase(name))
            i_var = name
        end
    end
    
    # If we couldn't find them by name, use the first two columns
    if isempty(v_var) && length(var_names) >= 1
        v_var = var_names[1]
    end
    if isempty(i_var) && length(var_names) >= 2
        i_var = var_names[2]
    end
    
    # Create plot
    if log_result.simulation_type == "iv" || log_result.simulation_type == "dc"
        plt = plot(title="$device_type IV Characteristic Curves",
                  xlabel="Drain Voltage (V)",
                  ylabel="Drain Current (A)")
        
        # Plot the data
        if !isempty(data[v_var]) && !isempty(data[i_var])
            plot!(plt, data[v_var], data[i_var], label="IV Curve", linewidth=2)
        else
            @warn "No data points found in $found_file"
        end
        
        return plt
    else
        # Generic plot for other simulation types
        return plot(data[v_var], data[i_var], 
                  title="$device_type NGSpice Simulation Data",
                  xlabel="Voltage (V)", 
                  ylabel="Current (A)",
                  label="Simulation", 
                  linewidth=2)
    end
end

"""
    parse_ngspice_iv_curves(datafile::String; is_nmos::Bool=true)

Parse IV curve data from an NGSpice data file that contains multiple sweeps for different VGS values.
This function handles the specific format where each VGS value has its own block of VDS sweeps.

# Arguments
- `datafile::String`: Path to the NGSpice data file
- `is_nmos::Bool`: Whether the device is an NMOS (true) or PMOS (false)

# Returns
- `Dict{String, Any}`: Contains parsed data with the following keys:
  - "vgs_values": Vector of VGS values
  - "vds": Matrix where each column corresponds to a VGS value
  - "ids": Matrix where each column corresponds to a VGS value
"""
function parse_ngspice_iv_curves(datafile::String; is_nmos::Bool=true)
    if !isfile(datafile)
        error("Data file not found: $datafile")
    end
    
    # Read all lines from the file
    lines = readlines(datafile)
    
    # Get header line and parse column names
    header = strip(lines[1])
    column_names = split(header)
    
    # Identify block start lines (where VDS = 0)
    block_starts = Int[]
    for i in 2:length(lines)
        line = strip(lines[i])
        if startswith(line, "0.00000000e+00")
            push!(block_starts, i)
        end
    end
    
    if isempty(block_starts)
        push!(block_starts, 2)  # If no explicit blocks, start at line 2 (after header)
    end
    
    # If first block doesn't start at line 2, add line 2 as the first block
    if block_starts[1] > 2
        pushfirst!(block_starts, 2)
    end
    
    # Extract data for each VGS value block
    all_vds = []
    all_ids = []
    
    for i in 1:length(block_starts)
        start_line = block_starts[i]
        end_line = i < length(block_starts) ? block_starts[i+1] - 1 : length(lines)
        
        vds = Float64[]
        ids = Float64[]
        
        for j in start_line:end_line
            line = strip(lines[j])
            if isempty(line) || line == header
                continue
            end
            
            values = split(line)
            if length(values) >= 3
                # Column 2 is v(d) - drain voltage
                # Column 3 is i(VDS) - drain current
                push!(vds, parse(Float64, values[2]))
                
                # Apply sign convention based on device type
                current = parse(Float64, values[3])
                if is_nmos
                    # For NMOS, drain current is typically negative, make it positive
                    push!(ids, abs(current))
                else
                    # For PMOS, keep negative convention for VDS, positive for IDS
                    push!(ids, abs(current))
                end
            end
        end
        
        if !isempty(vds) && !isempty(ids)
            push!(all_vds, vds)
            push!(all_ids, ids)
        end
    end
    
    # Check if we have any valid data
    if isempty(all_vds) || isempty(all_ids)
        error("No valid IV curve data found in the file")
    end
    
    # Determine VGS values - they are not directly in the data file
    # So we need to infer them from the simulation setup (usually 0.2V steps from VT)
    # For a typical 65nm NMOS, threshold is around 0.4V, so we might start at 0.2V
    vgs_start = is_nmos ? 0.2 : -0.2
    vgs_step = is_nmos ? 0.2 : -0.2
    vgs_values = [vgs_start + (i-1) * vgs_step for i in 1:length(all_vds)]
    
    # Convert to matrices for easier plotting
    # Find the maximum length of any VDS/IDS array
    max_length = maximum(length.(all_vds))
    
    # Create matrices and fill with NaN for missing values
    vds_matrix = fill(NaN, max_length, length(all_vds))
    ids_matrix = fill(NaN, max_length, length(all_vds))
    
    for i in 1:length(all_vds)
        vds_data = all_vds[i]
        ids_data = all_ids[i]
        n = length(vds_data)
        vds_matrix[1:n, i] = vds_data
        ids_matrix[1:n, i] = ids_data
    end
    
    return Dict{String, Any}(
        "vgs_values" => vgs_values,
        "vds" => vds_matrix,
        "ids" => ids_matrix
    )
end

"""
    plot_iv_family_curves(data::Dict{String, Any}; 
                         title::String="IV Characteristics",
                         is_nmos::Bool=true)

Plot a family of IV curves for a MOSFET at different VGS values.

# Arguments
- `data::Dict{String, Any}`: Data dictionary from parse_ngspice_iv_curves
- `title::String`: Plot title
- `is_nmos::Bool`: Whether the device is an NMOS (true) or PMOS (false)

# Returns
- `Plots.Plot`: The generated plot
"""
function plot_iv_family_curves(data::Dict{String, Any}; 
                              title::String="IV Characteristics",
                              is_nmos::Bool=true)    
    vgs_values = data["vgs_values"]
    vds_matrix = data["vds"]
    ids_matrix = data["ids"]
    
    device_type = is_nmos ? "NMOS" : "PMOS"
    
    plt = plot(
        title="$device_type $title",
        xlabel="Drain-Source Voltage (V)",
        ylabel="Drain Current (A)",
        legend=:bottomright,
        linewidth=2
    )
    
    # Plot each VGS curve
    for i in 1:length(vgs_values)
        # Get valid data points (not NaN)
        valid_indices = .!isnan.(vds_matrix[:, i])
        
        if any(valid_indices)
            vds = vds_matrix[valid_indices, i]
            ids = ids_matrix[valid_indices, i]
            
            plot!(plt, vds, ids, label="VGS = $(vgs_values[i])V")
        end
    end
    
    return plt
end

"""
    parse_and_plot_iv_curves(datafile::String; 
                            title::String="IV Characteristics",
                            is_nmos::Bool=true)

Parse and plot a family of IV curves from an NGSpice data file in one step.

# Arguments
- `datafile::String`: Path to the NGSpice data file
- `title::String`: Plot title
- `is_nmos::Bool`: Whether the device is an NMOS (true) or PMOS (false)

# Returns
- `Plots.Plot`: The generated plot
"""
function parse_and_plot_iv_curves(datafile::String; 
                                 title::String="IV Characteristics",
                                 is_nmos::Bool=true)
    data = parse_ngspice_iv_curves(datafile, is_nmos=is_nmos)
    return plot_iv_family_curves(data, title=title, is_nmos=is_nmos)
end

end # module