module NGSpiceTools

export NGSpiceResult

export simulate_vds_vgs_ids
export parse_ngspice_iv_curves, plot_iv_family_curves, parse_and_plot_iv_curves

export simulate_vgs_sweep, parse_and_plot_vgs_id, extract_vgs_id_data


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
    simulate_vds_vgs_ids(spice_file::String="./examples/sp_and_data/nmos_iv.sp")

Runs the NGSpice simulation defined in `spice_file`. It generates a `nmos_iv_data.txt` or `pmos_iv_data.txt` file containing the IV curve data for the mosfet.

# Arguments
- `spice_file::String`: Path to the NGSpice input file (.sp).

# Returns
- Nothing. Data is saved in `nmos_iv_data.txt` or `pmos_iv_data.txt`.

"""

function simulate_vds_vgs_ids(spice_file::String="./examples/sp_and_data/nmos_iv.sp")
    cmd = `ngspice -b $spice_file`
    run(cmd)
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

"""
    simulate_vgs_sweep(spice_file::String="./examples/sp_and_data/simulate_vgs_sweep.sp")

Runs the NGSpice simulation defined in `spice_file`. It generates a `vgs_id_data.txt` file containing
Vgs and Id values.

# Arguments
- `spice_file::String`: Path to the NGSpice input file (.sp).

# Returns
- Nothing. Data is saved in `vgs_id_data.txt`.
"""
function simulate_vgs_sweep(spice_file::String="./examples/sp_and_data/simulate_vgs_sweep.sp")
    cmd = `ngspice -b $spice_file`
    run(cmd)
end

"""
    extract_vgs_id_data(datafile::String="./examples/sp_and_data/vgs_id_data.txt")

Extracts Vgs and Id data from the NGSpice-generated data file and returns them in a dictionary.

# Arguments
- `datafile::String`: Path to the NGSpice-generated data file.

# Returns
- Dictionary with keys `"Vgs"` and `"Id"` containing respective data arrays.
"""
function extract_vgs_id_data(datafile::String="./examples/sp_and_data/vgs_id_data.txt")
    vgs = Float64[]
    id = Float64[]

    open(datafile, "r") do file
        for (idx, line) in enumerate(eachline(file))
            if idx == 1 && occursin("v(gate)", line)
                continue  # Skip header
            end
            vals = split(strip(line))
            if length(vals) >= 3
                push!(vgs, parse(Float64, vals[2]))
                push!(id, abs(parse(Float64, vals[3])))
            end
        end
    end

    return Dict("Vgs" => vgs, "Id" => id)
end

"""
    parse_and_plot_vgs_id(datafile::String="./examples/sp_and_data/vgs_id_data.txt")

Parses the NGSpice-generated Vgs vs. Id data file and plots Id as a function of Vgs.

# Arguments
- `datafile::String`: Path to the NGSpice-generated data file.

# Returns
- Plot object showing Id vs. Vgs.
"""
function parse_and_plot_vgs_id(datafile::String="./examples/sp_and_data/vgs_id_data.txt")
    data = extract_vgs_id_data(datafile)

    if isempty(data["Vgs"]) || isempty(data["Id"])
        error("No data found in file or incorrect data format")
    end

    plot(data["Vgs"], data["Id"], xlabel="V_{GS} (V)", ylabel="I_{D} (A)",
         title="Drain Current vs. Gate Voltage", linewidth=2, legend=false)
end


end # module