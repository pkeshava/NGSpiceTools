module NGSpiceTools

export NGSpiceResult, NGSpiceSimulation
export SimulationType, DeviceType
export DC_SWEEP, AC_ANALYSIS, TRANSIENT, OPERATING_POINT, TRANSFER_CURVE, IV_CURVE
export NMOS, PMOS, RESISTOR, CAPACITOR, INDUCTOR, BJT, DIODE, OTHER
export run_simulation, read_simulation_output
export create_simulation, simulate_iv_curves, simulate_transfer_curve
export plot_iv_family_curves, plot_transfer_curve

using Plots

"""
    SimulationType

Enumeration of supported simulation types.
"""
@enum SimulationType begin
    DC_SWEEP
    AC_ANALYSIS
    TRANSIENT
    OPERATING_POINT
    TRANSFER_CURVE
    IV_CURVE
end

"""
    DeviceType

Enumeration of supported device types.
"""
@enum DeviceType begin
    NMOS
    PMOS
    RESISTOR
    CAPACITOR
    INDUCTOR
    BJT
    DIODE
    OTHER
end

"""
    NGSpiceResult

Structure containing parsed NGSpice simulation results.

# Fields
- `simulation_type::SimulationType`: Type of simulation that was performed
- `device_type::DeviceType`: Type of device being simulated
- `variables::Vector{String}`: Names of variables in the simulation
- `data::Dict{String, Vector{Float64}}`: Data values for each variable
- `metadata::Dict{String, Any}`: Additional simulation metadata
"""
struct NGSpiceResult
    simulation_type::SimulationType
    device_type::DeviceType
    variables::Vector{String}
    data::Dict{String, Vector{Float64}}
    metadata::Dict{String, Any}
end

"""
    NGSpiceSimulation

Structure representing an NGSpice simulation configuration.

# Fields
- `spice_file::String`: Path to the NGSpice input file (.sp)
- `output_file::String`: Path where simulation output should be saved
- `simulation_type::SimulationType`: Type of simulation to perform
- `device_type::DeviceType`: Type of device being simulated
- `parameters::Dict{String, Any}`: Simulation parameters
"""
struct NGSpiceSimulation
    spice_file::String
    output_file::String
    simulation_type::SimulationType
    device_type::DeviceType
    parameters::Dict{String, Any}
end

"""
    create_simulation(;
        spice_file::String,
        output_file::String="",
        simulation_type::SimulationType=DC_SWEEP,
        device_type::DeviceType=OTHER,
        parameters::Dict{String, Any}=Dict{String, Any}()
    ) -> NGSpiceSimulation

Create an NGSpiceSimulation object with the specified parameters.

# Arguments
- `spice_file::String`: Path to the NGSpice input file (.sp)
- `output_file::String`: Path where simulation output should be saved (inferred from spice_file if empty)
- `simulation_type::SimulationType`: Type of simulation to perform
- `device_type::DeviceType`: Type of device being simulated
- `parameters::Dict{String, Any}`: Simulation parameters

# Returns
- `NGSpiceSimulation`: Configured simulation object
"""
function create_simulation(;
    spice_file::String,
    output_file::String="",
    simulation_type::SimulationType=DC_SWEEP,
    device_type::DeviceType=OTHER,
    parameters::Dict{String, Any}=Dict{String, Any}()
) :: NGSpiceSimulation
    
    # If output file not specified, infer from spice file name
    if isempty(output_file)
        # Get directory and base name
        dir_path = dirname(spice_file)
        base_name = splitext(basename(spice_file))[1]
        
        # Create output file path in the same directory
        if simulation_type == IV_CURVE
            output_file = joinpath(dir_path, "$(base_name)_data.txt")
        elseif simulation_type == TRANSFER_CURVE
            output_file = joinpath(dir_path, "vgs_id_data.txt")  # Default for transfer curves
        else
            output_file = joinpath(dir_path, "$(base_name)_output.txt")
        end
    end
    
    return NGSpiceSimulation(
        spice_file,
        output_file,
        simulation_type,
        device_type,
        parameters
    )
end

"""
    run_simulation(simulation::NGSpiceSimulation; modify_file::Bool=false, working_dir::String="")

Run an NGSpice simulation with the given configuration.

# Arguments
- `simulation::NGSpiceSimulation`: The simulation configuration
- `modify_file::Bool=false`: Whether to modify the spice file before running
- `working_dir::String=""`: Working directory for the simulation

# Returns
- `Bool`: True if simulation completed successfully
"""
function run_simulation(simulation::NGSpiceSimulation; modify_file::Bool=false, working_dir::String="")
    # Use the specified working directory or extract from spice file path
    if isempty(working_dir)
        working_dir = dirname(simulation.spice_file)
        if isempty(working_dir) || working_dir == "."
            working_dir = pwd()
        end
    end
    
    # Use just the filename when calling ngspice from the working directory
    spice_filename = basename(simulation.spice_file)
    
    try
        # Change to the working directory before running ngspice
        cd(working_dir) do
            # Run ngspice with just the filename since we're already in the correct directory
            cmd = `ngspice -b $spice_filename`
            run(cmd)
        end
        
        return true
    catch e
        @error "Simulation failed" exception=(e, catch_backtrace()) working_dir=working_dir spice_file=simulation.spice_file
        return false
    end
end

"""
    read_simulation_output(simulation::NGSpiceSimulation) -> NGSpiceResult

Read and parse the output from an NGSpice simulation.

# Arguments
- `simulation::NGSpiceSimulation`: The simulation configuration

# Returns
- `NGSpiceResult`: Structured simulation results
"""
function read_simulation_output(simulation::NGSpiceSimulation)
    # Get the expected output file
    output_file = simulation.output_file
    
    # If the output file isn't specified or doesn't exist, try to find it based on common patterns
    if isempty(output_file) || !isfile(output_file)
        # Get the directory and base name
        dir_path = dirname(simulation.spice_file)
        base_name = splitext(basename(simulation.spice_file))[1]
        
        # List of possible output file patterns to try
        file_patterns = [
            # Pattern from spice file
            joinpath(dir_path, "$(base_name)_data.txt"),
            # Common NGSpice output file names
            joinpath(dir_path, "vgs_id_data.txt"),
            joinpath(dir_path, "iv_data.txt"),
            joinpath(dir_path, "nmos_iv_data.txt"),
            joinpath(dir_path, "pmos_iv_data.txt")
        ]
        
        # Try each pattern
        for pattern in file_patterns
            if isfile(pattern)
                output_file = pattern
                break
            end
        end
    end
    
    # Check if output file exists after trying all patterns
    if !isfile(output_file)
        error("Output file not found: $(output_file). Tried common patterns in $(dirname(simulation.spice_file))")
    end
    
    # Update the simulation object with the found output file
    simulation = NGSpiceSimulation(
        simulation.spice_file,
        output_file,
        simulation.simulation_type,
        simulation.device_type,
        simulation.parameters
    )
    
    # Parse the file based on simulation type
    if simulation.simulation_type == IV_CURVE
        return parse_iv_curve_data(simulation)
    elseif simulation.simulation_type == TRANSFER_CURVE
        return parse_transfer_curve_data(simulation)
    else
        # Generic parser for other simulation types
        return parse_generic_data(simulation)
    end
end

"""
    parse_generic_data(simulation::NGSpiceSimulation) -> NGSpiceResult

Generic parser for NGSpice output data.

# Arguments
- `simulation::NGSpiceSimulation`: The simulation configuration

# Returns
- `NGSpiceResult`: Structured simulation results
"""
function parse_generic_data(simulation::NGSpiceSimulation)
    lines = readlines(simulation.output_file)
    
    # Parse header to get variable names
    header = strip(lines[1])
    variables = split(header)
    
    # Initialize data dictionary
    data = Dict{String, Vector{Float64}}()
    for var in variables
        data[var] = Float64[]
    end
    
    # Parse data lines
    for i in 2:length(lines)
        line = strip(lines[i])
        if isempty(line)
            continue
        end
        
        values = split(line)
        if length(values) >= length(variables)
            for (j, var) in enumerate(variables)
                push!(data[var], parse(Float64, values[j]))
            end
        end
    end
    
    # Create and return NGSpiceResult
    return NGSpiceResult(
        simulation.simulation_type,
        simulation.device_type,
        variables,
        data,
        Dict{String, Any}("file" => simulation.output_file)
    )
end

"""
    parse_iv_curve_data(simulation::NGSpiceSimulation) -> NGSpiceResult

Parse NGSpice output data for IV curve simulations with multiple VGS values.

# Arguments
- `simulation::NGSpiceSimulation`: The simulation configuration

# Returns
- `NGSpiceResult`: Structured simulation results with IV curve data
"""
function parse_iv_curve_data(simulation::NGSpiceSimulation)
    is_nmos = simulation.device_type == NMOS
    
    # Read all lines from the file
    lines = readlines(simulation.output_file)
    
    # Get header line and parse column names
    header = strip(lines[1])
    variables = split(header)
    
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
                # Column 2 is v(d) - drain voltage (or v(d,s) for relative voltage)
                # Column 3 is i(VDS) - drain current
                push!(vds, parse(Float64, values[2]))
                
                # Apply sign convention based on device type
                current = parse(Float64, values[3])
                push!(ids, abs(current)) # Always store as positive
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
    
    # Extract or infer VGS values from simulation parameters or filename
    vgs_values = get(simulation.parameters, "vgs_values", [])
    
    if isempty(vgs_values)
        # Infer VGS values from the number of blocks
        vgs_start = is_nmos ? 0.2 : -0.2
        vgs_step = is_nmos ? 0.2 : -0.2
        vgs_values = [vgs_start + (i-1) * vgs_step for i in 1:length(all_vds)]
    end
    
    # Create matrices for the data
    max_length = maximum(length.(all_vds))
    vds_matrix = fill(NaN, max_length, length(all_vds))
    ids_matrix = fill(NaN, max_length, length(all_vds))
    
    for i in 1:length(all_vds)
        vds_data = all_vds[i]
        ids_data = all_ids[i]
        n = length(vds_data)
        vds_matrix[1:n, i] = vds_data
        ids_matrix[1:n, i] = ids_data
    end
    
    # Create structured data dictionary
    data = Dict{String, Vector{Float64}}()
    for var in variables
        data[var] = Float64[]  # Initialize as empty
    end
    
    # Create metadata dictionary
    metadata = Dict{String, Any}(
        "vgs_values" => vgs_values,
        "vds_matrix" => vds_matrix,
        "ids_matrix" => ids_matrix,
        "file" => simulation.output_file
    )
    
    return NGSpiceResult(
        simulation.simulation_type,
        simulation.device_type,
        variables,
        data,
        metadata
    )
end

"""
    parse_transfer_curve_data(simulation::NGSpiceSimulation) -> NGSpiceResult

Parse NGSpice output data for Vgs-Id (transfer characteristic) simulations.

# Arguments
- `simulation::NGSpiceSimulation`: The simulation configuration

# Returns
- `NGSpiceResult`: Structured simulation results with transfer curve data
"""
function parse_transfer_curve_data(simulation::NGSpiceSimulation)
    if !isfile(simulation.output_file)
        error("Output file not found: $(simulation.output_file)")
    end
    
    # Read data file
    lines = readlines(simulation.output_file)
    
    # Get header line and parse column names
    header = strip(lines[1])
    variables = split(header)
    
    # Initialize data vectors
    vgs = Float64[]
    id = Float64[]
    
    # Parse data lines
    for i in 2:length(lines)
        line = strip(lines[i])
        if isempty(line)
            continue
        end
        
        values = split(line)
        if length(values) >= 3
            # For transfer curves, we expect columns: v-sweep, v(gate), i(Vds)
            push!(vgs, parse(Float64, values[2]))  # v(gate)
            push!(id, abs(parse(Float64, values[3])))  # i(Vds), use abs to make current positive
        end
    end
    
    # Create structured data dictionary
    data = Dict{String, Vector{Float64}}(
        "Vgs" => vgs,
        "Id" => id
    )
    
    # Create metadata dictionary
    metadata = Dict{String, Any}(
        "file" => simulation.output_file
    )
    
    return NGSpiceResult(
        simulation.simulation_type,
        simulation.device_type,
        ["Vgs", "Id"],
        data,
        metadata
    )
end

"""
    plot_iv_family_curves(result::NGSpiceResult; title::String="IV Characteristics")

Plot a family of IV curves from simulation results.

# Arguments
- `result::NGSpiceResult`: Simulation results containing IV curve data
- `title::String`: Plot title

# Returns
- `Plots.Plot`: The generated plot
"""
function plot_iv_family_curves(result::NGSpiceResult; title::String="IV Characteristics")
    # Check that this is an IV curve result
    if result.simulation_type != IV_CURVE
        error("Result is not from an IV curve simulation")
    end
    
    # Extract data from the metadata
    vgs_values = result.metadata["vgs_values"]
    vds_matrix = result.metadata["vds_matrix"]
    ids_matrix = result.metadata["ids_matrix"]
    
    # Determine device type for proper labeling
    device_type = result.device_type == NMOS ? "NMOS" : "PMOS"
    
    # Create plot
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
    plot_transfer_curve(result::NGSpiceResult; 
                       title::String="Drain Current vs. Gate Voltage",
                       scale::Symbol=:linear)

Plot a transfer characteristic curve (Id vs. Vgs) from simulation results.

# Arguments
- `result::NGSpiceResult`: Simulation results containing transfer curve data
- `title::String`: Plot title
- `scale::Symbol`: Y-axis scale, either :linear or :log

# Returns
- `Plots.Plot`: The generated plot
"""
function plot_transfer_curve(result::NGSpiceResult; 
                            title::String="Drain Current vs. Gate Voltage",
                            scale::Symbol=:linear)
    # Check that this is a transfer curve result
    if result.simulation_type != TRANSFER_CURVE
        error("Result is not from a transfer curve simulation")
    end
    
    # Extract data
    vgs = result.data["Vgs"]
    id = result.data["Id"]
    
    # Determine device type for proper labeling
    device_type = result.device_type == NMOS ? "NMOS" : "PMOS"
    
    # Create plot with appropriate scale
    if scale == :log
        plt = plot(vgs, id, 
             xlabel="V_{GS} (V)", 
             ylabel="I_{D} (A)",
             title="$device_type $title", 
             linewidth=2, 
             legend=false,
             yaxis=:log10)
    else
        plt = plot(vgs, id, 
             xlabel="V_{GS} (V)", 
             ylabel="I_{D} (A)",
             title="$device_type $title", 
             linewidth=2, 
             legend=false)
    end
    
    return plt
end

"""
    simulate_iv_curves(spice_file::String; 
                      device_type::DeviceType=NMOS,
                      output_file::String="",
                      vgs_values::Vector{Float64}=Float64[],
                      working_dir::String="")

Run an IV curve simulation and return the results.

# Arguments
- `spice_file::String`: Path to the NGSpice input file (.sp)
- `device_type::DeviceType`: Type of device being simulated (NMOS or PMOS)
- `output_file::String`: Path where simulation output should be saved
- `vgs_values::Vector{Float64}`: Optional explicit VGS values for the simulation
- `working_dir::String`: Working directory for the simulation (default: directory of spice_file)

# Returns
- `NGSpiceResult`: Structured simulation results
"""
function simulate_iv_curves(spice_file::String; 
                           device_type::DeviceType=NMOS,
                           output_file::String="",
                           vgs_values::Vector{Float64}=Float64[],
                           working_dir::String="")
    
    # Determine working directory if not provided
    if isempty(working_dir)
        working_dir = dirname(spice_file)
        if isempty(working_dir) || working_dir == "."
            working_dir = pwd()
        end
    end
    
    # Resolve paths - make sure we don't duplicate directory information
    # For the spice file, use just the filename if it's in the working directory
    base_spice_file = basename(spice_file)
    
    # Create full paths based on whether the input paths are absolute
    resolved_spice_file = isabspath(spice_file) ? spice_file : joinpath(working_dir, base_spice_file)
    resolved_output_file = isempty(output_file) ? "" : 
                          (isabspath(output_file) ? output_file : joinpath(working_dir, basename(output_file)))
    
    parameters = Dict{String, Any}()
    if !isempty(vgs_values)
        parameters["vgs_values"] = vgs_values
    end
    
    # Store the working directory in parameters
    parameters["working_dir"] = working_dir
    
    simulation = create_simulation(
        spice_file=resolved_spice_file,
        output_file=resolved_output_file,
        simulation_type=IV_CURVE,
        device_type=device_type,
        parameters=parameters
    )
    
    success = run_simulation(simulation, working_dir=working_dir)
    if !success
        error("Failed to run IV curve simulation")
    end
    
    return read_simulation_output(simulation)
end

"""
    simulate_transfer_curve(spice_file::String; 
                           device_type::DeviceType=NMOS,
                           output_file::String="",
                           working_dir::String="")

Run a transfer curve (Vgs-Id) simulation and return the results.

# Arguments
- `spice_file::String`: Path to the NGSpice input file (.sp)
- `device_type::DeviceType`: Type of device being simulated (NMOS or PMOS)
- `output_file::String`: Path where simulation output should be saved
- `working_dir::String`: Working directory for the simulation (default: directory of spice_file)

# Returns
- `NGSpiceResult`: Structured simulation results
"""
function simulate_transfer_curve(spice_file::String; 
                                device_type::DeviceType=NMOS,
                                output_file::String="",
                                working_dir::String="")
    
    # Determine working directory if not provided
    if isempty(working_dir)
        working_dir = dirname(spice_file)
        if isempty(working_dir) || working_dir == "."
            working_dir = pwd()
        end
    end
    
    # Resolve paths - make sure we don't duplicate directory information
    # For the spice file, use just the filename if it's in the working directory
    base_spice_file = basename(spice_file)
    
    # Create full paths based on whether the input paths are absolute
    resolved_spice_file = isabspath(spice_file) ? spice_file : joinpath(working_dir, base_spice_file)
    resolved_output_file = isempty(output_file) ? "" : 
                          (isabspath(output_file) ? output_file : joinpath(working_dir, basename(output_file)))
    
    # Create a parameters dictionary to store the working directory
    parameters = Dict{String, Any}("working_dir" => working_dir)
    
    # For transfer curves, we want to look for vgs_id_data.txt by default, since that's what the
    # simulate_vgs_sweep.sp produces
    if isempty(resolved_output_file)
        resolved_output_file = joinpath(working_dir, "vgs_id_data.txt")
    end
    
    simulation = create_simulation(
        spice_file=resolved_spice_file,
        output_file=resolved_output_file,
        simulation_type=TRANSFER_CURVE,
        device_type=device_type,
        parameters=parameters
    )
    
    success = run_simulation(simulation, working_dir=working_dir)
    if !success
        error("Failed to run transfer curve simulation")
    end
    
    return read_simulation_output(simulation)
end

end # module