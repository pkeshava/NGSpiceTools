using NGSpiceTools

# Define the data directory where your files are located
data_dir = "./examples/sp_and_data"

println("Running NMOS IV curve simulation...")
# Just use the filename, not the full path, when specifying the working directory
nmos_result = simulate_iv_curves("nmos_iv.sp", 
                               device_type=NMOS, 
                               working_dir=data_dir)

# Plot the IV curves
nmos_plot = plot_iv_family_curves(nmos_result, title="NMOS Output Characteristics")
display(nmos_plot)

# Extract device parameters
nmos_params = extract_parameters(nmos_result)
if haskey(nmos_params, "output_resistance")
    println("NMOS Output Resistance Values: ", nmos_params["output_resistance"])
end

println("\nRunning PMOS IV curve simulation...")
pmos_result = simulate_iv_curves("pmos_iv.sp", 
                               device_type=PMOS, 
                               working_dir=data_dir)

# Plot the IV curves
pmos_plot = plot_iv_family_curves(pmos_result, title="PMOS Output Characteristics")

println("\nRunning NMOS transfer curve simulation...")
# For the transfer curve, specify the exact output filename we expect
transfer_result = simulate_transfer_curve("simulate_vgs_sweep.sp", 
                                        device_type=NMOS,
                                        output_file="vgs_id_data.txt",
                                        working_dir=data_dir)

# Plot linear and log scale transfer curves
linear_plot = plot_transfer_curve(transfer_result, title="NMOS Transfer Characteristic")

log_plot = plot_transfer_curve(transfer_result, title="NMOS Transfer Characteristic (Log Scale)", scale=:log)
display(log_plot)

# Extract threshold voltage
vt_params = extract_parameters(transfer_result)
if haskey(vt_params, "threshold_voltage")
    println("Threshold Voltage: ", vt_params["threshold_voltage"], " V")
end
if haskey(vt_params, "max_transconductance")
    println("Maximum Transconductance: ", vt_params["max_transconductance"], " S")
end

println("\nAlternatively, read existing transfer curve data without running a new simulation:")
# To read an existing data file without running a simulation
existing_data_file = joinpath(data_dir, "vgs_id_data.txt")
existing_simulation = create_simulation(
    spice_file=joinpath(data_dir, "simulate_vgs_sweep.sp"),
    output_file=existing_data_file,
    simulation_type=TRANSFER_CURVE,
    device_type=NMOS
)
existing_result = read_simulation_output(existing_simulation)
existing_plot = plot_transfer_curve(existing_result, title="Transfer Curve from Existing Data")
