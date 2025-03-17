using NGSpiceTools
using Plots
using Printf
"""
This script demonstrates how to parse and plot IV curves from NGSpice data
that contains multiple VGS sweeps. It uses the specific functions designed 
to handle the multi-block format that ngspice produces.
"""

simulate_vds_vgs_ids("./examples/sp_and_data/nmos_iv.sp")


# Path to the data file
data_file = "./examples/sp_and_data/nmos_iv_data.txt"

# Simple approach - parse and plot in one step
println("Parsing and plotting IV curves...")
plt = parse_and_plot_iv_curves(data_file, 
                              title="IV Characteristics (65nm PTM)",
                              is_nmos=true)

# Save the figure
savefig(plt, "nmos_iv_family_curves.png")
println("Plot saved to nmos_iv_family_curves.png")

# Display some information about the IV curves
data = parse_ngspice_iv_curves(data_file, is_nmos=true)
vgs_values = data["vgs_values"]
ids_matrix = data["ids"]

println("\nMOSFET IV Curve Analysis")
println("------------------------")
println("Device type: NMOS")
println("Number of VGS values: $(length(vgs_values))")
println("VGS values: $vgs_values")

# Calculate maximum current for each VGS
println("\nVGS (V) | Max IDS (A)")
println("--------|------------")
for i in 1:length(vgs_values)
    ids_column = ids_matrix[:, i]
    max_current = maximum(filter(!isnan, ids_column))
    println(@sprintf("%.1f      | %.4e", vgs_values[i], max_current))
end

# More advanced usage - if you want to customize the plot further
plt2 = plot_iv_family_curves(data, title="NMOS IV Characteristics (65nm PTM)")
xlabel!("Drain Voltage (V)")
ylabel!("Drain Current (A)")
plot!(size=(800, 600))  # Adjust plot size
plot!(dpi=300)  # High resolution for publishing
savefig(plt2, "nmos_iv_curves_advanced.png")
println("\nEnhanced plot saved to nmos_iv_curves_advanced.png")

data_file_pmos = "./examples/sp_and_data/pmos_iv_data.txt"
plt_pmos = parse_and_plot_iv_curves(data_file_pmos, 
                              title="PMOS IV Characteristics (65nm PTM)",
                              is_nmos=false)

#%% Vgs ID sweep

simulate_vgs_sweep("./examples/sp_and_data/simulate_vgs_sweep.sp")

vgs_id_dict = extract_vgs_id_data("./examples/sp_and_data/vgs_id_data.txt")
plt = parse_and_plot_vgs_id("./examples/sp_and_data/vgs_id_data.txt")
