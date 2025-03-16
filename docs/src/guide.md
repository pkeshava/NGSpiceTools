# User Guide

## Sign Conventions

NGSpiceTools.jl applies standard device physics sign conventions:

- **NMOS devices**: 
  - Drain voltage (VDS) is shown as positive
  - Drain current (ID) is shown as positive

- **PMOS devices**: 
  - Drain-source voltage (VDS) is shown as negative
  - Drain current (ID) is shown as positive

## Basic Usage

### Plotting IV Curves from NGSpice Output

```julia
using NGSpiceTools
using Plots

# Parse a log file
log_result = parse_ngspice_log("nmos_results.log")

# Plot IV curve from log file and its associated data
plt = simulate_iv_from_log("nmos_results.log")
title!("NMOS IV Characteristics")
xlabel!("Drain Voltage (V)")
ylabel!("Drain Current (A)")
savefig(plt, "nmos_iv.png")
```

### Plotting Multiple IV Curves (Family)

```julia
using NGSpiceTools
using Plots

# Parse and plot a family of IV curves
plt = parse_and_plot_iv_curves("nmos_iv_data.txt", 
                              title="NMOS IV Characteristics (65nm PTM)",
                              is_nmos=true)
savefig(plt, "nmos_iv_family.png")
```

## NMOS Example

Here's a complete example workflow for analyzing NMOS devices:

1. First, set up your NGSpice script with data export:

```text
* 65nm NMOS IV Curve Simulation
.include 65nm_bulk.pm

* Circuit Definition
M1 d g s b nmos W=1u L=65n

* DC Voltage Sources
VDS d 0 0.05V
VGS g 0 0V
VS  s 0 0V
VB  b 0 0V

* DC Sweep Analysis
.dc VDS 0 1.2 0.01 VGS 0.2 1.2 0.2

* Output Control
.control
run
set wr_vecnames
set wr_singlescale
wrdata nmos_iv_data.txt v(d) i(VDS)
.endc

.end
```

2. Run the simulation in NGSpice:

```bash
julia --project
using NGSpiceTools
; cd examples/sp_and_data
ngspice -b nmos_script.sp -o nmos_results.log
```

3. Analyze and plot the results in Julia:

```julia
using NGSpiceTools
using Plots

# Plot the IV curves
plt = parse_and_plot_iv_curves("nmos_iv_data.txt", 
                             title="NMOS IV Characteristics (65nm PTM)",
                             is_nmos=true)

# Customize the plot
xlabel!("Drain Voltage (V)")
ylabel!("Drain Current (A)")

# Save the figure
savefig(plt, "nmos_iv_family_curves.png")
```

## PMOS Example

```bash
julia --project
using NGSpiceTools
; cd examples/sp_and_data
ngspice -b pmos_script.sp -o nmos_results.log
```

For PMOS devices, the process is similar but with appropriate sign conventions:

1. NGSpice script:

```text
* 65nm PMOS IV Curve Simulation
.include 65nm_bulk.pm

* Circuit Definition
M1 d g s b pmos W=1u L=65n

* DC Voltage Sources
VDS d s 0.05V
VGS g s 0V
VS  s 0 0V
VB  b 0 0V

* DC Sweep Analysis
.dc VDS 0 -1.2 -0.01 VGS 0 -1.2 -0.2

* Output Control
.control
run
set wr_vecnames
set wr_singlescale
wrdata pmos_iv_data.txt v(d,s) i(VDS)
.endc

.end
```

2. Analysis:

```julia
using NGSpiceTools
using Plots

# Plot the IV curves with is_nmos=false for PMOS
plt = parse_and_plot_iv_curves("pmos_iv_data.txt", 
                             title="PMOS IV Characteristics (65nm PTM)",
                             is_nmos=false)

# Customize the plot
xlabel!("Drain-Source Voltage (V)")
ylabel!("Drain Current (A)")

# Save the figure
savefig(plt, "pmos_iv_family_curves.png")
```

## Advanced Usage

### Custom IV Curve Analysis

For more control over the analysis process:

```julia
using NGSpiceTools
using Plots

# Parse the IV curve data
data = parse_ngspice_iv_curves("nmos_iv_data.txt", is_nmos=true)

# Access the data
vgs_values = data["vgs_values"]
vds_matrix = data["vds"]
ids_matrix = data["ids"]

# Perform custom analysis
max_currents = [maximum(filter(!isnan, ids_matrix[:, i])) for i in 1:length(vgs_values)]

# Custom plot
plt = plot(vgs_values, max_currents, 
          title="Maximum Current vs. Gate Voltage",
          xlabel="Gate Voltage (V)",
          ylabel="Maximum Drain Current (A)",
          marker=:circle,
          linewidth=2)

savefig(plt, "nmos_max_current.png")
```

## Troubleshooting

### Common Issues

1. **Data files not found**: Make sure to check that your data files are in the right location and have the expected format.

2. **Incorrect sign conventions**: If your plots show currents or voltages with incorrect signs, make sure to set the `is_nmos` parameter appropriately.

3. **Missing headers in data files**: For best results, always include the headers in your NGSpice output by using the `set wr_vecnames` command.

### Debugging

To get more information about the parsed data:

```julia
data = parse_ngspice_iv_curves("nmos_iv_data.txt", is_nmos=true)
println("VGS values: ", data["vgs_values"])
println("Data dimensions: ", size(data["vds"]), " ", size(data["ids"]))
```