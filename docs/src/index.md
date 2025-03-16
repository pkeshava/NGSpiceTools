# NGSpiceTools.jl

A Julia package for processing and visualizing NGSpice simulation output.

## Overview

This package provides tools to read and visualize NGSpice simulation results, with a particular focus on MOSFET IV characteristics. It includes functions for parsing NGSpice output files, applying appropriate sign conventions, and creating plots.

## Features

- Parse NGSpice log files
- Read data files from NGSpice simulations
- Generate IV curve plots with appropriate axis labeling
- Apply correct sign conventions for NMOS and PMOS devices
- Automatic device type detection

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yourusername/NGSpiceTools.jl")
```

## Quick Start

```julia
using NGSpiceTools
using Plots

# Plot IV curves from a log file and its associated data file
plt = simulate_iv_from_log("nmos_results.log")
savefig(plt, "nmos_iv_curves.png")
```

## NGSpice Data Export

When running NGSpice simulations, make sure to add these commands to your `.control` section to save data in an easily parsable text format:

```spice
* Save data in text format
set wr_vecnames          * Include variable names as headers
set wr_singlescale       * Use a single scale for all variables
wrdata your_filename.txt v(d) i(VDS)  * Save specific vectors to a text file
```

## Sign Conventions

This package applies standard device physics sign conventions:

- **NMOS devices**: 
  - Drain voltage (VDS) is shown as positive
  - Drain current (ID) is shown as positive

- **PMOS devices**: 
  - Drain-source voltage (VDS) is shown as negative
  - Drain current (ID) is shown as positive

## Example: NMOS IV Curves

```julia
using NGSpiceTools
using Plots

# NMOS device (explicitly set is_nmos=true)
plt = simulate_iv_from_log("nmos_results.log", is_nmos=true)
title!("NMOS IV Characteristics (65nm PTM)")
xlabel!("Drain Voltage (V)")
ylabel!("Drain Current (A)")
savefig(plt, "nmos_iv_curves.png")
```

## Example: PMOS IV Curves

```julia
using NGSpiceTools
using Plots

# PMOS device (set is_nmos=false)
plt = simulate_iv_from_log("pmos_results.log", is_nmos=false)
title!("PMOS IV Characteristics (65nm PTM)")
xlabel!("Drain-Source Voltage (V)")
ylabel!("Drain Current (A)")
savefig(plt, "pmos_iv_curves.png")
```

## Package Functions

### Main Functions

- `parse_ngspice_log(logfile)`: Parse an NGSpice log file
- `read_iv_data(datafile; is_nmos=true)`: Read IV curve data from a text file
- `simulate_iv_from_log(logfile; datafile="", is_nmos=true)`: Plot IV curves from a log file

### Utilities

- `parse_simulation_type(lines)`: Determine the simulation type from log file lines

## License

This package is released under the MIT License.