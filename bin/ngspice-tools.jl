#!/usr/bin/env julia

using NGSpiceTools
using Plots
using ArgParse

"""
    parse_commandline()

Parse command line arguments for the NGSpice IV curve plotting script.
"""
function parse_commandline()
    s = ArgParseSettings(
        description="Plot IV curves from NGSpice simulation results",
        prog="plot_iv_curves.jl"
    )
    
    @add_arg_table! s begin
        "--log", "-l"
            help = "NGSpice log file"
            arg_type = String
            required = true
        "--data", "-d"
            help = "Data file (optional, will be auto-detected if not provided)"
            arg_type = String
            default = ""
        "--output", "-o"
            help = "Output image file"
            arg_type = String
            default = "iv_curves.png"
        "--device", "-t"
            help = "Device type (nmos or pmos)"
            arg_type = String
            default = "auto"
        "--title"
            help = "Plot title"
            arg_type = String
            default = "IV Characteristics"
        "--width", "-w"
            help = "Plot width in pixels"
            arg_type = Int
            default = 800
        "--height", "-h"
            help = "Plot height in pixels"
            arg_type = Int
            default = 600
    end
    
    return parse_args(s)
end

"""
    main()

Main function to plot IV curves from NGSpice simulation.
"""
function main()
    # Parse command line arguments
    args = parse_commandline()
    
    # Determine if device is NMOS or PMOS
    is_nmos = if args["device"] == "auto"
        # Auto-detect from filename
        !occursin(r"pmos|pfet|pchan", lowercase(args["log"]))
    else
        args["device"] == "nmos"
    end
    
    # Print information
    device_type = is_nmos ? "NMOS" : "PMOS"
    println("Processing $(args["log"]) as $device_type device")
    
    # Try to plot the IV curves
    try
        # Generate the plot
        plt = simulate_iv_from_log(
            args["log"], 
            datafile=args["data"],
            is_nmos=is_nmos
        )
        
        # Customize plot
        title!(plt, args["title"])
        xlabel!(plt, is_nmos ? "Drain Voltage (V)" : "Drain-Source Voltage (V)")
        ylabel!(plt, "Drain Current (A)")
        
        # Save the plot
        savefig(plt, args["output"])
        println("Plot saved to $(args["output"])")
    catch e
        # Handle errors
        println("Error: $(e.msg)")
        if occursin("No usable data file found", e.msg)
            println("\nHint: Make sure your NGSpice script includes data export commands:")
            println("set wr_vecnames")
            println("set wr_singlescale")
            println("wrdata filename.txt v(d) i(VDS)")
        end
        exit(1)
    end
end

# Run the main function if this script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end