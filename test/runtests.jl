using NGSpiceTools
using Test

@testset "NGSpiceTools.jl" begin
    # Test the parsing of simulation type
    @testset "Simulation Type Detection" begin
        lines = [
            "Circuit: * Test circuit",
            ".dc VDS 0 1.2 0.01 VGS 0.2 1.2 0.2",
            "No. of Data Rows : 726",
            "IV Curves for NMOS with varying Vgs values"
        ]
        @test parse_simulation_type(lines) == "dc"
    end
    
    # Test IV curve data parsing
    @testset "IV Data Parsing" begin
        # Create a temp data file for testing
        test_data = """
         v-sweep         v(d)            i(VDS)         
         0.00000000e+00  0.00000000e+00  4.21287615e-12 
         1.00000000e-02  1.00000000e-02 -1.06350992e-07 
        """
        
        test_file = tempname()
        open(test_file, "w") do io
            write(io, test_data)
        end
        
        # Test the parsing function
        result = read_iv_data(test_file, is_nmos=true)
        @test haskey(result, "v")
        @test haskey(result, "i")
        
        # Test family curves parsing
        if isfile("../examples/nmos_iv_data.txt")
            data = parse_ngspice_iv_curves("../examples/nmos_iv_data.txt")
            @test haskey(data, "vgs_values")
            @test haskey(data, "vds")
            @test haskey(data, "ids")
        end
        
        # Clean up
        rm(test_file)
    end
end