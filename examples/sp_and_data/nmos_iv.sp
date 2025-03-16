* 65nm NMOS IV Curve Simulation
* Include the PTM model
.include 65nm_bulk.pm

* Circuit Definition
* Format: M{name} {drain} {gate} {source} {body} {model} W={width} L={length}
M1 d g s b nmos W=1u L=65n

* DC Voltage Sources
VDS d 0 0.05V     ; Drain-Source voltage
VGS g 0 0V        ; Gate-Source voltage
VS  s 0 0V        ; Source tied to ground
VB  b 0 0V        ; Body tied to ground

* DC Sweep Analysis for NMOS Id vs Vds (Output characteristics)
* Sweep VDS from 0 to 1.2V at different VGS values
.dc VDS 0 1.2 0.01 VGS 0.2 1.2 0.2

* Control Options
.options temp=27

* Output Control
.control
run
echo "IV Curves for NMOS with varying Vgs values"

* Save data in text format - this is the key part for easy parsing
set wr_vecnames          * Include variable names as headers
set wr_singlescale       * Use a single scale for all variables
wrdata nmos_iv_data.txt v(d) i(VDS)  * Save specific vectors to a text file

* Also create a plot (only works in interactive mode)
plot -i(VDS) vs v(d) title 'NMOS Output Characteristic' ylabel 'Drain Current (A)' xlabel 'Drain Voltage (V)'
echo "Data saved to nmos_iv_data.txt"
.endc

.end