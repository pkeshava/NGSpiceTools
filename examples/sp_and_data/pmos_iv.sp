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