* NGSPICE simulation of Vgs vs. Id sweep using PTM 65nm MOS model

.include '65nm_bulk.pm'

* Transistor definition
M1 drain gate source bulk nmos W=1u L=65n

* Voltage sources
Vds drain source DC 1.2
Vgs gate source DC 0
Vbs bulk source DC 0

* DC sweep of Vgs from 0 to 1.2V in 0.01V steps
.dc Vgs 0 1.2 0.01

* Save Id data
.control
  run
  set wr_vecnames
  set wr_singlescale
  wrdata ./examples/sp_and_data/simulate_vgs_sweep.txt v(gate) i(Vds)
  quit
.endc

.end
