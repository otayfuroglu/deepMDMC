#! /usr/bin/bas

SCRIPT_DIR=/path/to/mdmc
model_md_path=/path/to/nequip_model/
model_gcmc_path=/path/to/nequip_model/
struc_path=/path/to/mof_structure
molecule_path=//path/to/guest_molecule
T=298
P=1.0

nmdsteps=0
nmcswap=250000
nmcmoves=250000
timestep=0.0005 # ps


python $SCRIPT_DIR/runDeepMDMC.py\
       	-sim_type rigid\
       	-pressure $P\
       	-temperature $T\
       	-timestep $timestep\
       	-totalsteps 1000000\
       	-nmdsteps $nmdsteps\
       	-nmcswap $nmcswap\
       	-nmcmoves $nmcmoves\
       	-flex_ads no\
       	-opt no\
	-interval 50\
       	-model_gcmc_path $model_gcmc_path\
       	-model_md_path $model_md_path\
       	-struc_path $struc_path\
       	-molecule_path $molecule_path 
