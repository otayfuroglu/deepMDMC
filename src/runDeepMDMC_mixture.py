import sys
import os
import numpy as np
import torch
torch.set_num_threads(6)
from ase import Atoms
from ase.io import read, write
from ase.data import vdw_radii
from ase.build import make_supercell
from ase.optimize import BFGS, LBFGS
from ase.filters import UnitCellFilter
from ase.io.trajectory import Trajectory

from time import time
from deep_mdmc import DeepMDMC
from molmod.units import *
from time import time

from utilities import PREOS

from multiprocessing import Pool
import argparse

from utilities import calculate_fugacity_with_coolprop

def getBoolStr(string):
    string = string.lower()
    if "true" in string or "yes" in string:
        return True
    elif "false" in string or "no" in string:
        return False
    else:
        print("%s is bad input!!! Must be Yes/No or True/False" %string)
        sys.exit(1)


parser = argparse.ArgumentParser(description="Give something ...")
parser.add_argument("-sim_type", type=str, required=True, help="")
parser.add_argument("-pressure", type=float, required=True, help="")
parser.add_argument("-temperature", type=float, required=True, help="")
parser.add_argument("-timestep", type=float, required=True, help="")
parser.add_argument("-totalsteps", type=int, required=True, help="")
parser.add_argument("-nmdsteps", type=int, required=True, help="")
parser.add_argument("-nmcswap", type=int, required=True, help="")
parser.add_argument("-nmcmoves", type=int, required=True, help="")
parser.add_argument("-flex_ads", type=str, required=True, help="")
parser.add_argument("-opt", type=str, required=True, help="")
parser.add_argument("-model_gcmc_path", type=str, required=True, help="")
parser.add_argument("-model_md_path", type=str, required=True, help="")
parser.add_argument("-struc_path", type=str, required=True, help="")
parser.add_argument("-molecule_path", type=str, required=True, help="")
parser.add_argument("-molecule2_path", type=str, required=False, help="")
parser.add_argument("-interval", type=int, required=True, help="")
parser.add_argument("-binary_ratio", type=float, required=False, help="")
args = parser.parse_args()

torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False

#  torch.set_default_dtype(torch.float32)

#  from nequip.ase.nequip_calculator import NequIPCalculator
#
#  def load_model(model_path):
#      #  Modify species for Mg-MOF-74 (see training yaml file)
#      return NequIPCalculator.from_deployed_model(
#          model_path = model_path, #'./MgF2_nonbonded_v10_nnp1_e10.pth',
#          species_to_type_name = {"C" : "C", "H" : "H", "O" : "O", "Mg" : "Mg"}, device=device,
#          #  set_global_options=True
#      )


#  temperature = 273 * kelvin
sim_type = args.sim_type
temperature = args.temperature * kelvin
pressure = args.pressure * bar
timestep = args.timestep
totalsteps = args.totalsteps
nmdsteps = args.nmdsteps # invoke this fix every nmdsteps steps
nmcswap = args.nmcswap # average number of GCMC exchanges to attempt every nmdsteps steps
nmcmoves = args.nmcmoves # average number of MC moves to attempt every nmdsteps steps
model_gcmc_path = args.model_gcmc_path
model_md_path = args.model_md_path
struc_path = args.struc_path
molecule_path = args.molecule_path
molecule2_path = args.molecule2_path
interval = args.interval
binary_ratio = args.binary_ratio

flex_ads = getBoolStr(args.flex_ads)
opt = getBoolStr(args.opt)
# Preferably run on GPUs
device = 'cuda'

#  calc_gcmc = load_model(model_gcmc_path)
#  calc_md = load_model(model_md_path)


vdw_radii = vdw_radii.copy()
#  vdw_radii[1] = 1.0
#  vdw_radii[6] = 1.0
#  vdw_radii[8] = 1.25
vdw_radii[12] = 1.25

#  fugacity = calculate_fugacity_with_coolprop("HEOS", "CO2", temperature, pressure)
eos_co2 = PREOS.from_name('carbondioxide')
eos_ch4 = PREOS.from_name('methane')
fugacity_co2 = eos_co2.calculate_fugacity(temperature, pressure)
fugacity_ch4 = eos_ch4.calculate_fugacity(temperature, pressure)

flex_ads_keyword = ""
if flex_ads:
    flex_ads_keyword = "_flexAds"
results_dir = f"{sim_type}_results_N{nmdsteps}_X{nmcswap+nmcmoves}{flex_ads_keyword}_{pressure/bar}bar_{int(temperature)}K"
if not os.path.exists(results_dir):
    os.mkdir(results_dir)


#  atom_type_pairs = {"Mg": 1, "O": 2,  "C": 3, "H": 4}
#  atom_type_pairs = {"Mg": 1, "O": 2,  "C": 3, "H": 4, "C2": 5, "O2": 6}
#  atom_type_pairs_frame = {"Mg": [1, 24.3050], "O": [2, 15.9994],  "C": [3, 12.0107], "H": [4, 1.00794]}
#  atom_type_pairs_ads = {"C": [5, 12.0107], "O": [6, 15.9994]}
#  atom_type_pairs_ads = {"C1": [5, 12.0107], "O1": [6, 15.9994], "O2": [7, 15.9994], "H2": [8, 1.00794]}
#  print(atom_type_pairs_ads.keys())

atom_type_pairs_frame = {1: ["Mg", 24.3050], 2: ["O", 15.9994],  3: ["C", 12.0107], 4: ["H", 1.00794]}
atom_type_pairs_ads = {5: ["C", 12.0107], 6: ["O", 15.9994], 7: ["C", 12.0107], 8: ["H", 1.00794]}
#  print([value[0] for value in list(atom_type_pairs_frame.values()) if value[0] == "C"])
#  atom_type_pairs_ads = {"C": [5, 12.0107], "O": [6, 15.9994]}
#  masses = {
#      1: 24.3050,
#      2: 15.9994,
#      3: 12.0107,
#      4: 1.00794,
#      5: 12.0107,
#      6: 15.9994,
#  }
#  specorder=["Mg", "O", "C", "H"]

tdump = 500 * timestep
pdump = 5000 * timestep

if sim_type != "tmmcmd":
    atoms_frame = read(struc_path)
    replica = [1, 1, 1]
    P = [[0, 0, -replica[0]], [0, -replica[1], 0], [-replica[2], 0, 0]]
    atoms_frame = make_supercell(atoms_frame, P)
    #  if opt:
    #      atoms_frame.calc = calc_md
    #      #  write("framebeforeopt.cif", atoms_frame)
    #      # geom opt frame based on model
    #      # to account relaxation of cell
    #      #  traj = Trajectory('frames.traj', 'w', atoms_frame)
    #      ucf = UnitCellFilter(atoms_frame)
    #      optimizer = LBFGS(ucf)
    #      #  optimizer.attach(traj)
    #      optimizer.run(fmax=0.001)
    write("frame0.extxyz", atoms_frame)
    #  write("frame0.cif", atoms_frame)
    #  atoms_frame = read("frame0.cif")
    #  quit()
    # C and O were renamed to Co and Os to differentiate them from framework atoms during training
    #  atoms_ads = read('./co2_v2.xyz')
    atoms_ch4 = read(molecule_path)
    atoms_co2 = read(molecule2_path)
    Z_ch4 = 0
    Z_co2 = 0
    deep_mdmc = DeepMDMC(model_gcmc_path, model_md_path, results_dir, interval, atoms_frame,
                         temperature, pressure, device, vdw_radii)

# initialize MD

#  if nmdsteps == 0:
    #  sim_type = "rigid"


if sim_type.lower() == "binary_rigid":
    if binary_ratio is None:
        print("Error: binary ratio is none")
        sys.exit(1)
    deep_mdmc.init_gcmc(atoms_ch4, fugacity_ch4, Z_ch4, flex_ads,)
    deep_mdmc.init_gcmc_binary( atoms_ch4, fugacity_ch4, Z_ch4, atoms_co2, fugacity_co2, Z_co2, flex_ads,)
    deep_mdmc.run_gcmc_binary(nmcswap, nmcmoves, binary_ratio)
elif sim_type.lower() == "rigid":
    deep_mdmc.init_gcmc(atoms_ch4, fugacity_ch4, Z_ch4, flex_ads,)
    #  deep_mdmc.init_gcmc(atoms_co2, fugacity_co2, Z_co2, flex_ads,)
    deep_mdmc.run_gcmc(nmcswap, nmcmoves)
elif sim_type.lower() == "binary_gcmcmd":
    if binary_ratio is None:
        print("Error: binary ratio is none")
        sys.exit(1)
    deep_mdmc.init_md(timestep, atom_type_pairs_frame,
                      atom_type_pairs_ads, units_lmp="metal",
                      tdump=tdump, pdump=pdump, md_type="npt",
                      opt=True, equ_steps=50000)
                      #  opt=True, equ_steps=100)
    #  deep_mdmc.init_gcmc(atoms_ads, Z_ads, flex_ads)
    deep_mdmc.init_gcmc_binary( atoms_ch4, fugacity_ch4, Z_ch4, atoms_co2,
                               fugacity_co2, Z_co2, flex_ads,)
    deep_mdmc.run_gcmcmd_binary(totalsteps, nmdsteps, nmcswap, nmcmoves, binary_ratio)
elif sim_type.lower() == "gcmcmd":
    deep_mdmc.init_md(timestep, atom_type_pairs_frame,
                      atom_type_pairs_ads, units_lmp="metal",
                      tdump=tdump, pdump=pdump, md_type="npt",
                      opt=True, equ_steps=50000)
                      #  opt=True, equ_steps=100)
    #  deep_mdmc.init_gcmc(atoms_ads, Z_ads, flex_ads)
    deep_mdmc.init_gcmc(atoms_ch4, fugacity_ch4, Z_ch4, flex_ads,)
    #  deep_mdmc.init_gcmc(atoms_co2, fugacity_co2, Z_co2, flex_ads,)
    deep_mdmc.run_gcmcmd(totalsteps, nmdsteps, nmcswap, nmcmoves)
elif sim_type.lower() == "tmmcmd":
    #  pass
    #  atoms_frame0 = read(struc_path)
    #  replica = [1, 1, 1]
    #  P = [[0, 0, -replica[0]], [0, -replica[1], 0], [-replica[2], 0, 0]]
    #  atoms_frame = make_supercell(atoms_frame, P)
    #  write("frame0.extxyz", atoms_frame)

    atoms_ads = read(molecule_path)
    vdw_radii = vdw_radii.copy()
    # Mg radius is set to 1.0 A
    vdw_radii[12] = 1.0
    atoms_frame0 = read("frame0.extxyz")
    atoms_frame = read("loaded_frame1.extxyz")
    Z_ads = int((len(atoms_frame) - len(atoms_frame0)) / len(atoms_ads))
    #  print(N_ads)
    #  quit()
    deep_mdmc = DeepMDMC(model_gcmc_path, model_md_path, results_dir, interval, atoms_frame, atoms_ads, flex_ads,
                      Z_ads, temperature, pressure, fugacity, device, vdw_radii)

    #  deep_mdmc.init_gcmc()
    deep_mdmc.init_md(timestep, atom_type_pairs_frame,
                      atom_type_pairs_ads, units_lmp="metal",
                      tdump=tdump, pdump=pdump, md_type="npt",
                      opt=False, equ_steps=0)
    deep_mdmc.run_tmmcmd(nmdsteps=180000)
else:
    print("Error: invalid sim_type")



