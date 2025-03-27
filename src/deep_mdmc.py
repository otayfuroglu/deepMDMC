import numpy as np
import random

from ase import Atoms
from ase.cell import Cell
from ase.io import read, write
from molmod.units import *
from molmod.constants import *
from molmod.periodic import periodic
from utilities import _random_rotation, random_position, vdw_overlap, vdw_collision
import tqdm

#  from ase.md.langevin import Langevin
#  from ase.md.npt import NPT
#  from ase.md.nptberendsen import NPTBerendsen, Inhomogeneous_NPTBerendsen
#
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import fs, kB
from ase.io.trajectory import Trajectory

from ase.io.lammpsdata import read_lammps_data
from lammpsdata import write_lammps_data # modified 

from itertools import combinations

from lammps import lammps
from ase.calculators.lammps import Prism, convert



class DeepMDMC():
    def __init__(self, model_path_gcmc, model_path_md, results_dir, interval, atoms_frame, atoms_ads, flex_ads, Z_ads, T, P, fugacity, device, vdw_radii):
        self.model_path_gcmc = model_path_gcmc
        #  self.calc_md = calc_md
        self.model_path_md = model_path_md
        self.results_dir = results_dir
        self.interval = interval
        self.atoms = atoms_frame
        self.n_frame = len(atoms_frame) - 3 * Z_ads
        self.atoms_ads = atoms_ads
        self.flex_ads = flex_ads
        self.n_ads = len(self.atoms_ads)
        self.cell = np.array(atoms_frame.get_cell())
        self.V = np.linalg.det(self.cell) * angstrom**3
        self.T = T
        self.P = P
        self.fugacity = fugacity
        self.device = device
        self.beta = 1 / (boltzmann * T)
        self.Z_ads = Z_ads

        self.vdw = vdw_radii - 0.35
        self.lmp = None

    def _get_ads_atoms(self, atoms):

        ads_atoms = atoms[self.n_frame:]

        ads_molecules = []
        for i in range(0, len(ads_atoms), self.n_ads):
            ads_molecule = Atoms(
                symbols = ads_atoms[i:i+self.n_ads].get_chemical_symbols(),
                positions = ads_atoms[i:i+self.n_ads].get_positions(),
                cell = atoms.get_cell()
            )
            ads_molecules.append(ads_molecule)

        # Combine all CO₂ molecules into a single Atoms object
        #  ads_system = sum(ads_molecules, start=Atoms())
        #  write("test_ads_molecules.extxyz", ads_system)
        assert len(ads_molecules) == self.Z_ads

        return ads_molecules

    def _set_rigid_ads_atoms(self):

        from ase.constraints import FixBondLengths

        indices = [atom.index for atom in self.atoms]
        ads_indices = indices[self.n_frame:]

        assert len(ads_indices)/3 == self.Z_ads

        for ad_indices in self._split_ads_indices(ads_indices, self.n_ads):
            fix_bond_indices = list(combinations(ad_indices, 2))[:-1]
            c = FixBondLengths(fix_bond_indices)
            self.atoms.set_constraint(c)

    def _set_rigid_triatomic_ads_atoms(self):

        from ase.constraints import FixLinearTriatomic

        indices = [atom.index for atom in self.atoms]
        ads_indices = indices[self.n_frame:]

        assert len(ads_indices)/3 == self.Z_ads

        if self.Z_ads != 0:
            sub_ads_indices = self._split_ads_indices(ads_indices, self.n_ads)
            c = FixLinearTriatomic(triples=sub_ads_indices)
            self.atoms.set_constraint(c)

    def _split_ads_indices(self, lst, length):
        """Split a list into sub tuple of the specified length."""

        # Return a new tuple with the smallest value in the middle
        return [(sorted(lst[i:i + length])[1],
                 sorted(lst[i:i + length])[0],
                 sorted(lst[i:i + length])[2])
                    for i in range(0, len(lst), length)
               ]

    def _get_e_ads(self, atoms):
        e_ads = 0.0
        for ads_atoms in self._get_ads_atoms(atoms):
            ads_atoms.calc = self.calc_gcmc
            e_ads += ads_atoms.get_potential_energy()
        return e_ads

    def _setTqdm(self, start_step, totalsteps):
        self.pbar = tqdm.tqdm(total=totalsteps)
        self.pbar.n = start_step
        self.pbar.last_print_n = start_step  # Ensure the display updates correctly
        self.pbar.refresh()

    def _tqdmMD(self):
        self.pbar.update()

    def _atoms2lammpsdata(self, data_file, atoms, atom_type_pairs):
        #  print(atoms.cell)
        #  data_file = "data.lmp"
        write_lammps_data(data_file, atoms, units=self.units_lmp,
                          specorder=atom_type_pairs,
                          velocities=False,
                          atom_style="full"
                         )
        #  write("test_trj.extxyz", atoms, append=True)
        #NOTE: there is some problem velocity unit convertsioin with write_lammps_data.
        # we write manual below 
        with open(data_file, 'a') as f:
            f.write("\n\nVelocities\n\n")
            for i, vel in enumerate(atoms.get_velocities(), start=1):
                f.write(f"{i} {vel[0]:.6f} {vel[1]:.6f} {vel[2]:.6f}\n")

    def _lammps2ase_atoms(self, atom_type_pairs):
        """Retrieve the last frame of the simulation and convert it to ASE Atoms."""
        natoms = self.lmp.get_natoms()
        atom_types = np.array(self.lmp.gather_atoms("type", 0, 1))  # Gather atom types
        positions = np.array(self.lmp.gather_atoms("x", 1, 3))  # Gather positions
        velocities = np.array(self.lmp.gather_atoms("v", 1, 3))  # Gather velocities
        box_data = self.lmp.extract_box()

        xlo, xhi = box_data[0][0], box_data[1][0]
        ylo, yhi = box_data[0][1], box_data[1][1]
        zlo, zhi = box_data[0][2], box_data[1][2]
        xy = box_data[2]
        xz = box_data[3]
        yz = box_data[4]

        # Construct the cell matrix
        cell = Cell(np.array([
            [xhi - xlo, 0.0, 0.0],
            [xy, yhi - ylo, 0.0],
            [xz, yz, zhi - zlo]
        ]))

        symbols = []
        for atom_type in atom_types:
            symbols += [symbol for symbol, atom_type_pair in atom_type_pairs.items() if atom_type_pair == atom_type]

        atoms = Atoms(
            symbols=symbols,
            positions=positions.reshape((natoms, 3)),
            cell=cell,
            pbc=True
        )
        atoms.set_velocities(velocities.reshape((natoms, 3)))
        return atoms

    def _ase_to_lammpstrj(self, atoms_frame, atoms_ads, filename, timestep=0, units="metal"):
        """
        Converts an ASE Atoms object to a LAMMPS trajectory (lammpstrj) file.

        Parameters:
            atoms (ase.Atoms): ASE Atoms object containing atom data.
            filename (str): Output LAMMPS trajectory file.
            timestep (int): Timestep to write in the trajectory file.
        """
        with open(filename, 'w') as f:
            # Write timestep
            f.write("ITEM: TIMESTEP\n")
            f.write(f"{timestep}\n")

            # Write number of atoms
            f.write("ITEM: NUMBER OF ATOMS\n")
            f.write(f"{len(atoms_frame + atoms_ads)}\n")

            # Write box bounds
            # NOTE info about lammps dumped file cell orginization
            """
                xlo_bound xhi_bound xy
                ylo_bound yhi_bound xz
                zlo_bound zhi_bound yz

                xlo_bound = xy
                ylo_bound = yz
                zlo_bound = xz

                xhi_bound = xhi
                yhi_bound = yhi
                zhi_bound = zhi

            """
            f.write("ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
            prismobj = Prism(atoms_frame.get_cell(), reduce_cell=False)

            #  Get cell parameters and convert from ASE units to LAMMPS units
            xhi, yhi, zhi, xy, xz, yz = convert(
                prismobj.get_lammps_prism(), 'distance', 'ASE', units)

            f.write(f"{xy:.16e} {xhi:.16e} {xy:.16e}\n")
            f.write(f"{yz:.16e} {yhi:.16e} {xz:.16e}\n")
            f.write(f"{xz:.16e} {zhi:.16e} {yz:.16e}\n")

            # Write atomic data for frame
            i = 0
            mol_type = 0
            f.write("ITEM: ATOMS id mol type x y z vx vy vz\n")
            positions = convert(atoms_frame.positions, 'distance', 'ASE', units)
            velocities = convert(atoms_frame.get_velocities(), 'velocity', 'ASE', units)
            for atom in atoms_frame:
                pos = positions[i]
                vel = velocities[i] if velocities is not None else [0.0, 0.0, 0.0]
                f.write(f"{i+1} {mol_type} {self.atom_type_pairs_frame[atom.symbol][0]} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f} {vel[0]:.8f} {vel[1]:.8f} {vel[2]:.8f}\n")
                i += 1

            mol_type += 1 # next molecule type 
            # Write atomic data for ads
            positions = convert(atoms_ads.positions, 'distance', 'ASE', units)
            velocities = convert(atoms_ads.get_velocities(), 'velocity', 'ASE', units)
            for j, atom in enumerate(atoms_ads):
                pos = positions[j]
                vel = velocities[j] if velocities is not None else [0.0, 0.0, 0.0]
                f.write(f"{i+1} {mol_type} {self.atom_type_pairs_ads[atom.symbol][0]} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f} {vel[0]:.8f} {vel[1]:.8f} {vel[2]:.8f}\n")
                if j != 0 and (j+1) % self.n_ads == 0:
                    mol_type += 1 # next molecule type
                i += 1

    def init_md(self, timestep, atom_type_pairs_frame, atom_type_pairs_ads, units_lmp, tdump, pdump, md_type, opt, equ_steps):

        self.units_lmp = units_lmp
        self.tdump = tdump
        self.pdump = pdump
        data_file = f"{self.results_dir}/data.lmp"
        self.atom_type_pairs_frame = atom_type_pairs_frame
        self.atom_type_pairs_ads = atom_type_pairs_ads

        specorder = list(self.atom_type_pairs_frame.keys()) + list(self.atom_type_pairs_ads.keys())

        MaxwellBoltzmannDistribution(self.atoms, self.T * 0.8)

        # NOTE _atoms2lammpsdata only for inintia frame not system with gas
        # becauase of the same atom symbol not extracted from atom_type_pairs
        self._atoms2lammpsdata(data_file, self.atoms, atom_type_pairs_frame)

        #  if not self.flex_ads:
            #  self._set_rigid_ads_atoms()
            #  self._set_rigid_triatomic_ads_atoms()

        # Initialize LAMMPS
        self.lmp = lammps()
        self.lmp.command("clear")
        self.lmp.command(f"units {self.units_lmp}")
        self.lmp.command("dimension 3")
        self.lmp.command("boundary p p p")
        #  self.lmp.command("atom_style atomic")
        self.lmp.command("atom_style full") # NOTE rigid/small not valid for atom_style atomic
        self.lmp.command("newton off")
        #  self.lmp.command(f"read_data {data_file}")
        self.lmp.command(f"read_data {data_file} extra/atom/types 2")
        self.lmp.command(f"molecule co2mol CO2.txt")

        for values in atom_type_pairs_frame.values():
           self.lmp.command(f"mass {values[0]} {values[1]}")
        for values in atom_type_pairs_ads.values():
           self.lmp.command(f"mass {values[0]} {values[1]}")


        self.lmp.command("pair_style nequip")
        self.lmp.command(f"pair_coeff * * {self.model_path_md} {' '.join(specorder)}")
        #  self.lmp.command(f"pair_coeff * * {self.model_path_md} Mg O C H C O")

        self.lmp.command(f"thermo {self.interval}")
        self.lmp.command("compute moltemp all temp")
        self.lmp.command("compute_modify moltemp dynamic/dof yes")
        self.lmp.command("compute_modify thermo_temp dynamic/dof yes")
        self.lmp.command("neighbor 2.0 bin")
        self.lmp.command("neigh_modify every 1 delay 10 check yes")
        #  self.lmp.command(f"velocity all create {self.T} 54654")

        self.lmp.command(f"log {self.results_dir}/md_{self.T}K_{self.P/bar}bar.log")
        self.lmp.command(f"timestep {timestep}")

        if opt:
            self.lmp.command("min_style cg")
            self.lmp.command("minimize 1.0e-8 1.0e-4 1000 10000")

        if equ_steps != 0:
            self.lmp.command(f"fix init_nvt all nvt temp {self.T} {self.T} {self.tdump}")
            self.lmp.command(f"run {equ_steps}")
            #  self.lmp.command(f"run 50")
            self.lmp.command(f"unfix init_nvt")
            self.lmp.command(f"reset_timestep 0")


        self.lmp.command(f"fix 1 all recenter INIT INIT NULL")
        self.lmp.command(f"dump 1 all atom {self.interval} {self.results_dir}/md_{self.T}K_{self.P/bar}bar.lammpstrj")
        self.lmp.command(f"dump_modify 1 append yes")

        if md_type == "npt":
            self.lmp.command(f"fix mynpt all npt temp {self.T} {self.T} {self.tdump} iso {self.P/bar} {self.P} {self.pdump}")
        elif md_type == "nvt":
            self.lmp.command(f"fix mynvt all nvt temp {self.T} {self.T} {self.tdump}")
        else:
            raise ValueError("md_type is not valid")

    def set_rigid_ads_lapmms(self, unfix_flag):
        self.lmp.command(f"group frame type 1 2 3 4")
        self.lmp.command(f"group co2 type 5 6")

        if self.Z_ads > 0 and unfix_flag:
            self.lmp.command(f"compute mdtemp co2 temp")
        self.lmp.command(f"compute_modify  mdtemp dynamic/dof yes")
        self.lmp.command(f"fix rigidco2 co2 rigid/small molecule")
        self.lmp.command(f"fix_modify rigidco2 dynamic/dof yes")
        #  self.lmp.command(f"fix rigidco2 co2 rigid/nvt/small molecule mol co2mol temp {self.T} {self.T} {100*timestep}")
        #  self.lmp.command(f"fix_modify rigidco2 dynamic/dof yes")
        self.lmp.command(f"fix myframe frame npt temp {self.T} {self.T} {self.tdump} iso {self.P/bar} {self.P/bar} {self.pdump}")
        #  self.lmp.command(f"fix myframe frame nvt temp {self.T} {self.T} {self.tdump}")

    def _insertion_acceptance(self, e_trial):
        exp_value = self.beta * (self.e - e_trial)
        if exp_value > 100:
            return True
        elif exp_value < -100:
            return False
        else:
            acc = min(1, self.V*self.beta*self.fugacity/self.Z_ads * np.exp(exp_value))
            return np.random.rand() < acc

    def _deletion_acceptance(self, e_trial):
        exp_value = -self.beta * (e_trial - self.e)
        if exp_value > 100:
            return True
        else:
            acc = min(1, (self.Z_ads+1)/self.V/self.beta/self.fugacity * np.exp(exp_value))
            return np.random.rand() < acc

    def _insertion(self):

        self.n_trial_insertion += 1
        self.Z_ads += 1
        self.Z_ads_in_loop += 1

        atoms_trial = self.atoms.copy() + self.atoms_ads.copy()
        pos = atoms_trial.get_positions()
        pos[-self.n_ads:] = random_position(pos[-self.n_ads:], atoms_trial.get_cell())
        atoms_trial.set_positions(pos)
        if vdw_overlap(atoms_trial, self.vdw, self.n_frame, self.n_ads, self.Z_ads-1):
            e_trial = 10**10 * kjmol
        else:
            atoms_trial.calc = self.calc_gcmc
            if self.flex_ads:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self._get_e_ads(atoms_trial) * electronvolt
            else:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self.e_ads * self.Z_ads
        if self._insertion_acceptance(e_trial):
            self.n_succ_insertion += 1
            self.n_tot_succ_steps += 1

            self.atoms.extend(atoms_trial[-self.n_ads:].copy())
            self.e = e_trial
        else:
            self.Z_ads -= 1
            self.Z_ads_in_loop -= 1

    def _deletion(self, i_ads, atoms_trial, picked_atoms_idx):

        self.n_trial_deletion += 1
        self.Z_ads -= 1
        self.Z_ads_in_loop -= 1

        del atoms_trial[picked_atoms_idx]
        atoms_trial.calc = self.calc_gcmc
        if self.flex_ads:
            e_trial = atoms_trial.get_potential_energy() * electronvolt - self._get_e_ads(atoms_trial) * electronvolt
        else:
            e_trial = atoms_trial.get_potential_energy() * electronvolt - self.e_ads * self.Z_ads
        if self._deletion_acceptance(e_trial):
            self.n_succ_deletion += 1
            self.e = e_trial

            del self.atoms[picked_atoms_idx]
            self.n_tot_succ_steps += 1
        else:
            self.Z_ads += 1
            self.Z_ads_in_loop += 1

    def _translation(self, i_ads, atoms_trial, pos, picked_atoms_idx):
        #  pos[picked_atoms_idx] += 0.5 * (np.random.rand(3) - 0.5)

        # recoded according to  tmmc org. code from this paper: 10.1021/jacs.4c15287
        pos[picked_atoms_idx] = random_position(pos[picked_atoms_idx], atoms_trial.get_cell())

        atoms_trial.set_positions(pos)
        if vdw_overlap(atoms_trial, self.vdw, self.n_frame, self.n_ads, i_ads):
            e_trial = 10**10 * kjmol
        else:
            atoms_trial.calc = self.calc_gcmc
            if self.flex_ads:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self._get_e_ads(atoms_trial) * electronvolt
            else:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self.e_ads * self.Z_ads
        acc = min(1, np.exp(-self.beta*(e_trial-self.e)))
        if acc > np.random.rand():
            self.atoms.set_positions(pos)
            self.e = e_trial
            self.n_tot_succ_steps += 1

    def _rotation(self, i_ads, atoms_trial, pos, picked_atoms_idx):
        #  pos[picked_atoms_idx] = _random_rotation(pos[picked_atoms_idx], circlefrac = 0.1)

        # recoded according to  tmmc org. code from this paper: 10.1021/jacs.4c15287
        pos[picked_atoms_idx] = _random_rotation(pos[picked_atoms_idx], circlefrac=0.25)
        pos[picked_atoms_idx] = pos[picked_atoms_idx] + 2.0 * (np.random.rand(3) - 0.5)


        atoms_trial.set_positions(pos)
        if vdw_overlap(atoms_trial, self.vdw, self.n_frame, self.n_ads, i_ads):
            e_trial = 10**10 * kjmol
        else:
            atoms_trial.calc = self.calc_gcmc
            if self.flex_ads:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self._get_e_ads(atoms_trial) * electronvolt
            else:
                e_trial = atoms_trial.get_potential_energy() * electronvolt - self.e_ads * self.Z_ads
        acc = min(1, np.exp(-self.beta*(e_trial-self.e)))
        if acc > np.random.rand():
            #  self.atoms = atoms_trial.copy()
            self.atoms.set_positions(pos)
            self.e = e_trial
            self.n_tot_succ_steps += 1

    def _load_model(self, model_path):
        from nequip.ase.nequip_calculator import NequIPCalculator
        return NequIPCalculator.from_deployed_model(
            model_path = model_path,
            device="cuda")

    def init_gcmc(self):

        # load gcmc model as a ase calc
        self.calc_gcmc = self._load_model(self.model_path_gcmc)

        # setting for gcmc
        if not self.flex_ads:
            #  self.atoms.calc = self.calc_gcmc

            atoms_ads = self.atoms_ads.copy()
            atoms_ads_e = self.atoms_ads.copy()
            atoms_ads_e.calc = self.calc_gcmc
            self.e_ads = atoms_ads_e.get_potential_energy() * electronvolt # eV to Hartree
        #  e = atoms.get_potential_energy() * electronvolt # eV to Hartree

        self.uptake = []
        self.adsorption_energy = []

        self.fl_status = open(f"{self.results_dir}/status.csv", "w")
        print("Steps,trial_insertion, succ_insertion,trial_deletaion,succ_deletion,mcswap(%),mcmoves(%)", file=self.fl_status)
        self.fl_status.flush()

        self.n_trial_insertion = 0
        self.n_succ_insertion = 0
        self.n_trial_deletion = 0
        self.n_succ_deletion = 0
        self.n_tot_succ_steps = 0
        self.old_n_tot_succ_steps = 0

    def run_gcmc(self, nmcswap, nmcmoves):

        ncycles = nmcswap + nmcmoves
        n_trial_mcswap = 0
        n_trial_mcmoves = 0


        # NOTE flex ads not tested
        self.atoms.calc = self.calc_gcmc
        if self.flex_ads:
            self.e = self.atoms.get_potential_energy() * electronvolt - self._get_e_ads(self.atoms) * electronvolt
        else:
            self.e = self.atoms.get_potential_energy() * electronvolt - self.e_ads * self.Z_ads # eV to Hartree

        # set V after MD steps
        self.cell = np.array(self.atoms.get_cell())
        self.V = np.linalg.det(self.cell) * angstrom**3

        self.Z_ads_in_loop = 0
        #  print("Start GCMC for %d steps" %ncycles)
        for interation in tqdm.trange(ncycles):
            # run GCMC exchanhes and MC moves
            if self.Z_ads > 0:
                i_ads = np.random.randint(self.Z_ads)
                picked_atoms_idx = slice(self.n_frame + self.n_ads*i_ads, self.n_frame + self.n_ads*(i_ads+1))

            atoms_trial = self.atoms.copy()
            pos = atoms_trial.get_positions()

            ixm = random.randint(0, ncycles) + 1
            if ixm <= nmcmoves and self.Z_ads != 0:
                #  xmcmove = random.uniform(0, 1)
                n_trial_mcmoves += 1
                if np.random.rand() < 0.5:
                    # Translation
                    self._translation(i_ads, atoms_trial, pos, picked_atoms_idx)
                else:
                    # Rotation
                    self._rotation(i_ads, atoms_trial, pos, picked_atoms_idx)
            else:
                #  xgcmc = random.uniform(0, 1)
                n_trial_mcswap += 1
                if np.random.rand() < 0.5 and self.Z_ads != 0:
                    # Deletion
                    self._deletion(i_ads, atoms_trial, picked_atoms_idx)
                else:
                    # insertion
                    self._insertion()

            self.uptake.append(self.Z_ads)
            self.adsorption_energy.append(self.e)

            # save uptake and energies at every step in
            try:
                np.save(f'{self.results_dir}/uptake_{self.P/bar}bar.npy', np.array(self.uptake))
                np.save(f'{self.results_dir}/adsorption_energy_{self.P/bar}bar.npy', np.array(self.adsorption_energy))
            except:
                print("saving npz failed but continue GCMC")

            # save status and geoms
            if self.n_tot_succ_steps % self.interval == 0 and self.old_n_tot_succ_steps != self.n_tot_succ_steps:
                try:
                    write(f'{self.results_dir}/trajectory_{self.P/bar}bar.extxyz',
                          self.atoms, append=True,
                      write_results=False, write_info=False)
                except:
                    print("saving extxyz failed but continue GCMC")

                print(f"{self.n_tot_succ_steps},"\
                      f"{self.n_trial_insertion},"\
                      f"{self.n_succ_insertion},"\
                      f"{self.n_trial_deletion},"\
                      f"{self.n_succ_deletion},"\
                      f"{100*n_trial_mcswap/ncycles:.1f},"\
                      f"{100*n_trial_mcmoves/ncycles:.1f}", file=self.fl_status)
                self.fl_status.flush()
                self.old_n_tot_succ_steps = self.n_tot_succ_steps


        # assign velocites to new added atoms
        if self.Z_ads_in_loop > 0:
            velocities = self.atoms.get_velocities()
            add_ads_atoms = self.atoms[-self.n_ads*self.Z_ads_in_loop:].copy()
            MaxwellBoltzmannDistribution(add_ads_atoms, temperature_K=self.T * 0.8)
            velocities[-self.n_ads*self.Z_ads_in_loop:] = add_ads_atoms.get_velocities()
            self.atoms.set_velocities(velocities)

    def run_gcmcmd(self, totalsteps, nmdsteps, nmcswap, nmcmoves):

        #  cell_params = self.atoms.get_cell_lengths_and_angles()

        unfix_flag = True
        for iteration in range(int(totalsteps/nmdsteps)):
            print("iter ", iteration+1)

            # for Zero Surface Tension (𝜎𝑎 = 0):
            #Change the volume and/or shape of the simulation box during a dynamics run
            #  self.lmp.command(f"fix no_deform all deform 1 x final 0 {cell_params[0]} y final 0 {cell_params[1]} z final 0 {cell_params[2]}")

            if not self.flex_ads:
                if self.Z_ads > 0 and unfix_flag:
                    self.lmp.command(f"unfix mynpt")
                    #  self.lmp.command(f"unfix mynvt")
                    self.set_rigid_ads_lapmms(unfix_flag)
                    unfix_flag = False
                elif self.Z_ads > 0 and not unfix_flag:
                    self.lmp.command(f"unfix rigidco2")
                    self.lmp.command(f"unfix myframe")
                    self.set_rigid_ads_lapmms(unfix_flag)

            self.lmp.command(f"run {nmdsteps}")
            #  self.lmp.command(f"write_dump all custom last_frame.lammpstrj id type x y z vx vy vz")
            # Retrieve the last frame as ASE Atoms
            self.lmp.command(f"write_data {self.results_dir}/data.last_frame")
            self.atoms = read_lammps_data(f"{self.results_dir}/data.last_frame", units=self.units_lmp)
            #  write("test.extxyz", self.atoms)
            #  self.atoms = self._lammps2ase_atoms()
            # GCMC
            print("GMCM steps...")
            self.run_gcmc(nmcswap, nmcmoves)
            #  self._atoms2lammpsdata(self.atoms)
            atoms_frame = self.atoms[:self.n_frame]
            atoms_ads = self.atoms[self.n_frame:]

            self._ase_to_lammpstrj(atoms_frame, atoms_ads, f"{self.results_dir}/last_frame.lammpstrj", timestep=0, units=self.units_lmp)

            self.lmp.command(f"undump 1")
            self.lmp.command(f"read_dump {self.results_dir}/last_frame.lammpstrj 0 x y z vx vy vz replace no purge yes add yes timestep no")

            # Set the ids to obtain the rigid body of each molecule individually
            mol_id = 0
            atom_id = 1
            for i in range(len(atoms_frame)):
                self. lmp.command(f"set atom {atom_id} mol {mol_id}")
                atom_id += 1
            mol_id +=1
            for i in range(len(atoms_ads)):
                self. lmp.command(f"set atom {atom_id} mol {mol_id}")
                if i != 0 and (i+1) % self.n_ads == 0:
                    mol_id +=1
                atom_id += 1

            self.lmp.command(f"dump 1 all atom {self.interval} {self.results_dir}/md_{self.T}K_{self.P/bar}bar.lammpstrj")
            self.lmp.command(f"dump_modify 1 append yes")
            #  self.lmp.command(f"rerun last_frame.lammpstrj dump x y z vx vy vz")

    def run_tmmcmd(self, nmdsteps):

        nmcmoves = max(20000, 5000 * self.Z_ads)
        #  nmcmoves = 100
        print("NVT-MC equilibration...")
        self.init_gcmc()
        self.run_gcmc(nmcswap=0, nmcmoves=nmcmoves)

        atoms_frame = self.atoms[:-self.Z_ads*3]
        atoms_ads = self.atoms[-self.Z_ads*3:]

        self._ase_to_lammpstrj(atoms_frame, atoms_ads, f"{self.results_dir}/last_frame.lammpstrj", timestep=0, units=self.units_lmp)

        self.lmp.command(f"unfix 1 all recenter")
        self.lmp.command(f"undump 1")
        # Due to velecity units covert issue from ase to lammps and do not need update it we ingnore reading velecities
        # and recreate velecities with next line
        self.lmp.command(f"read_dump {self.results_dir}/last_frame.lammpstrj 0 x y z replace no purge yes add yes timestep no")
        self.lmp.command(f"velocity all create {self.T} 4928459 ")

        # Set the ids to obtain the rigid body of each molecule individually
        mol_id = 0
        atom_id = 1
        for i in range(len(atoms_frame)):
            self. lmp.command(f"set atom {atom_id} mol {mol_id}")
            atom_id += 1
        mol_id +=1
        for i in range(len(atoms_ads)):
            self. lmp.command(f"set atom {atom_id} mol {mol_id}")
            if i != 0 and (i+1) % self.n_ads == 0:
                mol_id +=1
            atom_id += 1

        # set rigid for ads
        unfix_flag = True
        if not self.flex_ads and self.Z_ads > 0 and unfix_flag:
            self.lmp.command(f"unfix mynpt")
            self.set_rigid_ads_lapmms(unfix_flag)

        self.lmp.command(f"fix 1 all recenter INIT INIT NULL")
        self.lmp.command(f"dump 1 all atom {self.interval} {self.results_dir}/md_{self.T}K_{self.P/bar}bar.lammpstrj")
        self.lmp.command(f"dump_modify 1 append yes")
        self.lmp.command(f"run {nmdsteps}")
        #  self.lmp.command(f"write_dump all custom last_frame.lammpstrj id type x y z vx vy vz")
        # Retrieve the last frame as ASE Atoms
        self.lmp.command(f"write_data {self.results_dir}/data.last_frame")
        self.atoms = read_lammps_data(f"{self.results_dir}/data.last_frame", units=self.units_lmp)

        #  write("test.extxyz", self.atoms)
        #  self.atoms = self._lammps2ase_atoms()
        # GCMC
        #  print("GMCM steps...")
        #  self.run_gcmc(nmcswap, nmcmoves)
        #  #  self._atoms2lammpsdata(self.atoms)
        #  atoms_frame = self.atoms[:self.n_frame]
        #  atoms_ads = self.atoms[self.n_frame:]
        #
        #  self._ase_to_lammpstrj(atoms_frame, atoms_ads, f"{self.results_dir}/last_frame.lammpstrj", timestep=0, units=self.units_lmp)
        #
        #  self.lmp.command(f"undump 1")
        #  self.lmp.command(f"read_dump {self.results_dir}/last_frame.lammpstrj 0 x y z vx vy vz replace no purge yes add yes timestep no")
        #
        #  # Set the ids to obtain the rigid body of each molecule individually
        #  mol_id = 0
        #  atom_id = 1
        #  for i in range(len(atoms_frame)):
        #      self. lmp.command(f"set atom {atom_id} mol {mol_id}")
        #      atom_id += 1
        #  mol_id +=1
        #  for i in range(len(atoms_ads)):
        #      self. lmp.command(f"set atom {atom_id} mol {mol_id}")
        #      if i != 0 and (i+1) % self.n_ads == 0:
        #          mol_id +=1
        #      atom_id += 1
        #
        #  self.lmp.command(f"dump 1 all atom {self.interval} {self.results_dir}/md_{self.T}K_{self.P/bar}bar.lammpstrj")
        #  self.lmp.command(f"dump_modify 1 append yes")
        #  self.lmp.command(f"rerun last_frame.lammpstrj dump x y z vx vy vz")

