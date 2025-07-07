# This script serves to call the MACE energy and gradients for a Caracal
#  calculation, where Caracal uses the respective MACE commands
from ase import units
from ase import build
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.md.npt import NPT
import numpy as np
import time
from timeit import default_timer as timer
import datetime
from datetime import timedelta
import os.path
import os
from mace.calculators import mace_mp
import sys

def init_mace:
#
#   Define the MACE model
#
   #macemp = mace_mp() # return the default medium ASE calculator equivalent to mace_mp(model="medium")
   #macemp = mace_mp(model="large") # return a larger model
   #macemp = mace_mp(model="https://tinyurl.com/y7uhwpje") # downlaod the model at the given url
   macemp = mace_mp(model="./mace_fine_tuning_run2_run-1.model") # return a model with D3 dispersion correction

#
#    Read in the initial geometry from the POSCAR file and initialize the atoms object
#
   atoms = read('POSCAR')

   natoms = len(atoms)
#
#    Setup the MACE foundation model
#
   calc = macemp

def ase_mace(coords,unitcell,natoms):


   xyz = np.zeros((natoms,3))
   xyz = coords
   atoms.set_cell(unitcell)
   atoms.set_positions(xyz)
#
#    Perform the actual MACE calculation for potential energy 
#     and the forces
#
   atoms.calc=calc
   energy=atoms.get_potential_energy()
   gradient=atoms.get_forces()

   return energy,gradient
