#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
#   CARACAL - Ring polymer molecular dynamics and rate constant calculations
#             on black-box generated potential energy surfaces
#
#   Copyright (c) 2023 by Julien Steffen (mail@j-steffen.org)
#                         Stefan Grimme (grimme@thch.uni-bonn.de) (QMDFF code)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#
#     Python-file call_mace: Initialization and usage of MACE, called via the 
#      C-wrapper routines
#

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

#
#    Predefine global objects that are preserved in memory during simulation
#

atoms=None
calc=None

#
#    The initialization routine: read in the MACE MLIP and define it for the 
#      current geometry/system
#
def init_mace(mlip_file,coord_file,set_disp):
#
#    Redefine global variables
#
   global atoms, calc
#
#    Define the MACE model
#
#
#    Returns a model with D3 dispersion correction
#
   if set_disp:
      macemp = mace_mp(model=mlip_file,dispersion=True)
#
#    Returns a model without D3 dispersion correction
#
   else:
      macemp = mace_mp(model=mlip_file,dispersion=False) 

#
#    Read in the initial geometry from the POSCAR file and initialize the atoms object
#
   atoms = read(coord_file)

   natoms = len(atoms)
#
#    Setup the MACE foundation model
#
   calc = macemp

   if atoms is None or calc is None:
      raise RuntimeError("init_mace must be called before ase_mace")

#
#    The energy+gradient calculation routine: call MACE to calculate energy 
#     and gradient for the current structure
#
def ase_mace(coords,unitcell,natoms):
#  
#    Redefine global variables
#

   global atoms, calc


   if atoms is None or calc is None:
      raise RuntimeError("init_mace must be called before ase_mace")

#
#    Update coordinates in global atoms object
#
   atoms.set_cell(unitcell)
   atoms.set_positions(coords)
#
#    Perform the actual MACE calculation for potential energy 
#     and the forces
#
   atoms.calc=calc
   energy=atoms.get_potential_energy()
   gradient=atoms.get_forces()
#
#    Return the gradient as usual list insteaf of a numpy array to avoid a 
#     Segmentation Fault on the C side!
#
   return energy,gradient.tolist()
