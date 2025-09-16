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
#     Python-file call_uma: Initialization and usage of UMA, called via the 
#      C-wrapper routines
#

from ase import units
from ase import build
from ase.io import read, write
from ase.io.trajectory import Trajectory
import urllib3
import numpy as np
from fairchem.core import pretrained_mlip, FAIRChemCalculator
import sys

#
#    Predefine global objects that are preserved in memory during simulation
#

atoms=None
calc=None

#
#    The initialization routine: read in the UMA MLIP and define it for the 
#      current geometry/system
#
def init_uma(mlip_file,coord_file,task,calc_device):
#
#    Redefine global variables
#
   global atoms, calc
#
#    Define the UMA model
#    The UMA model file shall be given without file ending!
#
#    Decide if calculations shall be done on CPU or GPU
#
   predictor = pretrained_mlip.get_predict_unit(mlip_file, device=calc_device)

   calc = FAIRChemCalculator(predictor, task_name=task)
#
#    Read in the initial geometry from the POSCAR file and initialize the atoms object
#
   atoms = read(coord_file)

   natoms = len(atoms)
#
#    Setup the MACE foundation model
#
   if atoms is None or calc is None:
      raise RuntimeError("init_uma must be called before ase_uma")

#
#    The energy+gradient calculation routine: call UMA to calculate energy 
#     and gradient for the current structure
#
def ase_uma(coords,unitcell,natoms):
#  
#    Redefine global variables
#

   global atoms, calc


   if atoms is None or calc is None:
      raise RuntimeError("init_uma must be called before ase_uma")
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
