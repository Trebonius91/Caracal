#This file contains the wrapper routines needed for a direct 
#communication between Fortran and Python for communication with ASE

#include <Python.h>
#include <stdio.h>

void init_mace() {
   printf("test!")
}

void ase_mace(double *coords, double energy, double *gradient, int natoms) {
    


}	
