!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CARACAL - Ring polymer molecular dynamics and rate constant calculations
!             on black-box generated potential energy surfaces
!
!   Copyright (c) 2023 by Julien Steffen (mail@j-steffen.org)
!                         Stefan Grimme (grimme@thch.uni-bonn.de) (QMDFF code)
!
!   Permission is hereby granted, free of charge, to any person obtaining a
!   copy of this software and associated documentation files (the "Software"),
!   to deal in the Software without restriction, including without limitation
!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!   and/or sell copies of the Software, and to permit persons to whom the
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be included in
!   all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!   DEALINGS IN THE SOFTWARE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     Subroutine egrad_tblite: Calculates energy and gradient of a (nonperiodic)
!      structure via either GFN1-xTB or GFN2-xTB over the tblite library 
!      included into Caracal
!
!     part of tblite
!

subroutine egrad_tblite(xyz_act,g_act,e_act)
use evb_mod
use general
use qmdff
use tblite_api_calculator
use tblite_api_context
use tblite_api_structure
use tblite_api_param
use tblite_api_result
use tblite_api_table
use tblite_api_error
use tblite_api_container
use, intrinsic :: iso_c_binding
implicit none
!     Loop indices
integer::i,j
!     The current coordinates (bohr)
real(kind=8),intent(in)::xyz_act(3,natoms)
!     The calculated gradient (hartree/bohr)
real(kind=8),intent(out)::g_act(3,natoms)
!     The calculated energy (hartree)
real(kind=8),intent(out)::e_act
!     Data for the communication with xtblite
!    The calculaton object (calculator)
type(c_ptr) :: calc_xtb
!    The context of the calculation
type(c_ptr) :: calc_context
!    All properties of the molecule to be calculated
type(c_ptr) :: calc_structure
!    The management of parametrization records (?)
type(c_ptr) :: calc_param
!    The results container
type(c_ptr) :: calc_results
!    The ALPB solvation container 
type(c_ptr) :: calc_alpb
!    The CPCM solvation container
type(c_ptr) :: calc_cpcm

!> Calculation context
!type(context_type) :: ctx
!> Molecular structure data
!type(structure_type) :: mol
!  Attributes for molecule class
type(c_ptr) :: verror
type(vp_error),pointer :: error
integer(c_int) :: numbers(natoms)
real(c_double) :: positions(3,natoms)
real(c_double) :: c_charge
integer(c_int) :: c_uhf
real(c_double) :: c_lattice(3,3)
logical(c_bool) :: c_periodic(3)
type(vp_result), pointer :: res


integer::ntest,natoms_test
real(kind=8)::energy_test

debug= .true.

!
!     Initialize calculation context instance
!
calc_context = new_context_api()

!
!     Translate element numbers and atomic positions 
!
numbers(:) = elem_index(:)

positions(:,:) = xyz_act(:,:)

!
!     Define total charge and if a unrestricted calculation shall be done
!
c_charge = xtb_charge
c_uhf = 0
!
!     For periodic systems (currently only cubic): initialize it
!
if (periodic) then
   c_periodic = .true. 
   c_lattice = 0.d0
   c_lattice(1,1)= boxlen_x
   c_lattice(2,2) = boxlen_y
   c_lattice(3,3) = boxlen_z
else 
   c_periodic = .false.
   c_lattice = 0.d0
   c_lattice(1,1)= 100.0
   c_lattice(2,2) = 100.0
   c_lattice(3,3) = 100.0
end if

!
!     Generate the structure instance: full definition of the system
!
calc_structure = new_structure_api(verror,natoms,numbers,positions,c_charge, &
          & c_uhf,c_lattice,c_periodic)


!
!     Depending on the chosen GFN-xTB Hamiltonian, initialize the calculator 
!     instance
!
if (trim(adjustl(hamil_string)) .eq. "GFN1-XTB") then
   calc_xtb = new_gfn1_calculator_api(calc_context,calc_structure)
else if (trim(adjustl(hamil_string)) .eq. "GFN2-XTB") then
   calc_xtb = new_gfn2_calculator_api(calc_context,calc_structure)
else if (trim(adjustl(hamil_string)) .eq. "IPEA1-XTB") then
   calc_xtb = new_ipea1_calculator_api(calc_context,calc_structure)
end if
!
!     Set further options for the calculator (details)
!
!     The SCF convergence criterion (electronic: given value, 
!     density(?): given value*20)
!
call set_calculator_accuracy_api(calc_context,calc_xtb,xtb_accuracy)
!
!     The maximum number of SCF iterations
!
call set_calculator_max_iter_api(calc_context,calc_xtb,xtb_maxiter)
!
!     The electronic temperature for the Fermi occupation smearing (Kelvin)
!
call set_calculator_temperature_api(calc_context,calc_xtb,xtb_el_temp)

!
!     If the ALPB solvation model was called, perform these two additional calls 
!     in order to modify the energy and gradient of the actual calculation
!     If a solvent species is given by a string, it will be preferred, else, a 
!     dielectric constant epsilon will be read
! 
if (trim(adjustl(solv_string)) .eq. "ALPB") then
   if (exist_spec) then
      calc_alpb = new_alpb_solvation_api(calc_context,calc_structure, &
                    & calc_xtb,trim(solv_spec))
   else 
      calc_alpb = new_alpb_solvation_api(calc_context,calc_structure, &
                    & calc_xtb,trim(solv_epsilon))
   end if
   call push_back_api(calc_context,calc_xtb,calc_alpb)
end if
!
!     If the CPCM solvation model was called, perform these two addditional calls 
!     in order to modify the energy and gradient of the actual calculation
!     If a solvent species is given by a string, it will be preferred, else, a 
!     dielectric constant epsilon will be read
! 

if (trim(adjustl(solv_string)) .eq. "CPCM") then
   if (exist_spec) then
      calc_cpcm = new_cpcm_solvation_api(calc_context,calc_structure, &
                    & calc_xtb,trim(solv_spec))
   else 
      calc_cpcm = new_cpcm_solvation_api(calc_context,calc_structure, &
                    & calc_xtb,trim(solv_epsilon))
   end if
   call push_back_api(calc_context,calc_xtb,calc_cpcm)
end if

calc_results = new_result_api()

write(84,'(a,i12,a)') " ------ Energy+gradient calculation No. ",xtb_calc_num," ------- "
write(84,*) " "
xtb_calc_num = xtb_calc_num + 1
!
!     Do the actual calculation
!
call get_singlepoint_api(calc_context,calc_structure,calc_xtb,calc_results)

!
!     Obtain the resulting energy and gradient
!
call get_result_energy_api(verror,calc_results,e_act)

call get_result_gradient_api(verror,calc_results,g_act)

!
!     Delete the relevant objects for this cycle
!

call delete_calculator_api(calc_xtb)
call delete_structure_api(calc_structure)
if (trim(adjustl(solv_string)) .eq. "ALPB") then
   call delete_container_api(calc_alpb)
else if (trim(adjustl(solv_string)) .eq. "CPCM") then
   call delete_container_api(calc_cpcm)
end if
call delete_result_api(calc_results)
call delete_context_api(calc_context)

return
end subroutine egrad_tblite
