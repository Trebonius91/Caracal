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
use pbc_mod
use tblite_api_calculator
use tblite_api_context
use tblite_api_structure
use tblite_api_param
use tblite_api_result
use tblite_api_table
use tblite_api_error
use tblite_api_container
use mctc_io_convert
use, intrinsic :: iso_c_binding
implicit none
!     Loop indices
integer::i,j,k
!     The current coordinates (bohr)
real(kind=8),intent(in)::xyz_act(3,natoms)
!     The calculated gradient (hartree/bohr)
real(kind=8),intent(out)::g_act(3,natoms)
!     The calculated energy (hartree)
real(kind=8),intent(out)::e_act
!     The calculated partial charges (e)
real(kind=8) :: charges_act(natoms)
!     The calculated Wiberg bond orders
real(kind=8) :: wbo_act(natoms,natoms)
real(kind=8) ::wbo_thr,xsum
integer::ibmax
integer, allocatable :: imem(:)
!     The calculated dipole and quadrupole moments 
real(kind=8) :: dipole_act(3),quadrupole_act(6,1),dip
!     The calculated virial tensor
real(kind=8) :: virial_act(3,3)
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
!     If the input was read from a POSCAR file, apply the VASP symmetry directly!
!
if (coord_vasp) then
   c_periodic = .true.
   c_lattice(:,1) = vasp_a_vec/bohr
   c_lattice(:,2) = vasp_b_vec/bohr
   c_lattice(:,3) = vasp_c_vec/bohr
else
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
!     Obtain further properties of the system
!
call get_result_charges_api(verror,calc_results,charges_act)
call get_result_bond_orders_api(verror,calc_results,wbo_act)
call get_result_dipole_api(verror,calc_results,dipole_act)
call get_result_quadrupole_api(verror,calc_results,quadrupole_act)
call get_result_virial_api(verror,calc_results,virial_act)

write(84,*) "   FURTHER RESULTS: "
write(84,*)
write(84,*) "   * Partial charges"
write(84,*)
write(84,*) " # atom No.     element    Z(e)         q(e) "
write(84,'(a)') "------------------------------------------------------------------"
do i=1,natoms
   write(84,'(i7,a,a,i7,f18.7)') i,"             ",name(i),numbers(i),charges_act(i)
end do
write(84,'(a)') "------------------------------------------------------------------"
write(84,'(a,f18.7)') " Total charge of the system: ",sum(charges_act)
write(84,'(a)') "------------------------------------------------------------------"
write(84,*)
write(84,*) "   * Wiberg-Mayer bond orders"
wbo_thr=0.1d0
allocate(imem(natoms))
write(84,'(a)')
write(84,'("  largest (>",f4.2,") Wiberg bond orders for each atom")') wbo_thr
write(84,'(a)')
write(84,'(75("-"))')
write(84,'(5x,"#",3x,"Z",1x,"sym",2x,"total",t25,3(5x,"#",1x,"sym",2x,"WBO",2x))')
write(84,'(75("-"))')
do i=1,natoms
   do j=1,natoms
      imem(j)=j
   enddo
   call wibsort(natoms,i,imem,wbo_act)
   ibmax=0
   xsum =0.0d0
   do j=1,natoms
      if (wbo_act(i,j).gt.wbo_thr) ibmax=j
      xsum=xsum+wbo_act(i,j)
   enddo
   if (ibmax > 0) then
      write(84,'(i6,1x,i3,1x,a4,f6.3,1x,"--")',advance='no') &
         & i,numbers(i),name(i),xsum
   else
      write(84,'(i6,1x,i3,1x,a4,f6.3)') &
         & i,numbers(i),name(i),xsum
   end if
   do j = 1, ibmax, 3
      if (j > 1) then
         write(84,'(t25)', advance='no')
      end if
      do k = j, min(ibmax, j+2)
         write(84,'(i6,1x,a4,f6.3)',advance='no') &
            & imem(k),name(imem(k)),wbo_act(i,k)
      enddo 
      write(84,'(a)')
   end do
enddo
write(84,'(75("-"))')
write(84,'(a)')

dip=norm2(dipole_act)
write(84,*) "   * Dipole moment, from electron density (a.u.)"
write(84,*)
write(84,'(1x,"    x          y          z      ")')
write(84,'(a)') "------------------------------------------------------------------"
write(84,'(3f11.5,"  total (Debye): ",f11.5)') &
     & dipole_act(1),   dipole_act(2),   dipole_act(3), dip*autod
write(84,'(a)') "------------------------------------------------------------------"
write(84,*)
write(84,*) "   * Quadrupole moment, from electron density (a.u.)"
write(84,*)
write(84,'(a)',advance='no')'     xx         xy         yy         '
write(84,'(a)',advance='yes')'xz         yz         zz'
write(84,'(a)') "------------------------------------------------------------------"
write(84,'(6f11.5)')  quadrupole_act(1,1),quadrupole_act(2,1),quadrupole_act(3,1), &
            & quadrupole_act(4,1),quadrupole_act(5,1),quadrupole_act(6,1)
write(84,'(a)') "------------------------------------------------------------------"
write(84,*)
write(84,*) "   * Virial tensor components (a.u.)"
write(84,*)
write(84,'(a)') "------------------------------------------------------------------"
write(84,'(6f12.6)') virial_act(1,:)
write(84,'(6f12.6)') virial_act(2,:)
write(84,'(6f12.6)') virial_act(3,:)
write(84,'(a)') "------------------------------------------------------------------"
write(84,*)
write(84,*)
!
!     For NPT calculations: set stress tensor to obtained virial tensor
!
if (npt) then
   vir_ten=virial_act
end if
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
