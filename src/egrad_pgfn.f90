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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine egrad_pgfn: calculate the energy, gradient (and stress tensor)
!      from pGFN-FF with the GULP program, if included
!     geometry (and from this, a specific bonding configuration)
!
!     part of EVB
! 
subroutine egrad_pgfn(natoms,xyz_act,grad_act,energy_act) 
#ifdef GULP
use datatypes
use m_gulp_interface
#endif
use pbc_mod
implicit none 
integer,intent(in)::natoms
real(kind=8),intent(in)::xyz_act(3,natoms)
real(kind=8),intent(out)::grad_act(3,natoms)
real(kind=8),intent(out)::energy_act
real(kind=8)::xyz_local(3,natoms)
real(kind=8)::charges(natoms)
real(kind=8),allocatable::coord_cell(:,:)
real(kind=8)::strainderivatives(6)
real(kind=8)::bohr,evolt
character(len=120)::keywords,libraryfile
logical::lgradients
integer::ndim
integer::ierror

parameter (bohr=0.52917721092d0)
parameter (evolt=27.21138503d0)
!
!    Parameters to be transferred to the GULP library
!
!ierror               : integer i4   : If the value is 0 then the calculation was 
!                                       successful. 1 => GULP wasn't initialised!
!energy               : real    dp   : The energy of the system (in eV)
!charges(natoms)      : real    dp   : On return this will contain the charges (if 
!                                       a variable charge model is being used then 
!                                       these are the latest charges)
!lgradients           : logical      : If true then GULP will compute the first 
!                                       derivatives of the energy
!gradients(3:natoms)  : real    dp   : The Cartesian first derivatives of the energy 
!                                       with respect to atom coordinates (in eV/Angstroms)
!strainderivatives(6) : real    dp   : The first derivatives of the energy with respect
!                                       to Voight strains (in eV). The return values 
!                                       depend on ndim:
!                     :              : 0 => no strain derivatives
!                     :              : 1 => 1 strain derivative for xx strain
!                     :              : 2 => 3 strain derivatives for xx, yy, xy strains 
!                                         (in this order)
!                     :              : 3 => 6 strain derivatives for xx, yy, zz, yz, 
!                                       xz, xy strains (in this order)
!

!
!     Convert structure from bohr to Angstrom
!
xyz_local=xyz_act*bohr
lgradients=.true.

!
!     Set the periodicity
!     Currently, only rectangular boxes are supported by Caracal
!     TODO: add support for arbitrary box shapes!
!
!     Extra case: a POSCAR file contains the coordinates, use 
!     its shape for the box (in Angstroms!)
!
if (coord_vasp) then
   ndim=3
   allocate(coord_cell(3,ndim))
   coord_cell(:,1)=vasp_a_vec
   coord_cell(:,2)=vasp_b_vec
   coord_cell(:,3)=vasp_c_vec
else    
   if (periodic) then
      ndim=3
      allocate(coord_cell(3,ndim))
      coord_cell(1:3,1:3) = 0.0d0
      coord_cell(1,1)=boxlen_x
      coord_cell(2,2)=boxlen_y
      coord_cell(3,3)=boxlen_z
   else 
      ndim=0
      allocate(coord_cell(3,ndim))
   end if
end if

charges=0.d0
!
!     Set GULP input keywords for pGFN-FF calculation
!
keywords="gfnff gwolf"
libraryfile=" "
#ifdef GULP
call gulp_energy(coord_cell,xyz_local,charges,energy_act,grad_act, &
                 strainderivatives,lgradients,ierror)
#endif
if (ierror .ne. 0) then
   write(*,*) "Something went wrong during calling the gulp_energy subroutine!"
   write(*,*) "Please check your settings or the GULP installation!"
   call fatal
end if
!
!     comvert energy and derivatives from eV / eV/Ang  to Hartree / Hartree/bohr
!
grad_act=grad_act/bohr/evolt
energy_act=energy_act/evolt

return
end subroutine egrad_pgfn
