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
!     subroutine mace_init: initialize the MACE MLIP via the ASE with a 
!      Python/C Wrapper mechanism
!
!     part of EVB
! 
subroutine mace_init(natoms,xyz_init,names_init) 
use pbc_mod
use inter_mace
implicit none 
integer,intent(in)::natoms
real(kind=8),intent(in)::xyz_init(3,natoms)
character(len=2),intent(in)::names_init(natoms)
real(kind=8)::charges(natoms)
real(kind=8),allocatable::coord_cell(:,:)
character(len=4)::atomicsymbols(natoms)
character(len=80)::keywords,libraryfile
logical::lgulpoutput
integer::ndim
real(kind=8),allocatable::cell(:,:)
real(kind=8)::xyz(3,natoms)

!
!     If a VASP POSCAR file is the input, take its box shape!
!

allocate(coord_cell(3,3))
coord_cell(1:3,1:3) = 0.0d0
if (coord_vasp) then
   ndim=3
   coord_cell(:,1)=vasp_a_vec
   coord_cell(:,2)=vasp_b_vec
   coord_cell(:,3)=vasp_c_vec
else
   if (periodic) then
      ndim=3
      coord_cell(1,1)=boxlen_x
      coord_cell(2,2)=boxlen_y
      coord_cell(3,3)=boxlen_z
   else
      ndim=0
   end if
end if
!
!  Set atomic coordinates
!
xyz(:,:)=xyz_init(:,:)
atomicsymbols(:)=names_init(:)

write(*,*) "Init routine!"
call init_mace(xyz,atomicsymbols,natoms)

return
end subroutine
