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
!     subroutine next_geo: The next geo from a concatenated xyz-file 
!     is read in and converted to bohr
! 
!     part of EVB
!
subroutine next_geo(coord,num_atoms,input_unit,has_next)
use evb_mod
use general
use pbc_mod

implicit none
integer::num_atoms,input_unit,atom_number,status,i
integer::readstat
real(kind=8),dimension(3,num_atoms)::coord
real(kind=8)::xr,yr,zr
real(kind=8)::x_tmp,y_tmp,z_tmp
character(len=80)::line,a80
character(len=3)::atom
real(kind=8),parameter::ang2bohr=0.52917720859d0
logical::has_next
has_next=.true.
!
!     If the main input file was a VASP file, assume that this 
!     file is also a VASP file (POSCAR or XDATCAR)
!
if (coord_vasp) then
   read(input_unit,*,iostat=status) a80
   if (status/=0) then
      has_next=.false.
      return
   end if
   do i=1,7
      read(input_unit,*,iostat=status)
      if (status .ne. 0) then
         write(*,*) "The structure file seems to be currupted!"
         call fatal
      end if
   end do
   if (vasp_selective) then
      read(input_unit,*,iostat=status)
      if (status .ne. 0) then
         write(*,*) "The structure file seems to be currupted!"
         call fatal
      end if
   end if
   do i=1,num_atoms
      read(input_unit,*,iostat=status) xr,yr,zr
      if (status .ne. 0) then
         write(*,*) "The structure file seems to be currupted!"
         call fatal
      end if

      if (vasp_direct) then
         x_tmp=xr
         y_tmp=yr
         z_tmp=zr
         xr=(x_tmp*vasp_a_vec(1)+y_tmp*vasp_b_vec(1)+z_tmp*vasp_c_vec(1))*vasp_scale
         yr=(x_tmp*vasp_a_vec(2)+y_tmp*vasp_b_vec(2)+z_tmp*vasp_c_vec(2))*vasp_scale
         zr=(x_tmp*vasp_a_vec(3)+y_tmp*vasp_b_vec(3)+z_tmp*vasp_c_vec(3))*vasp_scale
      end if
      coord(1,i)=xr
      coord(2,i)=yr
      coord(3,i)=zr
   end do
else

!
!     Read in a usual xyz file
! 
   read(input_unit,*,iostat=status) atom_number
   if (status/=0) then
      has_next=.false.
      return
   end if
   if (atom_number/=num_atoms) then
      write(*,*) atom_number,num_atoms
      write(*,*)"atom number mismatch! This should not happen..."
      stop 'atom number emergency!'
   end if
   read(input_unit,*,iostat=readstat) 
   if (readstat .ne. 0) then
      write(*,*) "The structure file seems to be currupted!"
      call fatal
   end if
!
!     if the PES TOPOL option is activated, give the element names     
! 
   if (pes_topol) then
      if (allocated(el_names)) deallocate(el_names)
      allocate(el_names(num_atoms))
   end if


   do i=1,num_atoms
      read(input_unit,*,iostat=readstat)atom,xr,yr,zr
      if (pes_topol) then
         el_names(i)=atom
      end if
      if (readstat .ne. 0) then
         write(*,*) "The structure file seems to be currupted!"
         call fatal
      end if

      coord(1,i)=xr
      coord(2,i)=yr
      coord(3,i)=zr
   end do
end if

return
end subroutine next_geo

