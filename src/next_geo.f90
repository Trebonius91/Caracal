!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   EVB-QMDFF - RPMD molecular dynamics and rate constant calculations on
!               black-box generated potential energy surfaces
!
!   Copyright (c) 2021 by Julien Steffen (steffen@pctc.uni-kiel.de)
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

implicit none
integer::num_atoms,input_unit,atom_number,status,i
integer::readstat
real(kind=8),dimension(3,num_atoms)::coord
real(kind=8)::xr,yr,zr
character(len=80)::line
character(len=3)::atom
real(kind=8),parameter::ang2bohr=0.52917720859d0
logical::has_next
has_next=.true.
!
!      special case for simple Mueller-Brown surface: one 2D particle
!      z-coordinate always zero!
!
if (mueller_brown) then
   read(input_unit,*,iostat=status) xr,yr
   if (status/=0) then
      has_next=.false.
      return
   end if
   coord(1,1)=xr
   coord(2,1)=yr
   coord(3,1)=0.d0
   return
end if
!
!      usual molecular structures
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

do i=1,num_atoms
   read(input_unit,*,iostat=readstat)atom,xr,yr,zr
   if (readstat .ne. 0) then
      write(*,*) "The structure file seems to be currupted!"
      call fatal
   end if

   coord(1,i)=xr
   coord(2,i)=yr
   coord(3,i)=zr
end do

return
end subroutine next_geo

