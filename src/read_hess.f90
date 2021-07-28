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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine read_hess: Read in cartesian hessian from ref.input 
!     file and convert to internals
!
!     part of EVB
!

subroutine read_hess(atom_num,hess_int,geo_int,geo_xyz1,grad_int)
use evb_mod
implicit none
character(len=80) a80
integer::atom_num,k,i,maxcycle,rest,nat,j
real(kind=8),dimension(atom_num*3-6)::geo_int,grad_int
real(kind=8),dimension(atom_num*3)::grad_xyz
real(kind=8),dimension(3,atom_num)::geo_xyz1
real(kind=8),dimension(atom_num*atom_num*9)::hess_xyz
real(kind=8),dimension(atom_num*3,atom_num*3)::hess_xyz1
real(kind=8),dimension(nat6,nat6)::hess_int
real(kind=8),dimension(atom_num*3)::freqs  ! TEST for hessian correction
logical::redun
redun=.false.
nat=atom_num
do
   read(ref_input_unit,'(a)')a80
   if(index(a80,'HESSIAN ').ne.0) then
      maxcycle=(nat*nat*9)/5
!      write(*,*) "maxcyle",maxcycle
      do i=1,maxcycle
         k=(i-1)*5
         read(ref_input_unit,*) hess_xyz(k+1),hess_xyz(k+2),hess_xyz(k+3),hess_xyz(k+4),hess_xyz(k+5)
      end do
      rest=(nat*nat*9)-maxcycle*5
      k=maxcycle*5
      if (rest .eq. 4) then
         read(ref_input_unit,*) hess_xyz(k+1),hess_xyz(k+2),hess_xyz(k+3),hess_xyz(k+4)
      else if (rest .eq. 3) then
         read(ref_input_unit,*) hess_xyz(k+1),hess_xyz(k+2),hess_xyz(k+3)
      else if (rest .eq. 2) then
         read(ref_input_unit,*) hess_xyz(k+1),hess_xyz(k+2)
      else if (rest .eq. 1) then
         read(ref_input_unit,*) hess_xyz(k+1)
      end if
      goto 50
   end if
end do
50 continue
!
!   For usage of the RP-EVB coupling term: correct the cartesian hessians for
!   imaginary frequencies; remove them
!
!write(*,*) rp_evb
if (rp_evb) then
!   call project_hess(hess_xyz)
end if
!
!   call the transformation routine to transform the hessian to internal coordinates
!
call hess2int(geo_xyz1,geo_int,hess_int,hess_xyz,grad_int) 
!
!   For usage of the RP-EVB coupling term: correct the internal hessians for
!   imaginary frequencies; remove them
!
if (rp_evb) then
   call project_hess(hess_int)
end if

return
end subroutine read_hess
