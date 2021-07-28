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
!     subroutine orca_grad: calculate orca energy and gradient 
!         for ab-initio MD
!
!     part of EVB
!

subroutine orca_grad(xyz2,pot_grad,e_evb)
use evb_mod 
use general
implicit none
real(kind=8)::xyz2(3,natoms)
real(kind=8)::pot_grad(3,natoms,1)
real(kind=8)::e_evb
integer::i,j


!
!     write orca input file
!
open (unit=19,file="egrad.inp",status="replace")
write (19,*) orca_com
write (19,*) "*xyz 0 1"
do i=1,natoms
   write(19,*) name(i),xyz2(:,i)*bohr
end do
write(19,*) "*"
close(19)

call system("orca egrad.inp > egrad.out")

!
!    read in orca output file
!
open(unit=29,file="egrad.engrad",status="old")
do i=1,7
   read(29,*)
end do
read(29,*) e_evb
do i=1,3
   read(29,*)
end do
do i=1,natoms
   do j=1,3
      read(29,*) pot_grad(j,i,1)
   end do
end do
close(29)
return

end subroutine orca_grad

