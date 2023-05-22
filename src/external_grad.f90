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
!     subroutine external_grad: calculate energy and gradient with a chosen
!         arbitrary external program for the current geometry!
!
!     part of EVB
!

subroutine external_grad(xyz2,pot_grad,e_evb)
use evb_mod 
use general
implicit none
real(kind=8)::xyz2(3,natoms)
real(kind=8)::pot_grad(3,natoms,1)
real(kind=8)::e_evb
integer::i,j


!
!     Write current structure that shall be calculated with the external program
!
open (unit=19,file="coords.xyz",status="replace")
write(19,*) natoms
write(19,*)
do i=1,natoms
   write(19,*) name(i),xyz2(:,i)*bohr
end do
close(19)

call system(trim(symlink_ext))
!
!    read in energy and gradient of the external calculation
!
open(unit=29,file="egrad_out.dat",status="old")
read(29,*) e_evb
do i=1,natoms
   read(29,*) pot_grad(:,i,1)
end do
close(29)
return

end subroutine external_grad

