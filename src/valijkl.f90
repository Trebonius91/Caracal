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
!     function valijkl: calculates torsional angle
!
!     part of QMDFF
!
real(kind=8) function valijkl(natoms,xyz,i,j,k,l)
use qmdff
use pbc_mod
implicit none
external::vecnorm_loc,valijk
integer::ic,i,j,k,l,natoms
real(kind=8)::xyz(3,natoms),eps,ra(3),rb(3),rc(3)
real(kind=8)::na(3),nb(3),rab,rbc,thab,thbc,valijk
real(kind=8)::vecnorm_loc,nan,nbn,rcn,snanb,deter

parameter (eps=1.0d-14)
!
!     get torsion coordinate: bond vectors
!
do ic=1,3
   ra(ic)=xyz(ic,j)-xyz(ic,i)
   rb(ic)=xyz(ic,k)-xyz(ic,j)
   rc(ic)=xyz(ic,l)-xyz(ic,k)
end do
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(ra)
   call box_image(rb)
   call box_image(rc)
end if

!
!     determinante of rb,ra,rc (the thee bond vectors)
!
deter= ra(1)*(rb(2)*rc(3)-rb(3)*rc(2)) &
  &   -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1)) &
  &   +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

thab=valijk(natoms,xyz,i,k,j)
thbc=valijk(natoms,xyz,j,l,k)
!
!     calculate cross product
!
call crossprod(ra,rb,na)
call crossprod(rb,rc,nb)
nan=vecnorm_loc(na,3,1)
nbn=vecnorm_loc(nb,3,1)

snanb=0.0d0
do ic=1,3
   snanb=snanb+na(ic)*nb(ic)
end do
if (abs(abs(snanb)-1.d0).lt.eps) then
   snanb=sign(1.d0,snanb)
end if

valijkl=acos(snanb)

!
!     the gradient dphir is only compatible with this subroutine
!     if the statement below is commented out. If not, opt. and
!     Hessian show large errors and imags. I don't understand
!     this entirely but thats how it is. 
!     SG, Sat May 24 11:41:42 CEST 2014
!
!     if (deter.lt.0) then 
!        valijkl=2.d0*pi-valijkl
!     end if
!

return
end function valijkl
