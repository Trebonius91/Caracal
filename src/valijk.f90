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
!     function valijk: take two cartesian vectors, normalize them
!       and calculate their dotproduct
!
!     part of QMDFF
!
real(kind=8) function valijk(natoms,xyz,j,k,i)
use qmdff
implicit none
external::vecnorm
integer::natoms,j,k,i,ic
real(kind=8)::ra(3),rb(3),rab,eps
real(kind=8)::xyz(3,natoms),vecnorm,ran,rbn

parameter (eps=1.d-14)

do ic=1,3
   ra(ic)=xyz(ic,j)-xyz(ic,i)
   rb(ic)=xyz(ic,k)-xyz(ic,i)
end do
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(ra)
   call box_image(rb)
end if


ran=vecnorm(ra,3,1)
rbn=vecnorm(rb,3,1)
rab=0.d0
do ic=1,3
   rab=rab+ra(ic)*rb(ic)
end do

if (abs(abs(rab)-1.d0).lt.eps) then
   rab=sign(1.d0,rab)
end if
valijk=acos(rab)

return
end function valijk
