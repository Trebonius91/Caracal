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
!     function omega: calculates the inversion angle
!
!     part of QMDFF
!

real(kind=8) function omega (natoms,xyz,i,j,k,l)
use qmdff
implicit none
external::vecnorm
integer::ic,natoms,i,j,k,l
real(kind=8)::xyz(3,natoms),rd(3),re(3),rn(3)
real(kind=8)::rv(3),rnv,vecnorm,rkjn,rljn,rnn,rvn

do ic=1,3
   re(ic)=xyz(ic,i)-xyz(ic,j)
   rd(ic)=xyz(ic,k)-xyz(ic,j)
   rv(ic)=xyz(ic,l)-xyz(ic,i)
end do
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(re)
   call box_image(rd)
   call box_image(rv)
end if


call crossprod(re,rd,rn)
rnn=vecnorm(rn,3,1)
rvn=vecnorm(rv,3,1)

rnv=rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
omega=asin( rnv )

return
end function omega
