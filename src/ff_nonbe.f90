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
!     subroutine ff_nonbe: calculate again energies for nonbonded parts of FF
!           (without gradient)
!
!     part of QMDFF
!
subroutine ff_nonbe(n,at,xyz,q,r0ab,zab,r0094,sr42,c66ab,epair,e)
use qmdff
implicit none
integer::n,at(*)
real(kind=8)::xyz(3,n),epair(*),q(n),e(3)
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),c66ab(n,n)
integer::i1,i2,iz1,iz2,k
real(kind=8)::r2,r,r4,r6,r06,R0,t6,t8,c6t6,c6t8,t27,drij
real(kind=8)::dx,dy,dz,c6,x,alpha,oner,aa,damp,damp2,e0,scl

e=0
!
!     if no term exists, return
!
if(nnci.le.1)return

!
!     loop over all nonbonded terms
!
do k=1,nnci
   i1=nci(1,k)
   i2=nci(2,k)
   dx=xyz(1,i1)-xyz(1,i2)
   dy=xyz(2,i1)-xyz(2,i2)
   dz=xyz(3,i1)-xyz(3,i2)
   r2=dx*dx+dy*dy+dz*dz
   r =sqrt(r2)
   oner =1.0d0/r
   iz1=at(i1)
   iz2=at(i2)
   R0=r0094(iz1,iz2)
   c6=c66ab(i2,i1)
   r4=r2*r2
   r6=r4*r2
   r06=R0**6
   t6=r6+r06
   t8=r6*r2+r06*R0*R0
   c6t6=c6/t6
   c6t8=c6/t8
   t27=sr42(iz1,iz2)*c6t8

   e0=-c6t6-t27
   epair(k)=e0*eps2(nci(3,k))
   e(1)=e(1)+e0*eps2(nci(3,k))

   e0=q(i1)*q(i2)*oner
   epair(k)=epair(k)+e0*eps1(nci(3,k))
   e(2)=e(2)+e0*eps1(nci(3,k))

   x    =zab(iz1,iz2)
   alpha=r0ab(iz1,iz2)
   t27  =x*dexp(-alpha*r)

   e0   =t27*oner
   epair(k)=epair(k)+e0*eps2(nci(3,k))
   e(3)=e(3)+e0*eps2(nci(3,k))
end do

return
end subroutine ff_nonbe
