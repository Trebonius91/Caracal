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
!     subroutine egrestrain: energy/gradient of harmonic restrain 
!          potential for hessian fit
!
!     part of QMDFF
!
subroutine egrestrain(n,xyz,g,ers)
use qmdff
implicit none
real(kind=8)::xyz(3,n),ers
real(kind=8)::g(3,n)
integer::n
real(kind=8)::ra(3),rij,kij,fac,r
integer::m,ic,i,j

ers=0
do m=1,nrs
   i=rest(1,m)
   j=rest(2,m)
!
!     distance between two atoms..
!
   r=sqrt((xyz(1,i)-xyz(1,j))**2 &
     &   +(xyz(2,i)-xyz(2,j))**2 &
     &   +(xyz(3,i)-xyz(3,j))**2)
   do ic=1,3
      ra(ic)=xyz(ic,i)-xyz(ic,j)
   end do
   rij=vrest(1,m)
   kij=vrest(2,m)
!
!     energy formula
!
   ers=ers+kij*(r-rij)*(r-rij)
   fac=2.0d0*kij*(r-rij)/r
!
!     gradient entries
!
   do ic=1,3
      g(ic,i)=g(ic,i)+fac*ra(ic)
      g(ic,j)=g(ic,j)-fac*ra(ic)
   end do
end do

return
end subroutine egrestrain
