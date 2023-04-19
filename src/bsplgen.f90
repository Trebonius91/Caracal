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
!     subroutine bsplgen: gets B-spline coefficients and derivatives for
!     a single PME atomic site along a particular direction
!
!     part of QMDFF 
!
subroutine bsplgen (w,thetai,bsorder)
!use qmdff
implicit none 
integer::i,j,k
integer::level
integer::bsorder 
real(kind=8)::w,denom
real(kind=8)::thetai(4,*)
real(kind=8),allocatable::bsbuild(:,:)

!
!     Allocate B-spline derivative coefficient storage
!
allocate(bsbuild(bsorder,bsorder))
!
!     set B-spline depth for partial charges or multipoles
!
level = 2
!
!     initialization to get to 2nd order recursion
!
bsbuild(2,2) = w
bsbuild(2,1) = 1.0d0 - w
!
!     perform one pass to get to 3rd order recursion
!
bsbuild(3,3) = 0.5d0 * w * bsbuild(2,2)
bsbuild(3,2) = 0.5d0 * ((1.0d0+w)*bsbuild(2,1) &
      &         +(2.0d0-w)*bsbuild(2,2))
bsbuild(3,1) = 0.5d0 * (1.0d0-w) * bsbuild(2,1)
!
!     compute standard B-spline recursion to desired order
!
do i = 4, bsorder
   k = i - 1
   denom = 1.0d0 / dble(k)
   bsbuild(i,i) = denom * w * bsbuild(k,k)
   do j = 1, i-2
      bsbuild(i,i-j) = denom * ((w+dble(j))*bsbuild(k,i-j-1) &
         &            +(dble(i-j)-w)*bsbuild(k,i-j))
   end do
   bsbuild(i,1) = denom * (1.0d0-w) * bsbuild(k,1)
end do
!
!     get coefficients for the B-spline first derivative
!
k = bsorder - 1
bsbuild(k,bsorder) = bsbuild(k,bsorder-1)
do i = bsorder-1, 2, -1
   bsbuild(k,i) = bsbuild(k,i-1) - bsbuild(k,i)
end do
bsbuild(k,1) = -bsbuild(k,1)
!
!     copy coefficients from temporary to permanent storage
!
do i = 1, bsorder
   do j = 1, level
      thetai(j,i) = bsbuild(bsorder-j+1,i)
   end do
end do
return


end subroutine bsplgen
