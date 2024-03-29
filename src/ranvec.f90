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
!     subroutine ranvec: generates a unit vector in 3-dimensional
!         space with uniformly distributed random orientation
!
!     literature references:
!
!     G. Marsaglia, Ann. Math. Stat., 43, 645 (1972)
!
!     R. C. Rapaport, The Art of Molecular Dynamics Simulation,
!     2nd Edition, Cambridge University Press, 2004, Section 18.4
!
!     part of EVB
!
subroutine ranvec (vector)
use general
implicit none
real(kind=8)::x1,y1,s
real(kind=8)::random
real(kind=8)::vector(3)
!
!     get a pair of appropriate components in the plane
!
s = 2.0d0
do while (s .ge. 1.0d0)
   x1 = 2.0d0 * random () - 1.0d0
   y1 = 2.0d0 * random () - 1.0d0
   s = x1**2 + y1**2
end do
!
!     construct the 3-dimensional random unit vector
!
vector(3) = 1.0d0 - 2.0d0*s
s = 2.0d0 * sqrt(1.0d0 - s)
vector(2) = s * y1
vector(1) = s * x1

return
end subroutine ranvec
