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
!     function erfinv: evalates the inverse of the error function 
!          in the interval [-1,1] 
!     adapted from the pseudocode for the Matlab function of the
!     same name; Matlab, version 4.2c, March 1995
!
!     part of EVB
!
function erfinv (x1)
use general
implicit none
real(kind=8)::erfinv,erf
real(kind=8)::x1,y1,z1
real(kind=8)::a1(4),b(4)
real(kind=8)::c(4),d(2)
real(kind=8)::sqrtpi
external::erf
!
!     coefficients for approximation to erfinv in central range
!
a1=(/ 0.886226899d0, -1.645349621d0, 0.914624893d0, -0.140543331d0/)
b=(/ -2.118377725d0,  1.442710462d0, -0.329097515d0, 0.012229801d0/)
!
!     coefficients for approximation to erfinv near endpoints
!
c=(/ -1.970840454d0, -1.624906493d0, 3.429567803d0,  1.641345311d0/)
d=(/  3.543889200d0,  1.637067800d0 /)
!
!     square root of pi
!
sqrtpi=1.772453850905516027d0
!
!     get an initial estimate for the inverse error function
!
if (abs(x1) .le. 0.7d0) then
   y1 = x1 * x1
   z1 = x1 * (((a1(4)*y1+a1(3))*y1+a1(2))*y1+a1(1))  &
       &     / ((((b(4)*y1+b(3))*y1+b(2))*y1+b(1))*y1+1.0d0)
else if (x1.gt.0.7d0 .and. x1.lt.1.0d0) then
   y1 = sqrt(-log((1.0d0-x1)/2.0d0))
   z1 = (((c(4)*y1+c(3))*y1+c(2))*y1+c(1)) / ((d(2)*y1+d(1))*y1+1.0d0)
else if (x1.lt.-0.7d0 .and. x1.gt.-1.0d0) then
   y1 = sqrt(-log((1.0d0+x1)/2.0d0))
   z1 = -(((c(4)*y1+c(3))*y1+c(2))*y1+c(1)) / ((d(2)*y1+d(1))*y1+1.0d0)
else
   write (iout,'(/," ERFINV  --  Illegal Argument to Inverse", &
     &      " Error Function")')
   call fatal
end if
!
!     use two steps of Newton-Raphson correction to increase accuracy
!
z1 = z1 - (erf(z1) - x1) / (2.0d0/sqrtpi * exp(-z1*z1))
z1 = z1 - (erf(z1) - x1) / (2.0d0/sqrtpi * exp(-z1*z1))
erfinv = z1

return
end function erfinv
