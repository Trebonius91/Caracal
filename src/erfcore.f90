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
!     subroutine erfcore: approximation of the error function from 
!         radional function
!     literature reference:
!
!     W. J. Cody, "Rational Chebyshev Approximations for the Error
!     Function", Mathematics of Computation, 631-638, 1969
!
!     adapted from an original program written by W. J. Cody,
!     Mathematics and Computer Science Division, Argonne National
!     Laboratory, Argonne, IL 60439
!
!     part of EVB
!
subroutine erfcore (arg,result,mode)
implicit none
integer::i,mode
real(kind=8)::arg,result
real(kind=8)::x,y,ysq
real(kind=8)::del,sqrpi
real(kind=8)::xnum,xden
real(kind=8)::xtiny,xbig
real(kind=8)::a(5),b(4)
real(kind=8)::c(9),d(8)
real(kind=8)::p(6),q(5)
!
!     mathematical and machine-dependent constants
!
sqrpi = 5.6418958354775628695d-1 
xtiny = 1.11d-16 
xbig = 26.543d0 
!
!     coefficients for approximation to erf in first interval
!
a =(/ 3.16112374387056560d0,  1.13864154151050156d2, &
     &          3.77485237685302021d2,  3.20937758913846947d3, &
     &          1.85777706184603153d-1 /)
b =(/ 2.36012909523441209d1,  2.44024637934444173d2, &
     &          1.28261652607737228d3,  2.84423683343917062d3 /)
!
!     coefficients for approximation to erfc in second interval
!
c =(/ 5.64188496988670089d-1, 8.88314979438837594d0, &
     &          6.61191906371416295d1,  2.98635138197400131d2, &
     &          8.81952221241769090d2,  1.71204761263407058d3, &
     &          2.05107837782607147d3,  1.23033935479799725d3, &
     &          2.15311535474403846d-8 /)
d =(/ 1.57449261107098347d1,  1.17693950891312499d2, &
     &          5.37181101862009858d2,  1.62138957456669019d3, &
     &          3.29079923573345963d3,  4.36261909014324716d3, &
     &          3.43936767414372164d3,  1.23033935480374942d3 /)
!
!     coefficients for approximation to erfc in third interval
!
p =(/ 3.05326634961232344d-1, 3.60344899949804439d-1, &
     &          1.25781726111229246d-1, 1.60837851487422766d-2, &
     &          6.58749161529837803d-4, 1.63153871373020978d-2 /)
q =(/ 2.56852019228982242d0,  1.87295284992346047d0, &
     &          5.27905102951428412d-1, 6.05183413124413191d-2, &
     &          2.33520497626869185d-3 /)
!
!     store the argument and its absolute value
!
x = arg
y = abs(x)
!
!     evaluate error function for |x| less than 0.46875
!
if (y .le. 0.46875d0) then
   ysq = 0.0d0
   if (y .gt. xtiny)  ysq = y * y
   xnum = a(5) * ysq
   xden = ysq
   do i = 1, 3
      xnum = (xnum + a(i)) * ysq
      xden = (xden + b(i)) * ysq
   end do
   result = x * (xnum + a(4)) / (xden + b(4))
   if (mode .ne. 0)  result = 1.0d0 - result
!
!     get complementary error function for 0.46875 <= |x| <= 4.0
!
else if (y .le. 4.0d0) then
   xnum = c(9) * y
   xden = y
   do i = 1, 7
      xnum = (xnum + c(i)) * y
      xden = (xden + d(i)) * y
   end do
   result = (xnum + c(8)) / (xden + d(8))
   ysq = aint(16.0d0*y) / 16.0d0
   del = (y-ysq) * (y+ysq)
   result = exp(-ysq*ysq-del) * result
   if (mode .eq. 0) then
      result = 1.0d0 - result
      if (x .lt. 0.0d0)  result = -result
   else
      if (x .lt. 0.0d0)  result = 2.0d0 - result
   end if
!
!     get complementary error function for |x| greater than 4.0
!
else
   result = 0.0d0
   if (y .lt. xbig) then
      ysq = 1.0d0 / (y * y)
      xnum = p(6) * ysq
      xden = ysq
      do i = 1, 4
         xnum = (xnum + p(i)) * ysq
         xden = (xden + q(i)) * ysq
      end do
      result = ysq * (xnum + p(5)) / (xden + q(5))
      result = (sqrpi - result) / y
      ysq = aint(16.0d0*y) / 16.0d0
      del = (y-ysq) * (y+ysq)
      result = exp(-ysq*ysq-del) * result
   end if
   if (mode .eq. 0) then
      result = 1.0d0 - result
      if (x .lt. 0.0d0)  result = -result
   else
      if (x .lt. 0.0d0)  result = 2.0d0 - result
   end if
end if

return
end subroutine erfcore
