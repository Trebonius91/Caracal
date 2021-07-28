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
!     subroutine moment: calculate moment for the reference structure from 
!           energy and gradient
!
!     part of QMDFF
!

subroutine moment(n,xyz,g,f1,f2)
implicit none
integer::i,n
real(kind=8)::xyz(3,n)
real(kind=8)::g(3,n)
real(kind=8)::f1(3),f2(3)

f1=0
f2=0
do i=1,n
   f1(1) = f1(1)+ g(1,i)
   f1(2) = f1(2)+ g(2,i)
   f1(3)=  f1(3)+ g(3,i)
   f2(1) = f2(1) + xyz(2,i)*g(3,i) - xyz(3,i)*g(2,i)
   f2(2) = f2(2) + xyz(3,i)*g(1,i) - xyz(1,i)*g(3,i)
   f2(3) = f2(3) + xyz(1,i)*g(2,i) - xyz(2,i)*g(1,i)
end do

return
end subroutine moment

