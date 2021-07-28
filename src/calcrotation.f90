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
!     subroutine calcrotation: Calculate rotation defined via rotfrag routine
!       store new coordinates in array for single atom (x,y,z)
!
!     part of QMDFF
!
subroutine calcrotation (x,ori,vec,phi)
implicit none
integer::i,j
real(kind=8)::x(3),ori(3),vec(3),d(3)
real(kind=8)::phi
real(kind=8)::xtmp(3)
real(kind=8)::absd

d(1:3) = vec(1:3)

xtmp(1:3) = x(1:3) - ori(1:3)

absd = sqrt(d(1)**2 + d(2)**2 + d(3)**2)

d(1:3) = d(1:3) / absd

x(1) = ( (d(2)**2+d(3)**2)*cos(phi) + d(1)**2 ) * xtmp(1) + &
  &     ( d(1)*d(2)*(1-cos(phi))-d(3)*sin(phi) ) * xtmp(2) + &
  &     ( d(1)*d(3)*(1-cos(phi))+d(2)*sin(phi) ) * xtmp(3)
x(2) = ( d(1)*d(2)*(1-cos(phi))+d(3)*sin(phi) ) * xtmp(1) + &
  &    ( (d(1)**2+d(3)**2)*cos(phi) + d(2)**2 ) * xtmp(2) + &
  &    ( d(2)*d(3)*(1-cos(phi))-d(1)*sin(phi) ) * xtmp(3)
x(3) = ( d(1)*d(3)*(1-cos(phi))-d(2)*sin(phi) ) * xtmp(1) + &
  &    ( d(2)*d(3)*(1-cos(phi))+d(1)*sin(phi) ) * xtmp(2) + &
  &    ( (d(1)**2+d(2)**2)*cos(phi) + d(3)**2 ) * xtmp(3)
x(1:3) = x(1:3) + ori(1:3)

return
end subroutine calcrotation
