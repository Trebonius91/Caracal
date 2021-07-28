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
!     subroutine bangle: calculates a bond angle between I-J-K
!
!     part of QMDFF
!
subroutine bangle(xyz,i,j,k,angle)
implicit none
real(kind=8)::xyz(3,*),angle
integer::i,j,k
real(kind=8)::d2ij,d2jk,d2ik,xy,temp

d2ij = (xyz(1,i)-xyz(1,j))**2+ &
  &    (xyz(2,i)-xyz(2,j))**2+ &
  &    (xyz(3,i)-xyz(3,j))**2
d2jk = (xyz(1,j)-xyz(1,k))**2+ &
  &    (xyz(2,j)-xyz(2,k))**2+ &
  &    (xyz(3,j)-xyz(3,k))**2
d2ik = (xyz(1,i)-xyz(1,k))**2+ &
  &    (xyz(2,i)-xyz(2,k))**2+ &
  &    (xyz(3,i)-xyz(3,k))**2
xy = sqrt(d2ij*d2jk)
temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
if (temp .gt. 1.0d0)  temp= 1.0d0
if (temp .lt. -1.0d0) temp=-1.0d0
angle = acos( temp )

end subroutine bangle
