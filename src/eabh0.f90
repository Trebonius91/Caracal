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
!     subroutine eabh0: calculate only energies for hydrogen/halogem bonds
!
!     part of QMDFF
!

subroutine eabh0(n,A,B,H,rab,xyz,ca,cb,energy)
use qmdff
implicit none

integer::A,B,H,n
real(kind=8)::xyz(3,n),energy,rab
real(kind=8)::ang,xy,cosabh,d2ij,d2ik,d2jk,term,temp,fa,fb
real(kind=8)::aterm,dampm,dampl,xm,ym,zm,rhm,rab2,ca,cb,da
REAL(kind=8)::longcut=8.
REAL(kind=8)::alp7=12
REAL(kind=8)::alp3= 6
real(kind=8)::rbond(3)
!
!     AB distance
!
rab2=rab*rab
!
!     long  damping
!
dampl=1./(1.+(rab/longcut)**alp7)
!
!     cos angle A-B and H (or B-H H) is term
!
rbond(1)=XYZ(1,A)-XYZ(1,H)
rbond(2)=XYZ(2,A)-XYZ(2,H)
rbond(3)=XYZ(3,A)-XYZ(3,H)
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(rbond)
end if

D2IK=rbond(1)**2+rbond(2)**2+rbond(3)**2

rbond(1)=XYZ(1,H)-XYZ(1,B)
rbond(2)=XYZ(2,H)-XYZ(2,B)
rbond(3)=XYZ(3,H)-XYZ(3,B)
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(rbond)
end if

D2JK=rbond(1)**2+rbond(2)**2+rbond(3)**2


if (d2ik.gt.d2jk) then
     XY = SQRT(rab2*D2JK)
     term = 0.5D0 * (rab2+D2JK-D2IK) / XY
else
     XY = SQRT(rab2*D2IK)
     term = 0.5D0 * (rab2+D2IK-D2JK) / XY
end if
!
!     angle damping term
!
aterm = (0.5d0*(term+1.0d0))**alp3
!
!     donor-acceptor term
!
fa=d2jk/d2ik
fb=d2ik/d2jk
da=(ca*fb+cb*fa)/(fa+fb)
energy=-da*dampl*aterm/(rab2*rab)

return
end subroutine eabh0
