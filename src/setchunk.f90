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
!     "setchunk" marks a chunk in the PME spatial table which is
!     overlapped by the B-splines for an electrostatic site
!  
!     part of QMDFF
!
subroutine setchunk (i,cid,off1,off2,off3,pmetable,n)
implicit none
integer::i,k,n
integer::off1,off2,off3
integer::cid(3),temp(3)
integer::pmetable(1,n)
!
!
!     mark neighboring chunk overlapped by an electrostatic site
!
temp(1) = cid(1) + off1
if (temp(1) .lt. 1)  temp(1) = 1
if (temp(1) .gt. 1)  temp(1) = 1
temp(2) = cid(2) + off2
if (temp(2) .lt. 1)  temp(2) = 1
if (temp(2) .gt. 1)  temp(2) = 1
temp(3) = cid(3) + off3
if (temp(3) .lt. 1)  temp(3) = 1
if (temp(3) .gt. 1)  temp(3) = 1
k = (temp(3)-1)*1 + (temp(2)-1)*1 + temp(1)
pmetable(i,k) = 1
return
end subroutine setchunk

