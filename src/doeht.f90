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
!     function doeht: Find cases where EHT does not yield reliable results
!
!     part of QMDFF
!

logical function doeht(n,at,i,j,r,hyb,cring,ringsize,wbo)
use qmdff
implicit none
integer::n,at(n),hyb(n),i,j
integer::cring(8,n),ringsize(n)
real(kind=8)::wbo(n,n),r
logical::samering

doeht=.true.
!
!     both atoms are in the same ring
!
if (samering(n,i,j,cring,ringsize)) doeht=.false.
!
!     bond order > some cutoff
!
if (wbo(j,i) .gt. 1.5) doeht=.false.
!
!     double bond
!
if (hyb(i).eq.2.and.hyb(j).eq.2) return
!
!     e.g. Si-Si barriers are not good with TB so take the Hessian
!     because this might hold for many long bonds we add the distance
!     as a criterion      
!
if (r.gt.4) doeht=.false.
!
!     partially conjugated 
!
if (ringsize(i).gt.0.and.wbo(j,i).gt.1.3) doeht=.false.
if (ringsize(j).gt.0.and.wbo(j,i).gt.1.3) doeht=.false.
!
!     metal atoms contained?  
! 
if (metal(at(i)).eq.1.or.metal(at(j)).eq.1) doeht=.false.

return
end function doeht
