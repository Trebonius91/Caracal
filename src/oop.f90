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
!     function oop: calculate out of plane angle for at1 out of 
!        plane on at3 and at2,at4
!     The local situation needs to be (4 in the middle):
!     1     3
!       \ /
!        4
!        |
!        2
!
!     part of EVB
!
function oop(at1,att2,at3,at4,xyz5)
use evb_mod
implicit none
integer::at1,att2,at3,at4
integer::i,j
real(kind=8)::oop
real(kind=8)::v41(3),v42(3),v43(3)
real(kind=8)::v41_len,v42_len,v43_len
real(kind=8)::v4142(3),v4243(3),v4341(3)
real(kind=8)::normvec(3)
real(kind=8)::oopnorm,vlen
real(kind=8)::xyz5(3,natoms)
!
!     calculate bond vectors of involved bonds 
!
v41=xyz5(:,at4)-xyz5(:,at1)
v42=xyz5(:,at4)-xyz5(:,att2)
v43=xyz5(:,at4)-xyz5(:,at3)
!
!     calculate bond lengths 
!
v41_len=vlen(v41)
v42_len=vlen(v42)
v43_len=vlen(v43)
!
!     normalize bond vectors 
!
v41=v41/v41_len
v42=v42/v42_len
v43=v43/v43_len
!
!     calculate norm vector for coordinate (cross products)
!
call crossprod(v41,v42,v4142)
call crossprod(v42,v43,v4243)
call crossprod(v43,v41,v4341)

normvec=v4142+v4243+v4341

oopnorm=vlen(normvec)
!
!     calculate the oop coordinate (normalized triple product)
!
oop=dot_product(v41,v4243)/oopnorm
!oop=sqrt(3d0-cosphi1**2-cosphi2**2-cosphi3**2+2d0*(cosphi1*(cosphi2-1d0)+&
!             & cosphi2*(cosphi3-1d0)+cosphi3*(cosphi1-1d0)))

return
end function oop
