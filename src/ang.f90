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
!    function ang: calculate angle between at1-at2-at3
!
!    part of EVB
!
function ang(at1,att2,at3,xyz5)
use evb_mod
integer::at1,att2,at3
integer::i,j
real(kind=8)::ang
real(kind=8)::uvec(3),vvec(3)
real(kind=8)::angvec(3)
real(kind=8)::uabs,vabs,vlen
real(kind=8)::xyz5(3,natoms)
real(kind=8)::ax,ay,az,bx,by,bz
!
!    ATTENTION: ANGLE CALCULATION CHANGED 21.03.2023
!
!    calculate bond vectors of both bonds 
!  
!uvec(:)=xyz5(:,at1)-xyz5(:,att2)
!vvec(:)=xyz5(:,at3)-xyz5(:,att2)
!
!    calculate lengths of these bonds 
!
!uabs=vlen(uvec)
!vabs=vlen(vvec)
! TEST
!ang=dot_product(uvec/uabs,vvec/vabs)


!
!    calculate the angle coordinate
!

!angvec=0.5d0*(uvec(:)/uabs-vvec(:)/vabs)

!ang=vlen(angvec)

ax=xyz5(1,at1)-xyz5(1,att2)
ay=xyz5(2,at1)-xyz5(2,att2)
az=xyz5(3,at1)-xyz5(3,att2)
bx=xyz5(1,at3)-xyz5(1,att2)
by=xyz5(2,at3)-xyz5(2,att2)
bz=xyz5(3,at3)-xyz5(3,att2)
ang=acos((ax*bx+ay*by+az*bz)/(sqrt(ax*ax+ay*ay+az*az)*&
          &sqrt(bx*bx+by*by+bz*bz)))


return
end function ang
