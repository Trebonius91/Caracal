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
!     subroutine nneighbor: special for EHT torsional profile: saturate bonds
!
!     part of QMDFF
!
subroutine nneighbor(natoms,xyz,iz,wbo,nb)
use qmdff
implicit none
integer::iz(*),natoms,nb(50,2)
real(kind=8)::xyz(3,*),wbo(100,100)

logical::da
integer::iat,i,j,k,nn,ni
real(kind=8)::dx,dy,dz,r,damp,xn,rr,rco,r2,f,aa1

do i=1,natoms
   f=1.3
   k=0
   100 continue
   do iat=1,natoms
      if(iat.ne.i)then
!
!     determine distance and compare with Wiberg-Mayer bond order
!
         dx=xyz(1,iat)-xyz(1,i)
         dy=xyz(2,iat)-xyz(2,i)
         dz=xyz(3,iat)-xyz(3,i)
         r2=dx*dx+dy*dy+dz*dz
         r=sqrt(r2)*0.52917726
         rco=rad(iz(i))+rad(iz(iat))
         if (r.lt.f*rco.and.wbo(iat,i).gt.0.5) then
            k=k+1
            nb(i,2)=iat
         end if
      end if
   end do
   if (k.lt.1.and.f.lt.1.5) then
      f=f*1.1
      goto 100
   end if
   nb(i,1)=k
end do

return
end subroutine nneighbor
