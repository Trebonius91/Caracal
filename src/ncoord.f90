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
!     ncoord: compute coordination numbers by adding an 
!     inverse damping function
!
!     part of QMDFF
!

subroutine ncoord(natoms,rcov,iz,xyz,cn,cn_thr)
implicit none
real(kind=8)::k1,k3
parameter (k1     =16)
parameter (k3     =-4)
integer::iz(*),natoms,i,max_elem
real(kind=8)::xyz(3,*),cn(*),rcov(94),input
real(kind=8)::cn_thr

integer::iat
real(kind=8)::dx,dy,dz,r,damp,xn,rr,rco,r2

do i=1,natoms
   xn=0.0d0
   do iat=1,natoms
      if (iat.ne.i)then
         dx=xyz(1,iat)-xyz(1,i)
         dy=xyz(2,iat)-xyz(2,i)
         dz=xyz(3,iat)-xyz(3,i)
         r2=dx*dx+dy*dy+dz*dz
!     go back to start of loop if number is too big
         if (r2.gt.cn_thr) cycle
         r=sqrt(r2)
!
!     covalent distance in Bohr between all atoms
!
         rco=rcov(iz(i))+rcov(iz(iat))
         rr=rco/r
!
!     counting function: exponential has a better long-range 
!     behavior than MHGs inverse damping
!     decides if a bond is built or not
!
         damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
         xn=xn+damp
       end if
   end do
   cn(i)=xn
end do

end subroutine ncoord
