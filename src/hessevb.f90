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
!     subroutine hessevb: calculates the EVB-QMDFF hessian matrix of an input-structure
!     numerically
!     

subroutine hessevb(xyz2,h)
use general
use evb_mod

implicit none
integer mode
real(kind=8)::xyz2(3,natoms)
real(kind=8)::h(3*natoms,3*natoms)
real(kind=8)::freq(3*natoms)
real(kind=8)::outval
real(kind=8)::ams(107)
real(kind=8)::step,amu2au,au2cm,zpve,e,e_two,dumi,dumj,zpve2,vmad,hstep
integer::n3,i,j,k,ic,jc,ia,ja,ii,jj,info,lwork
integer::mincol,maxcol,lin,bvib,kk
real(kind=8),allocatable::hs(:)
real(kind=8),allocatable::aux (:)
real(kind=8),allocatable::isqm(:)
real(kind=8),allocatable::gr  (:,:),gr_two(:,:)
real(kind=8),allocatable::gl  (:,:),gl_two(:,:)
real(kind=8),allocatable::xsave(:,:)
character(len=3)::sym   
character(len=10)::a10  
character(len=80)::fname

n=natoms
n3=3*n

allocate(gr(3,n),gl(3,n),isqm(n3),xsave(3,n))

xsave = xyz2

! step length
step=0.0000001
!
! compute the EVB-QMDFF-Hessian by moving the coordinates 
!
do ia = 1, natoms
   do ic = 1, 3
      ii = (ia-1)*3+ic
      xyz2(ic,indi(ia))=xyz2(ic,indi(ia))+step

      call gradient(xyz2,e,gr,1)
      xyz2(ic,indi(ia))=xyz2(ic,indi(ia))-2.*step

      call gradient(xyz2,e,gl,1)
      xyz2(ic,indi(ia))=xyz2(ic,indi(ia))+step
      do ja = 1, natoms
         do jc = 1, 3
            jj = (ja-1)*3 + jc
            h (ii,jj) =(gr(jc,indi(ja)) - gl(jc,indi(ja))) / (2.*step)
         end do
      end do 
   end do 
end do

xyz2=xsave
return
end subroutine hessevb
