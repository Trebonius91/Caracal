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
!     subroutine qmdffhess: calculates the hessian matrices of an input-structure
!     This is done for the two QMDFFs 
!     
!     part of EVB
!

subroutine hessqmdff(xyz2,h,h2)
use general
use evb_mod

implicit none
integer mode
real(kind=8)::xyz2(3,n)
real(kind=8)::h(3*natoms,3*natoms),h2(3*natoms,3*natoms)
real(kind=8)::freq(3*n)
real(kind=8)::outval
real(kind=8)ams(107)
real(kind=8)step,amu2au,au2cm,zpve,e,e_two,dumi,dumj,zpve2,vmad,hstep
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
allocate(gr_two(3,n),gl_two(3,n))

xsave = xyz2

! step length
step=0.00001d0
!
! compute the first QMDFF-Hessian by moving the coordinates 
!
do ia = 1, n
   do ic = 1, 3
      ii = (ia-1)*3+ic

      xyz2(ic,ia)=xyz2(ic,ia)+step

      call ff_eg(n,at,xyz2,e,gr)
      call ff_nonb(n,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e,gr)
      call ff_hb(n,at,xyz2,e,gr)

      xyz2(ic,ia)=xyz2(ic,ia)-2.*step

      call ff_eg(n,at,xyz2,e,gl)
      call ff_nonb(n,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e,gl)
      call ff_hb(n,at,xyz2,e,gl)
      xyz2(ic,ia)=xyz2(ic,ia)+step
      do ja = 1, n
         do jc = 1, 3
            jj = (ja-1)*3 + jc
            h (ii,jj) =(gr(jc,ja) - gl(jc,ja)) / (2.*step)
         end do
      end do 

   end do 
end do
!
!     compute the second QMDFF-Hessian by moving the coordinates 
!

xyz2 = xsave 

do ia = 1, n
   do ic = 1, 3
      ii = (ia-1)*3+ic

      xyz2(ic,ia)=xyz2(ic,ia)+step

      call ff_eg_two(n,at,xyz2,e_two,gr_two)
      call ff_nonb_two(n,at,xyz2,q,r0ab,zab,r094_mod, &
        &     sr42,c6xy,e_two,gr_two)
      call ff_hb_two(n,at,xyz2,e_two,gr_two)

      xyz2(ic,ia)=xyz2(ic,ia)-2.*step

      call ff_eg_two(n,at,xyz2,e_two,gl_two)
      call ff_nonb_two(n,at,xyz2,q,r0ab,zab,r094_mod, &
        &     sr42,c6xy,e_two,gl_two)
      call ff_hb_two(n,at,xyz2,e_two,gl_two)

      xyz2(ic,ia)=xyz2(ic,ia)+step

      do ja = 1, n
         do jc = 1, 3
            jj = (ja-1)*3 + jc
            h2 (ii,jj) =(gr_two(jc,ja) - gl_two(jc,ja)) / (2.*step)
         end do
      end do

   end do
end do

return
end subroutine hessqmdff
