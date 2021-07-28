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
!     subroutine ffhess: calculate the numerical QMDFF hessian based 
!               on the analytical gradients for parameter optimization
!
!     part of QMDFF
!
subroutine ffhess(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy, &
    &    n3,iadr,hessatoms,ffh,ip)
implicit none
integer::n,n3,ip,at(n),iadr(2,n3)
integer::hessatoms(5,*)
integer::i,j,ii,jj,kk,ia,jc,ja,ic,ih
real(kind=8)::xyz(3,n),ffh(n3)
real(kind=8)::step,step2,gl(3,n),gr(3,n),e
real(kind=8)::r0ab(94,94),zab(94,94),r094(94,94)
real(kind=8)::sr42(94,94),c6xy(n,n),q(n)
real(kind=4),allocatable :: h(:,:)               

allocate(h(3*n,3*n))

h=0

step=0.001 
step2=1./(2.*step)
!
!     for the standard case: hessian of QMDFF minimum
!
if (ip.eq.0) then
   do ia = 1, n
      do ic = 1, 3
         ii = (ia-1)*3 + ic

         xyz(ic,ia)=xyz(ic,ia)+step
         call ff_eg (n,at,xyz,e,gr)
         call ff_nonb(n,at,xyz,q,r0ab,zab,r094, &
           &  sr42,c6xy,e,gr)
         call ff_hb(n,at,xyz,e,gr)

         xyz(ic,ia)=xyz(ic,ia)-2.*step

         call ff_eg (n,at,xyz,e,gl)
         call ff_nonb(n,at,xyz,q,r0ab,zab,r094, &
           &  sr42,c6xy,e,gl)
         call ff_hb(n,at,xyz,e,gl)
         xyz(ic,ia)=xyz(ic,ia)+step

         do ja = 1, n
            do jc = 1, 3
               jj = (ja-1)*3 + jc
               h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2
            end do
         end do 
      end do 
   end do
else
!
!     for Levenberg-Marquard derivations after force constants 
!     to fill the matrix aa (called by pderiv-routine)
!
   do ih=1,hessatoms(5,ip)
      ia=hessatoms(ih,ip)
      do ic = 1, 3
         ii = (ia-1)*3 + ic

         xyz(ic,ia)=xyz(ic,ia)+step

         call ff_egh(n,at,xyz,e,gr,ia)

         xyz(ic,ia)=xyz(ic,ia)-2.*step

         call ff_egh(n,at,xyz,e,gl,ia)
         xyz(ic,ia)=xyz(ic,ia)+step

         do ja = 1, n
            do jc = 1, 3
               jj = (ja-1)*3 + jc
               h(ii,jj) =(gr(jc,ja) - gl(jc,ja)) * step2 
            end do
         end do
      end do
   end do
end if
do kk=1,n3
   i=iadr(1,kk)
   j=iadr(2,kk)
   ffh(kk)=(h(i,j)+h(j,i))*0.5
end do

deallocate(h)

return
end subroutine ffhess
