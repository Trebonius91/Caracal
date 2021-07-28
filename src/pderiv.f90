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
!     subroutine pderiv: calculate the numerical derivative of the 
!                 QMDFF hessian for the force constants (aa)
!
!     part of QMDFF
!
subroutine pderiv(i1,i2,n,at,xyzin,q,r0ab,zab,r094,sr42,c6xy, &
   &   n3,np,iadr)
implicit none
integer::n,n3,at(n),iadr(2,n3),np
real(kind=8)::xyzin(3,n)
real(kind=8)::r0ab(94,94),zab(94,94),r094(94,94)
real(kind=8)::sr42(94,94),c6xy(n,n),q(n)
integer::i1,i2,i,j,ii,k,TID,OMP_GET_THREAD_NUM
real(kind=8)::step,step2
character(len=20)::atmp
real(kind=8),allocatable::xyz(:,:)
real(kind=8),allocatable::p(:)
real(kind=8),allocatable::rr(:),ll(:)
real(kind=8),allocatable::aa(:)
integer,allocatable::hessatoms(:,:)

allocate(xyz(3,n),p(np),aa(n3),hessatoms(5,np),rr(n3),ll(n3))
!
!     store the matrix aa into an unformatted file (qmdfftmp..)!
!

TID = OMP_GET_THREAD_NUM()
if (tid.lt.100) write(atmp,'(''qmdfftmp.'',i2)') tid
if (tid.lt.10 ) write(atmp,'(''qmdfftmp.'',i1)') tid

open (unit=44+tid,file=atmp)!,form='unformatted')

call getpar(p,.false.,hessatoms)
!
!     elongate all force constants and calculate the derivative matrix aa
!
xyz = xyzin

do ii=i1,i2
   step=abs(0.01*p(ii))
   step=max(step,0.00001d0)

   p(ii)=p(ii)+step
   call putpar(p,.false.)
   call ffhess(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,n3,iadr, &
      &   hessatoms,rr,ii)
   p(ii)=p(ii)-2.*step
   call putpar(p,.false.)
   call ffhess(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,n3,iadr, &
      &   hessatoms,ll,ii)
   p(ii)=p(ii)+step
   call putpar(p,.false.)

   step2=1.0d0/(2.*step)
   do j=1,n3
      aa(j)=(rr(j)-ll(j))*step2
   end do
!   write(*,*) aa
   write(44+tid,*) aa
end do

close(44+tid)

deallocate(xyz,p,aa,hessatoms,rr,ll)

return
end subroutine pderiv
