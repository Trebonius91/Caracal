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
!     subroutine blckmgs: Modified Gram-Schmidt orthonormalization 
!               of a real matrix
!
!     part of QMDFF
!
subroutine blckmgs(m,n,ndim,darray)
implicit none

integer,intent(in)::m,n,ndim
real(kind=8),dimension(ndim,n),intent(out)::darray
integer::ii,jj,kk,ll,ibsize,nblcks,istrt,jstrt,iend,ncol,ierr
real(kind=8)::tmp
real(kind=8),dimension(:,:),allocatable::smat
real(kind=8)::thr
real(kind=8)::ddot

thr = epsilon(1.0d0) !Threshold for zero vectors
!
!     Block size optimized for Athlon 1200 MHz with 2.0GB memory for 
!     matrices up to 5000x5000
ibsize = 60
!
!     Allocate overlap matrix
!
allocate(smat(ibsize,ibsize),stat=ierr)
if(ierr /= 0)  stop 'Mamory allocation error in blckmgs'
!
!     Calculate the number of blocks
!
nblcks = (n+ibsize-1)/ibsize
ibsize = min(n,ibsize)
!
!     Orthogonalize the first block using modified schmidt
!
do ii=1,ibsize

   tmp = ddot(m,darray(1,ii),1,darray(1,ii),1)
!
!     Linear dependence
!
   if (tmp < thr) then
      darray(1:m,ii) = 0.0d0
      cycle
   end if

   tmp = 1.0d0/sqrt(tmp)
   call dscal(m,tmp,darray(1,ii),1)

   do jj=ii+1,ibsize
      tmp = ddot(m,darray(1,ii),1,darray(1,jj),1)
      call daxpy(m,-tmp,darray(1,ii),1,darray(1,jj),1)
   end do
end do
!
!     Loop over remaining blocks
!
do ii=1,nblcks-1
!
!     Initial and final column and number of columns in the block ii+1
!
   istrt = ii*ibsize+1
   iend  = (ii+1)*ibsize
   iend  = min(n,iend)
   ncol  = iend - istrt + 1
!
!     Orthogonalize the block ii+1 against the previous ones
!
   do jj=1,ii
!
!     Initial index of the block jj
!
      jstrt = (jj-1)*ibsize+1
!
!     multiply the matrices/vectors
!
      call dgemm('t','n',ibsize,ncol,m,1.0d0,darray(1,jstrt),ndim,darray(1,istrt),ndim,0.0d0,smat,ibsize)
      call dgemm('n','n',m,ncol,ibsize,-1.0d0,darray(1,jstrt),ndim,smat,ibsize,1.0d0,darray(1,istrt),ndim)
   end do
!
!     Othogonalize vectors on the block ii+1 among themself 
!       using modified schmidt
!
   do kk=istrt,iend
      tmp = ddot(m,darray(1,kk),1,darray(1,kk),1)
!
!     Linear dependence
!
      if(tmp < thr) then
         darray(1:m,kk) = 0.0d0
         cycle
      end if
      tmp = 1.0d0/sqrt(tmp)
      call dscal(m,tmp,darray(1,kk),1)
      do ll=kk+1,iend
         tmp = ddot(m,darray(1,kk),1,darray(1,ll),1)
         call daxpy(m,-tmp,darray(1,kk),1,darray(1,ll),1)
      end do
   end do
end do
!
!     Clean up
!
deallocate(smat,stat=ierr)
if(ierr /= 0) stop 'Mamory deallocation error in blckmgs'

return
end subroutine blckmgs
