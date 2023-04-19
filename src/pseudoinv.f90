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
!     subroutine pseudoinv: Calculate the pseudoinverse 
!        of a nonquadratic matrix (e.g. for coordinate transformations
!        in grad2int and hess2int) by using Singular Value 
!        Decomposition (SVD) of the input matrix.
!        This routine might also be used for usual matrix
!        inversions of quadratic matrices!
!        for grad2int: N=nat6,M=3*natoms
!
!     part of EVB
!
subroutine pseudoinv(matrix,M,N,inverse)
use evb_mod
implicit none
!     The matrix that shall be inverted
real(kind=8), intent(in) :: matrix(M,N)
!     First/leading dimension of the matrix
integer, intent(in) :: M
!     Second dimension of the matrix
integer, intent(in) :: N
!     The pseudoinverse of the matrix (result)
real(kind=8), intent(out) :: inverse(N,M)
!     Variables for Lapack routine
character(len=1)::JOBU,JOBVT
integer::INFO,LDA,LDU,LDVT,LWORK
real(kind=8)::lam_inv(N,M)
real(kind=8)::A(M,N),S(min(M,N)),U(M,M),VT(N,N)
real(kind=8),allocatable::WORK(:)
!     Transposed of the matrices resulting from SVD
real(kind=8)::V(N,N),Ut(M,M)
!     Loop indices
integer::i,j

!
!     call the lapack routine for SVD
!
JOBU='A'
JOBVT='A'
LDA=M
LWORK=MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N))
allocate(WORK(LWORK))
LDU=M
LDVT=N
A=matrix

call DGESVD(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, &
               & WORK, LWORK, INFO)
!
!     Transpose resulting matrices for pseudoinverse
!
V=transpose(Vt)
Ut=transpose(U)
!
!     fill the inverse singular value matrix
!
lam_inv=0.d0
do i=1,min(M,N)
   lam_inv(i,i)=1.d0/S(i)
   if (S(i) .lt. 1D-7) then
      lam_inv(i,i)=0.d0
   end if
!   write(*,*) i,lam_inv(i,i)
end do
!stop "Hph"
!
!     calculate the pseudoinverse: V^t*lambda^-1*U
!
inverse=matmul(V,matmul(lam_inv,Ut))


return
end subroutine pseudoinv

