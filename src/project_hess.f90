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
!     subroutine project_hess: This routine takes a hessian matrix in 
!     cartesian coordinates, performs a diagonalization and removes all 
!     negative eigenvalues of it. Then the matrix is transformed back
!     into its original form to serve as a corrected hessian.
!     This is needed for RP-EVB in order to avoid divergent behavior 
!     of modes with negative eigenvalues.
!     The hessian is taken as input and the modified hessian will be 
!     transformed back..
!

subroutine project_hess(hess)
use evb_mod
use general
implicit none 
integer::i,j  ! loop index
real(kind=8)::hess(nat6,nat6)
!real(kind=8)::testmat(3*natoms,3*natoms)
! for matrix diagonalization
character(len=1)::JOBZ,UPLO
integer::M,Nn,LDA,LDU,LDV,LWORK,INFO
integer,dimension(:),allocatable::ipiv
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::S_mat,coord
real(kind=8),dimension(:,:),allocatable::lam_mat ! diagonal matrix 
real(kind=8),dimension(:,:),allocatable::S_inv
real(kind=8),dimension(:,:),allocatable::U
integer::NRHS,LDB



!
!     First, diagonalize the matrix to obtain the eigenvalues and eigenvectors
!

JOBZ='V' !eigenvalues and eigenvectors(U)
UPLO='U' !upper triangle of a
Nn=nat6!3*natoms
LDA=Nn
INFO=0
LWORK=3*Nn-1

allocate(S_mat(Nn,Nn))
allocate(S_inv(Nn,Nn))
allocate(lam_mat(Nn,Nn))
allocate(W(Nn))
allocate(ipiv(Nn))
allocate(WORK(LWORK))
allocate(coord(nats,3))
S_mat=hess

call DSYEV(JOBZ,UPLO,Nn,S_mat,LDA,W,WORK,LWORK,INFO)
if (info .ne. 0) then
   write(*,*) "The diagonalization of the hessian matrix in project_hess.f90 failed!"
   call fatal
end if

!
!     Second, calculate the inverse of the eigenvector matrix
!
S_inv = S_mat

!
!    DGETRF computes an LU factorization of a general M-by-N matrix A
!    using partial pivoting with row interchanges.
!
deallocate(WORK)
allocate(WORK(Nn))
call DGETRF(Nn, Nn, S_inv, Nn, ipiv, info)

if (info /= 0) then
   write(*,*) 'The eigenvector Matrix in project_hess.f90 is numerically singular!'
   call fatal
end if
!
!    DGETRI computes the inverse of a matrix using the LU factorization
!    computed by DGETRF.
!
call DGETRI(Nn, S_inv, Nn, ipiv, work, Nn, info)

if (info .ne. 0) then
   write(*,*) "Inversion of eigenvector Matrix in project_hess.f90 failed!"
   call fatal
end if


!
!     Build the matrix of eigenvalues Lambda
!     Place the eigenvalues on the diagonals, set all negative eigenvalues to zero
!
lam_mat=0.d0
do i=1,nat6!3*natoms
   do j=1,nat6!3*natoms
      if (i .eq. j) then
         if (W(i) .ge. 0.d0) then
            lam_mat(i,j)=W(i)!/2.d0
         else
            lam_mat(i,j)=0d0!W(i) 
!            lam_mat(i,j)=0.d0 !W(i)
         end if
      end if
   end do
end do

!
!    Calculate the corrected hessian by basis transformation:
!    H_new=S*Lamda*S^-1
!
hess=matmul(S_mat,lam_mat)
hess=matmul(hess,S_inv)

!
!    TEST: look if new hessian has the desired eigenvalues
!
!S_mat=hess
!
!deallocate(WORK)
!allocate(WORK(LWORK))
!call DSYEV(JOBZ,UPLO,Nn,S_mat,LDA,W,WORK,LWORK,INFO)
!if (info .ne. 0) then
!   write(*,*) "The diagonalization of the hessian matrix in project_hess.f90 failed!"
!   call fatal
!end if
!write(155,*) "Eigenvalues:"
!do i=1,nat6
!   write(155,*) W(i)
!end do
!write(199,*) W
return
end subroutine project_hess
