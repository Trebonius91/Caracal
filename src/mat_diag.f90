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
!     subroutine mat_diag: diagonalize the DG-EVB coefficient matrix 
!        during coupling term optimizations
!
!     part of EVB
!
subroutine mat_diag(mat_size)
use evb_mod
character(len=1)::JOBA,JOBU,JOBV,JOBR,JOBT,JOBP
integer::M,N,LDA,LDU,LDV,LWORK,INFO
integer,dimension(:),allocatable::IWORK,IPIV
real(kind=8),dimension(:),allocatable::SVA,WORK,B,test
real(kind=8),dimension(:,:),allocatable::A,U,V
real(kind=8),dimension(:,:),allocatable::sigma,Ut,Vt,A_inv
integer::NRHS,LDB
!
!   for GMRES solver
!
integer(kind=4)::nz_num,itr_max,mr
real(kind=8),allocatable::matval(:),mat_buffer(:)
integer(kind=4),allocatable:: ia(:),ja(:),j_buffer(:),i_buffer(:)
real(kind=8),allocatable::x(:),w(:)
real(kind=8)::tol_abs,tol_rel

!
!   allocate parameters for lapck matrix diagonalizations
!

!   STANDARD routine
N=mat_size
NRHS=1
LDA=N
LDB=N
allocate(A(N,N))
allocate(B(N))
allocate(IPIV(N))
!   SINGULAR VALUE DECOMPISITION 

JOBA="C"
JOBU="F"
JOBV="V"
JOBR="N"
JOBT="T"
JOBP="N"
M=mat_size
LDU=N
LDV=N
LWORK=max(2*M+N,6*N+2*N*N)

allocate(IWORK(M+3*N))
allocate(SVA(N))
allocate(WORK(LWORK))
allocate(U(N,N))
allocate(V(N,N))
allocate(sigma(mat_size,mat_size))
allocate(Ut(mat_size,mat_size),Vt(mat_size,mat_size))
allocate(A_inv(mat_size,mat_size))

!
! For small matrices, solve the system of equations with the standard
! Lapack-routine
!
B=f_vec
A=d_mat
if (mat_size .le. 10000) then
   call dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
   if (INFO .ne. 0) then
      write(*,*) "ERROR!"
      write(*,*) "THE System of Equations F=D*X could not be solved!"
      write(*,*) "Please retry it with another set of reference informations"
      write(*,*) "The INFO value was",INFO
      write(*,*) "Retry new local optimization run..."
!      call fatal
   end if
 b_vec=B

!  --------------------------------------
!    SINGULAR VALUE DECOMPOSITION for larger matrices
!  For further information, see above on
!  http://www.netlib.org/lapack/explore-html/d1/d7e/group__double
!     _g_esing.html#gad8e0f1c83a78d3d4858eaaa88a1c5ab1    
!
!  For the theoretical foundation, see 
!  http://matheplanet.com/default3.html?call=article.php
!      ?sid=742&ref=https%3A%2F%2Fwww.google.de
!
!   ---> A = U * Sigma * V^T

else if (mat_size .eq. 0) then
   A=d_mat

   call dgejsv(JOBA,JOBU,JOBV,JOBR,JOBT,JOBP,M,N,A,LDA,SVA,U,&
              &LDU,V,LDV,WORK,LWORK,IWORK,INFO)
   if (INFO .ne. 0) then
      write(15,*) "ERROR!"
      write(15,*) "THE System of Equations F=D*X could not be solved!"
      write(15,*) "Please retry it with another set of reference informations"
      write(15,*) "The INFO value was",INFO
      write(*,*) "Retry new local optimization run..."
!      call fatal
   end if

!
!   The singular vaules Sigma
!
   do i=1,mat_size
      sigma(i,i)=1/SVA(i)
   end do
!
!   Transpose matrices
!

   do i=1,mat_size
      do j=1,mat_size
         Ut(i,j)=U(j,i)
      end do
   end do

!
!   Calculate Pseudoinversion: 
!   ---> A^+ = V * Sigma^+ U^T 
! 

   A_inv=matmul(V,sigma)
   A_inv=matmul(A_inv,Ut)

!
!   B = A^+ * F
!

   b_vec=matmul(A_inv,f_vec)

!  --------------------------------------
else
!
!   Call the subroutine for the GMRES solver
!
allocate(j_buffer(mat_size*mat_size))
allocate(i_buffer(mat_size*mat_size))
allocate(mat_buffer(mat_size*mat_size))
itr_max=100
mr=mat_size
tol_abs=1D-08
tol_rel=1D-08
!
!  Determine number of nonzero elements:
!
nz_num=0
do i=1,mat_size   ! row
   do j=1,mat_size   ! column
      if (abs(d_mat(i,j)) .ge. 1D-10) then
         nz_num=nz_num+1
         j_buffer(nz_num)=j
         i_buffer(nz_num)=i
!         write(*,*) i,j,nz_num,d_mat(i,j)
         mat_buffer(nz_num)=d_mat(i,j)
      end if 
   end do
end do
!write(*,*) "NZ",nz_num,mat_size
allocate(ja(nz_num))
allocate(ia(nz_num))
allocate(matval(nz_num))
do i=1,nz_num
   ja(i)=j_buffer(i)
   ia(i)=i_buffer(i)
   matval(i)=mat_buffer(i)
end do
!write(*,*) "Number of nonzero values:",nz_num
!write(*,*) "ia",ia
!write(*,*) "ja",ja
!write(*,*) matval
!b_vec(1)=-1
!b_vec(2)=2
open(unit=45,file="gmres.log",status="unknown")
call mgmres_st(mat_size,nz_num,ia,ja,matval,b_vec,f_vec,itr_max,mr,&
              &tol_abs,tol_rel)
close(45)
!write(*,*) "after call"
!write(*,*) b_vec
end if
deallocate(A,B,IPIV)
deallocate(IWORK)
deallocate(SVA)
deallocate(WORK)
deallocate(U,V)
deallocate(sigma)
deallocate(Ut,Vt)
deallocate(A_inv)

return
end subroutine mat_diag

