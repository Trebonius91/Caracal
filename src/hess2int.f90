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
!     subroutine hess2int: converts a cartesian hessian matrix to internal 
!     coordinates
!   
!     part of EVB
!

subroutine hess2int(xyz5,internal,H_q,H_x,g_q)
use evb_mod
implicit none
!     The actual structure in cartesian coordinates
real(kind=8), intent(in) :: xyz5(3,natoms)
!     The actual structure in internal coordinates
real(kind=8), intent(in) :: internal(nat6)
!     The actual gradient in internal coordinates
real(kind=8), intent(in) :: g_q(nat6)
!     The actual hessian matrix in cartesian coordinates
real(kind=8), intent(in) :: H_x(3*natoms,3*natoms)
!     The hessian matrix in internal coordinates (result)
real(kind=8), intent(out):: H_q(nat6,nat6)
!     The Wilson B matrix
real(kind=8) :: B_mat(nat6,3*natoms)
!     The transpose of the Wilson B matrix
real(kind=8) :: Bt(3*natoms,nat6)
!     The pseudo-inverted transpose of the Wilson B matrix
real(kind=8) :: Bt_inv(nat6,3*natoms)
!     The pseudo-inverted wilson matrix (transposed tranpose)
real(kind=8) :: B_inv(3*natoms,nat6)
!     Derivative of the Wilson matrix with respect to cartesians
real(kind=8) :: dB_mat(nat6,3*natoms,3*natoms)
!     Product of Wilson derivative and internal gradient
real(kind=8) :: dB_g_q(3*natoms,3*natoms)
!     The difference of cartesian hessian and product matrix
real(kind=8) :: mat_diff(3*natoms,3*natoms)
!     The projector matrix for nonredundant coordinates 
real(kind=8) :: P_mat(nat6,nat6)
!     Loop indices
integer :: i,j,k
!     Number of atoms in structure
integer::nat,nat3

!
!     Restore the number of atoms
!
nat=natoms
nat3=3*nat
!
!     first, build the Wilson B matrix
!

call calc_wilson(xyz5,internal,B_mat)
!
!     Transpose the Wilson Matrix B (--> Bt)
!
do i=1,nat6
   do j=1,3*nat
      Bt(j,i)=B_mat(i,j)
   end do
end do
!
!     Calculate the pseudoinverse of the wilson matrix via Singular value
!     decomposition(SVD)--> Bt^+
!
call pseudoinv(Bt,3*natoms,nat6,Bt_inv)
!
!     Transpose the pseudoinverse for last factor of resulting equation
!
B_inv=transpose(Bt_inv)
!
!     Call subroutine for calculation of cartesian derivative of the 
!     Wilson B-matrix
!

call calc_dwilson(xyz5,internal,dB_mat)


!
!     Calculate product of internal gradient times derivative of Wilson matrix
!

dB_g_q=0
do i=1,nat3
   do j=1,nat3
      do k=1,nat6
         dB_g_q(j,i)=dB_g_q(j,i)+g_q(k)*dB_mat(k,j,i)
      end do
   end do
end do

!
!     Calculate difference of cartesian hessian and the new product matrix
!     (H_x-Bijk*g_q)
!
mat_diff=0
mat_diff=H_x-dB_g_q

!
!     Now, calculate the internal hessian!
!
H_q=0.d0
H_q=matmul(Bt_inv,matmul(mat_diff,B_inv))
!
!     Calculate the projection matrix P in order to remove numerical noise 
!
P_mat=matmul(B_mat,B_inv)
!
!     Do projections with internal hessian
!
H_q=matmul(P_mat,matmul(H_q,P_mat))

!write(*,*) "Internal hessian:"
!do i=1,nat6
!   write(*,'(6e17.8)') H_q(i,:)
!end do

!stop "In hess2inbt!"
return
end subroutine hess2int
