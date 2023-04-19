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
!     subroutine grad2int: converts an array of cartesian gradients to internal 
!     coordinates
!
!     part of EVB
!

subroutine grad2int(xyz5,internal,g_q,g_x)
use evb_mod
implicit none
!     The actual structure in cartesian coordinates
real(kind=8), intent(in) :: xyz5(3,natoms)  
!     The actual structure in internal coordinates
real(kind=8), intent(in) :: internal(nat6)
!     The actual gradient in cartesian coordinates
real(kind=8), intent(in) :: g_x(3*natoms)
!     The returned actual gradient in internal coordinates (result)
real(kind=8), intent(out) :: g_q(nat6)
!     The Wilson B matrix
real(kind=8) :: B_mat(nat6,3*natoms)
!     The transpose of the Wilson B matrix
real(kind=8) :: Bt(3*natoms,nat6)
!     The pseudo-inverted transpose of the Wilson B matrix
real(kind=8) :: Bt_inv(nat6,3*natoms)
!     The pseudo-inverted wilson matrix (transposed tranpose)
real(kind=8) :: B_inv(3*natoms,nat6)
!     The projector matrix for nonredundant coordinates 
real(kind=8) :: P_mat(nat6,nat6)
!     Loop indices
integer::i,j
!     Number of atoms in structure
integer::nat
!
!     Restore the number of atoms
!
nat=natoms

!
!  first, build the Wilson B matrix
!
call calc_wilson(xyz5,internal,B_mat)
!
!  Transpose the Wilson Matrix B (--> Bt)
!
!write(*,*) "b_mat ------------------:"
!do i=1,3*nat
!   write(*,'(20f12.5)') B_mat(:,i)
!end do
!write(*,*) "-------------------"
Bt=0
do i=1,nat6
   do j=1,3*nat
      Bt(j,i)=B_mat(i,j)
   end do
end do

!
!    Calculate the pseudoinverse of the wilson matrix via Singular value
!    decomposition(SVD)--> Bt^+
!
call pseudoinv(Bt,3*natoms,nat6,Bt_inv)
!
!     Transpose the pseudoinverse for projection (see below)
!
B_inv=transpose(Bt_inv)
!
!    Transform the gradient to internals via g_q=((Bt)^+)*q_x
!
g_q=matmul(Bt_inv,g_x)
!
!     Calculate the projection matrix P in order to remove numerical noise 
!
P_mat=matmul(B_mat,B_inv)
!
!     Do projections with internal gradient
!
g_q=matmul(P_mat,g_q)




!stop "Hphpi"
return
end subroutine grad2int
