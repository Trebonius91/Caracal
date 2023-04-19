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
!     subroutine int2step: converts geometry step in internal
!        coordinates to cartesian coordinates. Due to the curvilineariy
!        of the internal coordinates, this will be done in an 
!        iterative process
!
!     part of EVB
!

subroutine int2step(xyz5,internal,s_q,s_x)
use evb_mod
implicit none
!     The actual structure in cartesian coordinates
real(kind=8), intent(in) :: xyz5(3,natoms)
!     The actual structure in internal coordinates
real(kind=8), intent(in) :: internal(nat6)
!     The actual geometry step in internal coordinates 
real(kind=8), intent(in) :: s_q(nat6)
!     The actual geometry step in cartesian coordinates (result)
real(kind=8), intent(out) :: s_x(3,natoms)
!     The cartesian coordinates of structure in vector format (start)
real(kind=8) :: xyz_1d(3*natoms)
!     The actual cartesian coordinates during cycle
real(kind=8) :: xyz_3d(3,natoms)
!     The actual internal coordinates during cycle
real(kind=8) :: int_act(nat6),int_old(nat6)
!     Partition of the needed change in internal coordinates
real(kind=8) :: deltq_k(nat6)
!     The cartesian coordinates of structure in vector format (actual)
real(kind=8) :: xyz_act(3*natoms)
!     The Wilson B matrix
real(kind=8) :: B_mat(nat6,3*natoms)
!     The local convergence criterion for this loop
real(kind=8) :: conv_crit
!     Maximal number of cycles that will be done
integer :: max_cycles
!     The pseudo-inverted Wilson B matrix
real(kind=8) :: B_inv(3*natoms,nat6)
!     Loop indices
integer :: i,j,k
!     Number of atoms in structure
integer :: nat
real(kind=8) :: cart_test(3*natoms)
!
!     Restore the number of atoms
!
nat=natoms
!
!     The local convergence criterion and maximal number of cycles for 
!     this loop
!
conv_crit=1E-13
max_cycles=50
!
!     Convert (3,nat) to (3*nat) vector format
!
do i=1,nat
   do j=1,3
      xyz_act((i-1)*3+j)=xyz5(j,i)
   end do 
end do
xyz_1d=xyz_act

!
!     first, build the Wilson B matrix
!
call calc_wilson(xyz5,internal,B_mat)
!
!     Build the inverse of the Wilson B matrix
!
call pseudoinv(B_mat,nat6,3*natoms,B_inv)
!
!     Calculate cartesian step
!

xyz_act=xyz_act+matmul(B_inv,s_q)
! 
!     Enter big loop till convergence 
!
k=0
do
!
!     Increment actual number of cycles
!
   k=k+1
!
!     Change formate of cartesian coordinates
!
   do i=1,nat
      do j=1,3
         xyz_3d(j,i)=xyz_act((i-1)*3+j)
      end do
   end do
!
!     Calculate internal coordinates of updated cartesian structure
!
   call xyz_2int(xyz_3d,int_act,natoms)
   
   deltq_k=s_q-(int_act-internal)
!
!     If absolute of new step correction is small enough, leave loop!
!
   if (dot_product(deltq_k,deltq_k) .le. conv_crit) exit
   if (k .gt. max_cycles) exit
!
!     Build the Wilson B matrix
!
   call calc_wilson(xyz_3d,int_act,B_mat)
!
!     Build the inverse of the Wilson B matrix
!
   call pseudoinv(B_mat,nat6,3*natoms,B_inv)
!
!     Calculate cartesian step
!
   xyz_act=xyz_act+matmul(B_inv,deltq_k)
end do
!
!     If convergence was reached: Determine total step in just building the 
!     difference of initial and final xyz structure
!
do i=1,nat
   do j=1,3
      xyz_3d(j,i)=xyz_act((i-1)*3+j)
   end do
end do


s_x=xyz_3d-xyz5

return
end subroutine int2step
