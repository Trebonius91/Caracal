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
!     subroutine umbrella: Calculate the actual value of the 
!       sampling coordinate for a structure and compare it with 
!       the reference value(s) for this sampling, then apply 
!       bias potential and gradient correction
!
!     part of EVB
!
subroutine umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,grad_xyz,mode)
use general
use evb_mod
implicit none
integer::i,j,i1,i2,j2,k ! loop indices
integer::mode  ! which type of calculation has invoked this subroutine:
               ! 0: usual umbrella sampling
               ! 1: constraint: calculate restrain xi values
               ! 2: pre equilibrate for tri-/tetramolecular reactions
               ! 3: child trajectory: calculate usual xi, no bias applied
real(kind=8)::centroid(3,natoms)  !  the centroid with all the 
                                  !  com's of all beads
real(kind=8)::int_ideal(nat6),int_real(nat6),int_tmp(nat6)
real(kind=8)::xi_ideal,xi_real
real(kind=8)::act_coord(3,natoms)
real(kind=8)::v_add,delta
real(kind=8)::pi
real(kind=8)::int_grad(nat6)
real(kind=8)::grad_xyz(3,natoms,nbeads),grad_bias(3,natoms)
real(kind=8)::grad_xyz2(3*natoms),grad_bias2(3*natoms)
!     for pre equilibration addition of bias potentials
real(kind=8)::r_eds(sum_eds,sum_eds) ! the actual distane between two educts
real(kind=8)::Red(3,sum_eds,sum_eds)  ! all possible educt-educt distances
real(kind=8)::com(sum_eds,3)
real(kind=8)::Rinv ! the inverse bond length/distance (for derivatives)
integer::atom  ! the current atom index
 
!
real(kind=8)::dxi_act(3,natoms) ! the calculated Xi derivative
real(kind=8)::d2xi_act(3,natoms,3,natoms)  ! the calculated Xi hessian
real(kind=8)::fs, fs2, log_fs, coeff1, coeff2, dhams ! for bias potential..
pi=3.1415926535897932384d0

!
!     calculate actual value of xi and compare it with
!     desired value
!
act_coord(1,:)=centroid(1,:)
act_coord(2,:)=centroid(2,:)
act_coord(3,:)=centroid(3,:)
call xyz_2int(act_coord,int_real,natoms)

grad_bias = 0
!
!      call routine to calculate Xi value and its derivatives
! 
if (mode .eq. 0) then
   call calc_xi(act_coord,xi_ideal,xi_real,dxi_act,d2xi_act,1) 
!
!      if a constrained trajectory is sampled, calculate directly the relative xi value 
!      and go back to verlet without applying the bias
!
else if (mode .eq. 1) then
      call calc_xi(act_coord,xi_ideal,xi_real,dxi_act,d2xi_act,2)
!
!      directly return to main routine!
!
   return
!
!      if a child trajectory is sampled, calculate usual Xi values and return directly
!      to main routine
!
else if (mode .eq. 3) then
   call calc_xi(act_coord,xi_ideal,xi_real,dxi_act,d2xi_act,1)
   return
end if
!
!      apply umbrella potential to gradient
!
delta=(xi_real-xi_ideal)
v_add=0.5*k_force*delta*delta 
do i=1,nbeads 
   grad_xyz(:,:,i)=grad_xyz(:,:,i)+k_force*delta*dxi_act
end do
!
!      apply bias potential (with hessian) to gradient (why???)
!      ---> temporarily deactivated!
fs2=0.0d0
do i = 1, 3
   do j = 1, Natoms
        fs2 = fs2 + dxi_act(i,j) * dxi_act(i,j) / mass(j)
   end do
end do
coeff1 = 2.0d0 * pi * beta
fs2 = fs2 / coeff1
fs = sqrt(fs2)
log_fs = log(fs)
coeff2 = -1.0d0 / beta

! Add bias term to potential
V_add = V_add + coeff2 * log_fs
! Add bias term to forces
!   if (k_force .ne. 0.d0) then
do i = 1, 3
   do j = 1, Natoms
      dhams = 0.0d0
      do i2 = 1, 3
         do j2 = 1, Natoms
            dhams = dhams + d2xi_act(i2,j2,i,j) * dxi_act(i2,j2) / mass(j2)
         end do
      end do
      dhams = dhams * coeff2 / (coeff1 * fs2)
      do k=1,nbeads
         grad_xyz(i,j,k) = grad_xyz(i,j,k) + dhams
      end do
   end do
end do

return
end subroutine umbrella
