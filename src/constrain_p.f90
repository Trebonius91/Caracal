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
!     subroutine constrain_p: constrain the momentum to the 
!       reaction coordinate, to ensure that the time derivative 
!       of the dividing surface is zero
!    
!     part of EVB
!
subroutine constrain_p(dxi_act)
use general
use evb_mod
implicit none
integer::i,j,k
integer::rank  ! the actual MPI processor ID
real(kind=8)::dt  ! the timestep in atomic units
real(kind=8)::xi_ideal ! the xi value of the PMF max, set as constraint
real(kind=8)::xi_real ! the real xi value, here =xi_ideal
real(kind=8)::dxi_act(3,natoms)  ! the derivative of xi
real(kind=8)::derivs(3,natoms) 
real(kind=8)::coeff1,coeff2,lambdas

coeff1 = 0.0d0
do i = 1, 3
   do j = 1, Natoms
      do k=1,nbeads
         coeff1 = coeff1 + dxi_act(i,j) * p_i(i,j,k) / mass(j)
      end do
   end do
end do

coeff2 = 0.0d0
do i = 1, 3
   do j = 1, Natoms
      coeff2 = coeff2 + dxi_act(i,j) * dxi_act(i,j) / mass(j)
   end do
end do

lambdas = -coeff1 / coeff2 /nbeads 
do i = 1, 3
   do j = 1, Natoms
      do k=1,nbeads
         p_i(i,j,k) = p_i(i,j,k) + lambdas * dxi_act(i,j)
      end do
   end do
end do


return
end subroutine constrain_p
