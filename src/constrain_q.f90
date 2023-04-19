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
!     subroutine constrain_q: constrain the position and the momentum
!       to the dividing surface, using the SHAKE/RATTLE algorithm
!    
!     part of EVB
!
subroutine constrain_q(centroid,xi_ideal,dxi_act,const_good,dt)
use general
use evb_mod
implicit none
integer::i,j,k
integer::rank  ! the actual MPI processor ID
integer::iter,maxiters
integer::const_good  ! if the constraint calculation was successful
real(kind=8)::xi_ideal ! the xi value of the PMF max, set as constraint
real(kind=8)::xi_real ! the real xi value, here =xi_ideal
real(kind=8)::centroid(3,natoms)
real(kind=8)::qtemp(3,natoms)
real(kind=8)::derivs(3,natoms) 
real(kind=8)::coeff,sigma,dsigma,dx,mult,dt
real(kind=8)::dxi_act(3,natoms)
real(kind=8)::xi_new,dxi_new(3,natoms),d2xi_new(3,natoms,3,natoms)
!
!     the Lagrange multiplier imposing the constraint
!
mult=0.d0
!
!     a temporary save of the structure 
!
qtemp=0.d0
!
!     maximum number of iterations
!
maxiters=200

do iter=1,maxiters

   coeff=mult*dt*dt/nbeads
   do i=1,3
      do j=1,natoms
         qtemp(i,j)=centroid(i,j)+coeff*dxi_act(i,j)/mass(j)
      end do
   end do

!
!     get the reaction coordinate and its derivatives for the 
!     actual structure 
!
   call calc_xi(qtemp,xi_ideal,xi_new,dxi_new,d2xi_new,2) 
   sigma=xi_new
   dsigma=0.d0
!
!     shift, calculated from gradient
!    
   do i=1,3
      do j=1,natoms
         dsigma=dsigma+dxi_new(i,j)*dt*dt*dxi_act(i,j)/(mass(j)*nbeads)
      end do
   end do 
!
!     check convergence critereon
!
   dx=sigma/dsigma
   mult=mult-dx
   if (dabs(dx) .lt. 1.0E-8 .or. dabs(sigma) .lt. 1.0E-10) exit
!
!     If too many iterations were needed: signal error and go back to main routine
!
   if (iter .eq. maxiters) then
      const_good=1
      return
   end if
end do
do i = 1, 3
   do j = 1, Natoms
      do k=1,nbeads
         q_i(i,j,k) = q_i(i,j,k) + coeff / mass(j) * dxi_act(i,j)
         p_i(i,j,k) = p_i(i,j,k) + mult * dt/nbeads * dxi_act(i,j)
      end do
   end do
end do


return
end subroutine constrain_q
