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
!     subroutine nhc: applies the Nose-Hoover chain (NHC) thermostat to a 
!       half step/full step in the Verlet MD integration routine 
!
!     part of EVB
!
subroutine nhc(dt,centroid)
use general
use evb_mod 
implicit none 
integer::i,j,k
integer::nc,ns
real(kind=8)::scale
real(kind=8)::rdum
real(kind=8)::dt,dt2,dt4,dt8,dts,dtc,ekt
real(kind=8)::kbt,expterm
real(kind=8)::eksum
real(kind=8)::centroid(3,natoms)
real(kind=8)::w(3)
!
!     Precalculate frequently used variables
!
ekt=1.380649E-23/4.3597447E-18*kelvin
nc = 5
ns = 3
dtc = dt / dble(nc)
w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
w(2) = 1.0d0 - 2.0d0*w(1)
w(3) = w(1)
scale = 1.0d0
eksum=0.d0

!
!     first, calculate the kinetic energy by applying the centroid-momenta approximation
!     For stick_coeff program: exclude the atoms of the gas phase molecule!
!
if (use_stick_coeff) then
   do i=1,natoms-2
      if (at_move(i)) then
         do j=1,nbeads
            eksum=eksum+dot_product(p_i(:,i,j),p_i(:,i,j))/(2d0*mass(i))/nbeads/nbeads
         end do
      end if
   end do
else
   do i=1,natoms
      if (at_move(i)) then
         do j=1,nbeads
            eksum=eksum+dot_product(p_i(:,i,j),p_i(:,i,j))/(2d0*mass(i))/nbeads/nbeads
         end do
      end if
   end do
end if
!
!     Calculate the actual kinetic energy via the virial theorem
!  
!do i=1,natoms
!   do j=1,nbeads
!      j_lower=j-1
!      if (j_lower .lt. 1) j_lower=nbeads
!      j_upper=j+1
!      if (j_upper .gt. nbeads) j_upper=1
!
!      eksum=eksum+1d0/(2d0*nbeads)*dot_product(q_i(:,i,j)-centroid(:,i),derivs(:,i,j)/nbeads)!+&
!                 !  & mass(i)/beta_n**2/nbeads*((2*q_i(:,i,j)-q_i(:,i,j_upper)-q_i(:,i,j_lower))))
!
!   end do
!end do
!
!     Apply the Nose-Hoover-chain equations
!
do i = 1, nc
   do j = 1, ns
      dts = w(j) * dtc
      dt2 = 0.5d0 * dts
      dt4 = 0.25d0 * dts
      dt8 = 0.125d0 * dts
      gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
      vnh(4) = vnh(4) + gnh(4)*dt4
      gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
      expterm = exp(-vnh(4)*dt8)
      vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
      gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
      expterm = exp(-vnh(3)*dt8)
      vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
      gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
      expterm = exp(-vnh(2)*dt8)
      vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
      expterm = exp(-vnh(1)*dt2)
      scale = scale * expterm
      eksum = eksum * expterm * expterm
      gnh(1) = (2.0d0*eksum-dble(nfree)*ekt) / qnh(1)
      expterm = exp(-vnh(2)*dt8)
      vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
      gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
      expterm = exp(-vnh(3)*dt8)
      vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
      gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
      expterm = exp(-vnh(4)*dt8)
      vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
      gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
      vnh(4) = vnh(4) + gnh(4)*dt4
   end do
end do
!
!    Rescale the momenta
!    For stick_coeff: exclude colliding molecule
!
if (use_stick_coeff) then
   do i = 1, natoms-2
      if (at_move(i)) then
         do j=1,nbeads
            do k = 1, 3
               p_i(k,i,j) = scale * p_i(k,i,j)
            end do
         end do
      else
         do j=1,nbeads
            do k = 1, 3
               p_i(k,i,j) = 0.d0
            end do
         end do
      end if
   end do
else
   do i = 1, natoms
      if (at_move(i)) then
         do j=1,nbeads
            do k = 1, 3
               p_i(k,i,j) = scale * p_i(k,i,j)
            end do
         end do
      else
         do j=1,nbeads
            do k = 1, 3
               p_i(k,i,j) = 0.d0
            end do
         end do
      end if
   end do
end if
return
end subroutine nhc 

