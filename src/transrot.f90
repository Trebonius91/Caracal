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
!     subroutine transrot: removes and translational or rotational kinetic 
!        energy of the overall systems center of mass during a MD simulation
!        ---> Code mainly taken from TINKER
!
!     part of EVB
!
subroutine transrot 
use general 
use evb_mod
implicit none 
integer::i,j,k
real(kind=8)::weigh,totmass,eps
real(kind=8)::xx,yy,zz,xy,xz,yz
real(kind=8)::xtot,ytot,ztot
real(kind=8)::xdel,ydel,zdel
real(kind=8)::mang(3),vang(3)
real(kind=8)::vtot(3),tensor(3,3)
real(kind=8),allocatable::vel(:,:,:)
real(kind=8)::erot,etrans
real(kind=8)::rubbish(3)

allocate(vel(3,natoms,nbeads))
!
!     Calculate velocity from momentum
!
do i=1,natoms
   vel(:,i,:)=p_i(:,i,:)/mass(i)
end do

!
!     Initialize variables
!
totmass = 0.0d0
do j = 1, 3
   vtot(j) = 0.0d0
end do

!
!     Compute linear velocity of the system center of mass
!
do i = 1, natoms
   weigh = mass(i)
   do k=1,nbeads
      totmass = totmass + weigh
      do j = 1, 3
         vtot(j) = vtot(j) + vel(j,i,k)*weigh
      end do
   end do
end do
totmass=totmass*nbeads
!
!     Compute translational kinetic energy of overall system
!
etrans = 0.0d0
do j = 1, 3
   vtot(j) = vtot(j) / totmass
   etrans = etrans + vtot(j)**2
end do
etrans = 0.5d0 * etrans * totmass

!
!     Find center of mass coordinates for overall system
!

xtot = 0.0d0
ytot = 0.0d0
ztot = 0.0d0
do i = 1, natoms
   weigh = mass(i)
   do j=1,nbeads
      xtot = xtot + q_i(1,i,j)*weigh
      ytot = ytot + q_i(2,i,j)*weigh
      ztot = ztot + q_i(3,i,j)*weigh
   end do
end do

xtot = xtot / totmass
ytot = ytot / totmass
ztot = ztot / totmass

!
!     Compute angular momentum for overall system
!
do j = 1, 3
   mang(j) = 0.0d0
end do

do i = 1, natoms
   weigh = mass(i)
   do k=1,nbeads
      mang(1) = mang(1) + (q_i(2,i,k)*vel(3,i,k)-q_i(3,i,k)*vel(2,i,k))*weigh
      mang(2) = mang(2) + (q_i(3,i,k)*vel(1,i,k)-q_i(1,i,k)*vel(3,i,k))*weigh
      mang(3) = mang(3) + (q_i(1,i,k)*vel(2,i,k)-q_i(2,i,k)*vel(1,i,k))*weigh
   end do
end do


mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass
mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass
mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass

!
!     calculate the moment of inertia tensor
!
xx = 0.0d0
xy = 0.0d0
xz = 0.0d0
yy = 0.0d0
yz = 0.0d0
zz = 0.0d0

do i = 1, natoms
   weigh = mass(i)
   do k=1,nbeads
      xdel = q_i(1,i,k) - xtot
      ydel = q_i(2,i,k) - ytot
      zdel = q_i(3,i,k) - ztot
      xx = xx + xdel*xdel*weigh
      xy = xy + xdel*ydel*weigh
      xz = xz + xdel*zdel*weigh
      yy = yy + ydel*ydel*weigh
      yz = yz + ydel*zdel*weigh
      zz = zz + zdel*zdel*weigh
   end do
end do
tensor(1,1) = yy + zz
tensor(2,1) = -xy
tensor(3,1) = -xz
tensor(1,2) = -xy
tensor(2,2) = xx + zz
tensor(3,2) = -yz
tensor(1,3) = -xz
tensor(2,3) = -yz
tensor(3,3) = xx + yy
!
!     avoid bad behavior (singularity) for diatomic systems
!

if (natoms .le. 2) then
   eps = 0.000001d0
   tensor(1,1) = tensor(1,1) + eps
   tensor(2,2) = tensor(2,2) + eps
   tensor(3,3) = tensor(3,3) + eps
end if
! 
!     invert the moment of inertia tensor
!
call invert (3,tensor)

!
!     compute angular velocity and rotational kinetic energy
!
erot = 0.0d0
do i = 1, 3
   vang(i) = 0.0d0
   do j = 1, 3
      vang(i) = vang(i) + tensor(i,j)*mang(j)
   end do
   erot = erot + vang(i)*mang(i)
end do
erot = 0.5d0 * erot 

!
!     eliminate any translation of the overall system
!
do i = 1, natoms
   do k=1,nbeads
      do j = 1, 3
         vel(j,i,k) = vel(j,i,k) - vtot(j)
      end do
   end do
end do

!
!     eliminate any rotation of the overall system
!
do i = 1, natoms
   do k=1,nbeads
      xdel = q_i(1,i,k) - xtot
      ydel = q_i(2,i,k) - ytot
      zdel = q_i(3,i,k) - ztot
      vel(1,i,k) = vel(1,i,k) - vang(2)*zdel + vang(3)*ydel
      vel(2,i,k) = vel(2,i,k) - vang(3)*xdel + vang(1)*zdel
      vel(3,i,k) = vel(3,i,k) - vang(1)*ydel + vang(2)*xdel
   end do
end do
!
!     Overwrite momentum with corrected velocity
!

do i=1,natoms
   p_i(:,i,:)=vel(:,i,:)*mass(i)
end do


!
!     set momentum to zero for fixed atoms
!
do i=1,natoms
   if (.not. at_move(i)) then
      p_i(:,i,:)=0.d0
   end if
end do




end subroutine transrot
