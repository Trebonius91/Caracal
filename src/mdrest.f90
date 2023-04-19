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
!     subroutine mdrest: finds and removes any translational or rotational
!     kinetic energy of the overall system center of mass
!
!
!     Based on:
!     TINKER molecular modeling package
!     COPYRIGHT (C)  1990  by  Jay William Ponder   
!     All Rights Reserved 
!
!     part of EVB
!
subroutine mdrest (istep)
use general
implicit none
integer::i,j,k,istep
real(kind=8)::etrans,erot
real(kind=8)::weigh,totmass1,eps
real(kind=8)::xx,yy,zz,xy,xz,yz
real(kind=8)::xtot,ytot,ztot
real(kind=8)::xdel,ydel,zdel
real(kind=8)::mang(3),vang(3)
real(kind=8)::vtot(3),tensor(3,3)
real(kind=8),allocatable::xcm(:)
real(kind=8),allocatable::ycm(:)
real(kind=8),allocatable::zcm(:)

!
!     check steps between center of mass motion removal
!
if (.not.dorest)  return
if (mod(istep,irest) .ne. 0)  return
!
!     zero out the total mass and overall linear velocity
!
totmass1 = 0.0d0
do j = 1, 3
   vtot(j) = 0.0d0
end do
!
!     compute linear velocity of the system center of mass
!
do i = 1, n
   weigh = mass(i)
   totmass1 = totmass1 + weigh
   do j = 1, 3
      vtot(j) = vtot(j) + v(j,i)*weigh
   end do
end do
!
!     compute translational kinetic energy of overall system
!
etrans = 0.0d0
do j = 1, 3
   vtot(j) = vtot(j) / totmass1
   etrans = etrans + vtot(j)**2
end do
etrans = 0.5d0 * etrans * totmass1 / convert

!
!     find the center of mass coordinates of the overall system
!

xtot = 0.0d0
ytot = 0.0d0
ztot = 0.0d0
   do i = 1, n
      weigh = mass(i)
      xtot = xtot + x(i)*weigh
      ytot = ytot + y(i)*weigh
      ztot = ztot + z(i)*weigh
   end do
xtot = xtot / totmass1
ytot = ytot / totmass1
ztot = ztot / totmass1
!
!     compute the angular momentum of the overall system
!
do j = 1, 3
   mang(j) = 0.0d0
end do
do i = 1, n
   weigh = mass(i)
   mang(1) = mang(1) + (y(i)*v(3,i)-z(i)*v(2,i))*weigh
   mang(2) = mang(2) + (z(i)*v(1,i)-x(i)*v(3,i))*weigh
   mang(3) = mang(3) + (x(i)*v(2,i)-y(i)*v(1,i))*weigh
end do
mang(1) = mang(1) - (ytot*vtot(3)-ztot*vtot(2))*totmass1
mang(2) = mang(2) - (ztot*vtot(1)-xtot*vtot(3))*totmass1
mang(3) = mang(3) - (xtot*vtot(2)-ytot*vtot(1))*totmass1
!
!     calculate the moment of inertia tensor
!
xx = 0.0d0
xy = 0.0d0
xz = 0.0d0
yy = 0.0d0
yz = 0.0d0
zz = 0.0d0
do i = 1, n
   weigh = mass(i)
   xdel = x(i) - xtot
   ydel = y(i) - ytot
   zdel = z(i) - ztot
   xx = xx + xdel*xdel*weigh
   xy = xy + xdel*ydel*weigh
   xz = xz + xdel*zdel*weigh
   yy = yy + ydel*ydel*weigh
   yz = yz + ydel*zdel*weigh
   zz = zz + zdel*zdel*weigh
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
!     fix to avoid singularity for one- or two-body systems
!
if (n .le. 2) then
   eps = 0.000001d0
   tensor(1,1) = tensor(1,1) + eps
   tensor(2,2) = tensor(2,2) + eps
   tensor(3,3) = tensor(3,3) + eps
end if
!
!     diagonalize the moment of inertia tensor
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
erot = 0.5d0 * erot / convert

!
!     eliminate any translation of the overall system
!
do i = 1, n
   do j = 1, 3
      v(j,i) = v(j,i) - vtot(j)
   end do
end do
!
!     eliminate any rotation about the system center of mass
!

do i = 1, n
   xdel = x(i) - xtot
   ydel = y(i) - ytot
   zdel = z(i) - ztot
   v(1,i) = v(1,i) - vang(2)*zdel + vang(3)*ydel
   v(2,i) = v(2,i) - vang(3)*xdel + vang(1)*zdel
   v(3,i) = v(3,i) - vang(1)*ydel + vang(2)*xdel
end do

return
end subroutine mdrest 
