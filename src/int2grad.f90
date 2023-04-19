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
!     subroutine int2grad: converts an array of internal gradients to cartesian 
!     coordinates
!
!     part of EVB
!

subroutine int2grad(xyz5,internal,g_q,g_x)
use evb_mod
implicit none

real(kind=8)::xyz5(3,natoms),internal(nat6)
real(kind=8)::B_mat(nat6,3*natoms)
real(kind=8)::Bt(3*natoms,nat6)
real(kind=8)::u(3*natoms,3*natoms)
real(kind=8)::Bu(3*natoms,3*natoms)
real(kind=8)::G(3*natoms,3*natoms),V(3*natoms,3*natoms)
real(kind=8)::xyza(3*natoms),eigen(3*natoms)
real(kind=8)::Vt(3*natoms,3*natoms)
real(kind=8),dimension(3*natoms,3*natoms)::Vlam,V1
real(kind=8)::G_min(3*natoms,3*natoms)
real(kind=8)::G_minB(3*natoms,nat6)
real(kind=8)::G_minBu(nat6,3*natoms)
real(kind=8)::g_q(nat6),g_x(3*natoms),g_q1(3*natoms)
real(kind=8)::ug_x(3*natoms),Bug_x(3*natoms),g_test(3*natoms)
real(kind=8)::coord,hi,lo,shift
real(kind=8)::dist,ang,dihed,oop
real(kind=8)::q_test(nat6)
! only for test reasons
real(kind=8)::g2_x(3,natoms),gu_q(nat6),gd_q(nat6)
real(kind=8)::method2
integer::n


integer::i,j,k,m,l,c,o
integer::nat
nat=natoms
n=nat
do i=1,nat
   do j=1,3
      xyza((i-1)*3+j)=xyz5(j,i)
   end do 
end do
!
!  first, build the Wilson B matrix
!
call calc_wilson(xyz5,internal,B_mat)
if (B_mat(1,1) .ne. B_mat(1,1)) then
write(*,*) "NAN in Wikson!"
end if
!
!  Transpose the Wilson Matrix B (--> Bt)
!
Bt=0
do i=1,nat6
   do j=1,3*nat
      Bt(j,i)=B_mat(i,j)
   end do
end do

!
!   calculate the cartesian gradient
!

g_x=matmul(Bt,g_q)

return
end subroutine int2grad
