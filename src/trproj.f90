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
!      subroutine trproj: drives projection of hessian out of the
!        space of translational and rotational motions:
!        first get xyz c.m.; second get transl. and rot. projection matrix
!
!      part of QMDFF
!
subroutine trproj(natoms,nat3,xyz,hess,ldebug)

logical,intent(in)::ldebug
integer,intent(in)::natoms,nat3
real(kind=8),dimension(3,natoms)::xyz
real(kind=8),dimension(nat3*(nat3+1)/2)::hess
integer::i
real(kind=8)::xm,ym,zm
real(kind=8),dimension(3,natoms)::xyzucm

xyzucm(:,:) = xyz(:,:)

xm = 0.0d0
ym = 0.0d0
zm = 0.0d0
!
!     move cartesian coordinates
!
do i=1,natoms
   xm = xm + xyzucm(1,i)
   ym = ym + xyzucm(2,i)
   zm = zm + xyzucm(3,i)
end do

xm = xm/natoms
ym = ym/natoms
zm = zm/natoms


do i=1,natoms
   xyzucm(1,i) = xyzucm(1,i) - xm
   xyzucm(2,i) = xyzucm(2,i) - ym
   xyzucm(3,i) = xyzucm(3,i) - zm
end do
!
!     get translational and rotational projection matrix
!
call gtrprojm(natoms,nat3,xyzucm,hess,ldebug)

return
end subroutine trproj
