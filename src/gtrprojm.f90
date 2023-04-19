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
!     subroutine gtrprojm: calculating the translational-rotational 
!              projection matrix (coordinate transformation)
!
!     part of QMDFF
!

subroutine gtrprojm(natoms,nat3,xyzucm,hess,ldebug)

logical, intent(in) :: ldebug
integer, intent(in) :: natoms,nat3
real*8, dimension(3,natoms) :: xyzucm
real*8, dimension(nat3*(nat3+1)/2) :: hess
integer :: i,ii
real*8, dimension(nat3,6) :: fmat 

fmat(:,:) = 0.0d0
!
!     translation vectors
!
do i=1,natoms
   do ii=1,3
     fmat(3*(i-1)+ii,ii) = 1.0d0
   end do
!
!     rotational vectors
!
   fmat(3*(i-1)+1,4)   =  0.0d0
   fmat(3*(i-1)+2,4) = -xyzucm(3,i)
   fmat(3*(i-1)+3,4) =  xyzucm(2,i)

   fmat(3*(i-1)+1,5)   =  xyzucm(3,i)
   fmat(3*(i-1)+2,5) =  0.0d0
   fmat(3*(i-1)+3,5) = -xyzucm(1,i)

   fmat(3*(i-1)+1,6)   =  -xyzucm(2,i) 
   fmat(3*(i-1)+2,6) =   xyzucm(1,i)
   fmat(3*(i-1)+3,6) =   0.0d0
end do

if (ldebug) then
   write(*,'(a)')
   write(*,'(a)') ' Basis vectors before orthonormalization'
   write(*,'(3e22.14)') fmat
end if
!
!     do orthogonalization 
!
call blckmgs(nat3,6,nat3,fmat)
!     
!     do projection
!
call dsyprj(nat3,6,fmat,nat3,hess)

return
end subroutine gtrprojm
