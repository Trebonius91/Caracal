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
!     subroutine int_dist:  calculates the distance between two points 
!     in internal coordinate space. For angles and dihedrals the usual pi/2pi shifts 
!     are considered in order to avoid artificial large jumps
!
!     part of EVB
!
subroutine int_diff(internal1,internal2,internal3)
use general
use evb_mod
use qmdff
implicit none 
real(kind=8),intent(in)::internal1(nat6),internal2(nat6)
real(kind=8),intent(out)::internal3(nat6)
integer::i,j,k  ! loop indices
real(kind=8)::dist1,dist2   ! temporal angle distances

!
!     calculate differences for each kind of coordinates 
!
do i=1,nat6
! 
!     bondlength: usual difference
!
   if (coord_def(i,1) .eq. 1) then
      internal3(i)=internal1(i)-internal2(i)
!
!     angles: consider also the Pi periodicity!
!   
   else if (coord_def(i,1) .eq. 2) then
      if (internal1(i) .gt. internal2(i)) then
         dist1=internal1(i)-internal2(i)
         dist2=internal2(i)-(internal1(i)-pi)
         if (dist1 .gt. dist2) then
            internal3(i)=-dist2
     !       write(*,*) "case!"
         else 
            internal3(i)=dist1
         end if
      else   
         dist1=internal1(i)-internal2(i)
         dist2=internal1(i)-(internal2(i)-pi)
         if (dist1 .gt. dist2) then
            internal3(i)=dist2
     !       write(*,*) "case!"
         else
            internal3(i)=dist1
         end if
      end if
!
!     dihedral angles 
!
   else if (coord_def(i,1) .eq. 3) then
      internal3(i)=internal1(i)-internal2(i)
!      if (internal3(i) .gt. pi/2d0) then
!         write(*,*) "big!",internal3(i)
!      end if
!
!     out of plane angles: no periodicity needed!
!
   else if (coord_def(i,1) .eq. 4) then
      internal3(i)=internal1(i)-internal2(i)
   end if
!   write(*,*) coord_def(i,1)
!   write(*,*) internal1(i),internal2(i)

end do
!do i=1,nat6
!   write(*,*) "output",internal3(i)

!end do
!stop "HIphpui"
return
end subroutine int_diff
