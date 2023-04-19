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
!    subroutine isort: Sort an array of integer numbers
!        in ascending order
!
!    part of QMDFF
! 
subroutine isort(lab,ew)
implicit none
integer::ew(*),pp,lab,i,k,ii,j

do ii = 2, lab
   i = ii - 1
   k = i
   pp= ew(i)
   do j = ii, lab
      if (ew(j) .gt. pp) cycle
      k = j
      pp= ew(j)
   end do
   if (k .eq. i) cycle
   ew(k) = ew(i)
   ew(i) = pp
end do

return
end subroutine isort

!
!     subroutine isort2: takes an input list of integers and sorts it into
!     ascending order using the Heapsort algorithm, duplicate
!     values are removed from the final sorted list
!
subroutine isort2 (n,list)
implicit none
integer::i,j,k,n
integer::index
integer::lists
integer::list(*)
!
!     perform the heapsort of the input list
!
k = n/2 + 1
index = n
do while (n .gt. 1)
   if (k .gt. 1) then
      k = k - 1
      lists = list(k)
   else
      lists = list(index)
      list(index) = list(1)
      index = index - 1
      if (index .le. 1) then
         list(1) = lists
!
!     remove duplicate values from final list
!
         j = 1
         do i = 2, n
            if (list(i-1) .ne. list(i)) then
               j = j + 1
               list(j) = list(i)
            end if
         end do
         if (j .lt. n)  n = j
         return
      end if
   end if
   i = k
   j = k + k
   do while (j .le. index)
      if (j .lt. index) then
         if (list(j) .lt. list(j+1))  j = j + 1
      end if
      if (lists .lt. list(j)) then
         list(i) = list(j)
         i = j
         j = j + j
      else
         j = index + 1
      end if
   end do
   list(i) = lists
end do
return
end subroutine isort2

