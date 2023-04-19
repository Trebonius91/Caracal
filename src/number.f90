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
!     function "number" converts a text numeral into an integer value;
!     the input string must contain only numeric characters
!
!
!     Based on:
!     TINKER molecular modeling package
!     COPYRIGHT (C)  1990  by  Jay William Ponder   
!     All Rights Reserved 
!
!     part of others
!
function number (string)
use general
implicit none
integer i,j,number
integer first,last,trimtext
integer digit,place(10)
character*1 letter
character*(*) string
place=(/ 1, 10, 100, 1000, 10000, 100000, 1000000, &
         &    10000000, 100000000, 1000000000 /)
!
!     initialize the integer value of number to zero
!
number = 0
!
!     get the first and last nonblank characters
!
last = trimtext (string)
if (last .gt. 10) then
   write (iout,*) "function number: Text String is Too Long"
end if
first = 1
do i = 1, last
   letter = string(i:i)
   if (letter .ne. ' ') then
      first = i
      exit
   end if
end do
!
!     convert the text numeral into an integer number
!
j = 0
do i = last, first, -1
   j = j + 1
   letter = string(i:i)
   if (letter .eq. '0') then
      digit = 0
   else if (letter .eq. '1') then
      digit = 1
   else if (letter .eq. '2') then
      digit = 2
   else if (letter .eq. '3') then
      digit = 3
   else if (letter .eq. '4') then
      digit = 4
   else if (letter .eq. '5') then
      digit = 5
   else if (letter .eq. '6') then
      digit = 6
   else if (letter .eq. '7') then
      digit = 7
   else if (letter .eq. '8') then
      digit = 8
   else if (letter .eq. '9') then
      digit = 9
   else
      number = 0
      exit
   end if
   number = number + digit * place(j)
end do
return
end function number
