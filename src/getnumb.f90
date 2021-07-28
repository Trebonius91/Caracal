!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   EVB-QMDFF - RPMD molecular dynamics and rate constant calculations on
!               black-box generated potential energy surfaces
!
!   Copyright (c) 2021 by Julien Steffen (steffen@pctc.uni-kiel.de)
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
!     subroutine "getnumb" searches an input string from left to right for an
!     integer and puts the numeric value in "number"; returns zero
!     with "next" unchanged if no integer value is found
!
!     part of others
!
subroutine getnumb (string,number,next)
use general
implicit none
integer::i,j,length
integer::number,digit
integer::next,trimtext
integer::first,last,code
integer::initial,final
integer::place(10)
logical::negate,numeral
character(len=1)::letter
character(len=*)::string
place=(/ 1, 10, 100, 1000, 10000, 100000, 1000000, &
       &     10000000, 100000000, 1000000000 /)

!
!     initialize number and get the input text string length
!
number = 0
negate = .false.
numeral = .false.
length = trimtext(string(next:))
!
!     search the string for the first run of numeric characters
!
first = next
last = 0
initial = next
final = next + length - 1
do i = initial, final
   letter = string(i:i)
   code = ichar(letter)
   if (letter.ge.'0' .and. letter.le.'9') then
      if (.not. numeral) then
         numeral = .true.
         first = i
      end if
      if (i .eq. final) then
         last = final
         next = i + 1
      end if
   else if (code.eq.minus .and. .not.negate) then
      negate = .true.
   else if (numeral) then
      if (code.eq.space .or. code.eq.tab .or. &
         &    code.eq.comma .or. code.eq.semicolon .or. &
         &    code.eq.colon .or. code.eq.underbar) then
         last = i - 1
         next = i
      else
         numeral = .false.
      end if
      exit
   else if (negate) then
      numeral = .false.
      exit
   else if (code.ne.space .and. code.ne.tab) then
      numeral = .false.
      exit
   end if
end do
!
!     trim the actual number if it is too big to return
!
if (.not. numeral)  next = initial
last = min(last,first+9)
!
!     convert the text numeral into an integer number
!
j = 0
do i = last, first, -1
   j = j + 1
   if (string(i:i) .eq. '0') then
      digit = 0
   else if (string(i:i) .eq. '1') then
      digit = 1
   else if (string(i:i) .eq. '2') then
      digit = 2
   else if (string(i:i) .eq. '3') then
      digit = 3
   else if (string(i:i) .eq. '4') then
      digit = 4
   else if (string(i:i) .eq. '5') then
      digit = 5
   else if (string(i:i) .eq. '6') then
      digit = 6
   else if (string(i:i) .eq. '7') then
      digit = 7
   else if (string(i:i) .eq. '8') then
      digit = 8
   else if (string(i:i) .eq. '9') then
      digit = 9
   end if
   number = number + digit * place(j)
end do
if (negate)  number = -number
return
end subroutine getnumb 
