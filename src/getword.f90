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
!     subroutine getword: search input string for first alphabetic character
!          and return until blank space or separator
!
!
!     Based on:
!     TINKER molecular modeling package
!     COPYRIGHT (C)  1990  by  Jay William Ponder   
!     All Rights Reserved 
!
!     part of EVB
!
subroutine getword (string,word,next)
use general
implicit none
integer::i,j
integer::len,length
integer::size,next
integer::first,last
integer::code,extent
integer::initial,final
character(len=1)::letter
character(len=*)::string
character(len=*)::word
!
!     get the length of input string and output word
!
length = len(string(next:))
size = len(word)
!
!     move through the string one character at a time,
!     searching for the first alphabetic character
!
first = next
last = 0
initial = next
final = next + length - 1
do i = initial, final
   letter = string(i:i)
   if ((letter.ge.'A' .and. letter.le.'Z') .or. &
    &  (letter.ge.'a' .and. letter.le.'z')) then
      first = i
      do j = i+1, final
         code = ichar(string(j:j))
         if (code.eq.space .or. code.eq.tab .or. &
          &  code.eq.comma .or. code.eq.colon .or. &
          &  code.eq.semicolon) then
            last = j - 1
            next = j
            goto 10
         end if
      end do
      last = final
      next = last + 1
   end if
end do
10 continue
!
!     trim the actual word if it is too long to return
!
extent = next - first
final = first + size - 1
if (extent .gt. size)  last = final
!
!     transfer the word into the return string
!
j = 0
do i = first, last
   j = j + 1
   word(j:j) = string(i:i)
end do
do i = next, final
   j = j + 1
   word(j:j) = ' '
end do
!
!     skip over the next character when it is a separator
!
code = ichar(string(next:next))
if (code.eq.tab .or. code.eq.comma .or. &
 &  code.eq.colon .or. code.eq.semicolon) then
    next = next + 1
end if

return
end subroutine getword
