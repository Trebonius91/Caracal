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
!     "getstring" searches for a quoted text string within an input
!     character string; the region between the first and second
!     double quote is returned as the "text"; if the actual text is
!     too long, only the first part is returned
!
!     Based on:
!     TINKER molecular modeling package
!     COPYRIGHT (C)  1990  by  Jay William Ponder   
!     All Rights Reserved 
!
!
!     part of others
!
subroutine getstring (string,text,next)
use general
implicit none
integer::i,j,k,m
integer::len,length
integer::size,next
integer::code,extent
integer::initial,final
integer::first,last
integer::maxascii
character(len=*)::string
character(len=*)::text

!
!     get the length of input string and output text
!
length = len(string(next:))
size = len(text)
!
!     convert first two non-ascii regions to delimiting quotes
!
maxascii = 126
initial = next
final = next + length - 1
do i = initial, final
   code = ichar(string(i:i))
   if (code .gt. maxascii) then
      string(i:i) = ' '
      do j = i+1, final
         code = ichar(string(j:j))
         if (code .le. maxascii) then
            string(j-1:j-1) = '"'
            do k = j+1, final
               code = ichar(string(k:k))
               if (code .gt. maxascii) then
                  string(k:k) = '"'
                  do m = k+1, final
                     code = ichar(string(m:m))
                     if (code .gt. maxascii) then
                        string(m:m) = ' '
                     end if
                  end do
                  goto 10
               end if
            end do
         end if
      end do
   end if
end do
10 continue
!
!     search the string for quoted region of text characters
!
first = next
last = 0
do i = initial, final
   code = ichar(string(i:i))
   if (code .eq. quote) then
      first = i + 1
      do j = first, final
         code = ichar(string(j:j))
         if (code .eq. quote) then
            last = j - 1
            next = j + 1
            goto 20
         end if
      end do
   end if
end do
20 continue
!
!     trim the actual word if it is too long to return
!
extent = last - first + 1
final = first + size - 1
if (extent .gt. size)  last = final
!
!     transfer the text into the return string
!
j = 0
do i = first, last
   j = j + 1
   text(j:j) = string(i:i)
end do
do i = last+1, final
   j = j + 1
   text(j:j) = ' '
end do
return
end subroutine getstring
