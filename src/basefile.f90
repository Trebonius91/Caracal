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
!     subroutine basefile: extracts directory name and base filename 
!            from input filename
!
!     part of EVB
!
subroutine basefile (string)
use general
implicit none
integer::i,k,trimtext
character(len=1)::letter
character(len=120)::string
!
!     store the input filename and find its full length
!
filename = string
leng = trimtext (string)
!
!     count the number of characters prior to any extension
!
k = leng
do i = 1, leng
   letter = string(i:i)
   if (letter .eq. '/')  k = leng
   if (ichar(letter) .eq. backslash)  k = leng
   if (letter .eq. ']')  k = leng
   if (letter .eq. ':')  k = leng
   if (letter .eq. '~')  k = leng
   if (letter .eq. '.')  k = i - 1
end do
leng = min(leng,k)
!
!     find the length of any directory name prefix
!
k = 0
do i = leng, 1, -1
   letter = string(i:i)
   if (letter .eq. '/')  k = i
   if (ichar(letter) .eq. backslash)  k = i
   if (letter .eq. ']')  k = i
   if (letter .eq. ':')  k = i
   if (letter .eq. '~')  k = i
   if (k .ne. 0)  exit
end do
ldir = k
!
!     read and store the keywords from the keyfile
!
call getkey(0)
!
!     get the information level and output style
!
!call control  -->  not needed here?

return
end subroutine basefile
