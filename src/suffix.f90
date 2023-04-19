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
!     subroutine suffix: check filename for presence of an extension and 
!           appends such an extension if it couldn´t be found
!
!     part of EVB
!
subroutine suffix (filename1,extension,status)
use general
implicit none
integer::i,leng1,lext1
integer::last,trimtext
logical::exist
character(len=1)::letter
character(len=3)::status
character(len=*)::filename1
character(len=*)::extension
!
!     get the length of the current filename
!
leng1 = trimtext (filename1)
lext1 = trimtext (extension)
!
!     check for an extension on the current filename
!
last = leng1
do i = 1, leng1
   letter = filename1(i:i)
   if (letter .eq. '/')  last = leng1
   if (ichar(letter) .eq. backslash)  last = leng1
   if (letter .eq. ']')  last = leng1
   if (letter .eq. ':')  last = leng1
   if (letter .eq. '~')  last = leng1
   if (letter .eq. '.')  last = i - 1
end do
!
!     append an extension or version as appropriate
!
if (last .eq. leng1) then
   exist = .false.
   if (leng1 .ne. 0) then
      inquire (file=filename1(1:leng1),exist=exist)
   end if
   if (.not. exist) then
      filename = filename1(1:leng1)//'.'//extension(1:lext1)
      call version (filename1,status)
   end if
else if (status .eq. 'new') then
   call version (filename1,status)
end if

return
end subroutine suffix

