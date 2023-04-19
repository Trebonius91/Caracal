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
!     subroutine version: checks the name of a file to be opened; highest 
!         current version is returned 
!
!     part of EVB
!
subroutine version (filename1,status)
use general
implicit none
integer i,leng1,trimtext
integer thousand,hundred
integer tens,ones
logical exist,noversion
character(len=1)::digit(0:9)
character(len=3)::status
character(len=120)::filename1
character(len=120)::oldfile
character(len=120)::newfile

digit=(/'0','1','2','3','4','5','6','7','8','9'/)
!
!     process the filename and status variables
!
call lowcase (status)
leng1 = trimtext (filename1)
!
!     no change is needed if the file doesn't exist
!
exist = .false.
if (leng1 .ne. 0)  inquire (file=filename1(1:leng1),exist=exist)
if (.not. exist)  return
!
!     set initial values for the current and next versions
!
newfile = filename1
oldfile = filename1
!
!     append an artificial version number to the filename;
!     currently handles up to 10000 versions of a file
!
if (.not. noversion) then
   i = 1
   do while (exist)
      i = i + 1
      oldfile = newfile
      thousand = i / 1000
      hundred = (i - 1000*thousand) / 100
      tens = (i - 1000*thousand - 100*hundred) / 10
      ones = i - 1000*thousand - 100*hundred - 10*tens
      if (thousand .ne. 0) then
         newfile = filename1(1:leng1)//'_'//digit(thousand) &
              &    //digit(hundred)//digit(tens)//digit(ones)
      else if (hundred .ne. 0) then
         newfile = filename1(1:leng1)//'_'//digit(hundred) &
              &    //digit(tens)//digit(ones)
      else if (tens .ne. 0) then
         newfile = filename1(1:leng1)//'_'//digit(tens)//digit(ones)
      else
         newfile = filename1(1:leng1)//'_'//digit(ones)
      end if
      inquire (file=newfile,exist=exist)
   end do
end if
!
!     set the file name based on the requested status
!
if (status .eq. 'old') then
   filename1 = oldfile
else if (status .eq. 'new') then
   filename1 = newfile
   inquire (file=filename1,exist=exist)
   if (exist) then
      call nextarg (filename1,exist)
      if (exist) then
         inquire (file=filename1,exist=exist)
      else
         exist = .true.
      end if
      do while (exist)
         write (iout,'(/," Enter File Name for Coordinate Output :  ",$)')
         read (input,'(a120)')  filename1
         inquire (file=filename1,exist=exist)
      end do
   end if
end if

return
end subroutine version
