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
!     subroutine cutline: cuts the at blanks and tabstops and returns 
!     all floats and strings in order of occurence
!
!     part of QMDFF
!

subroutine cutline(line,floats,strings)
implicit none
real(kind=8)::floats(*),num
character(len=128)::line,str,stmp
character(len=80)::strings(3)
character(len=1)::digit
integer::i,ty,cs,cf

stmp=''
cs=1
cf=1
strings(:)=''
do i=1,len(trim(line))
   digit=line(i:i)
!
!     should exclude tabstops and blanks, 9 is ascii code for tab
!
   if (digit.ne.' '.and.digit.ne.char(9)) then
      stmp=trim(stmp)//trim(digit)
   else if (stmp.ne.'') then
!
!     get type of string, 0=number, 1=character
!
      call checktype(stmp,num,str,ty)   
      if (ty.eq.0) then
         floats(cf)=num
         cf=cf+1
      else if (ty.eq.1) then
         strings(cs)=str
         cs=cs+1
      else
         write(10,*)'Problem in checktype, must abort'
         exit
      end if
      stmp=''
   end if
!
!     special case: end of line
!
   if (i.eq.len(trim(line))) then 
      call checktype(stmp,num,str,ty)
      if (ty.eq.0) then
         floats(cf)=num
         cf=cf+1
      else if (ty.eq.1) then
         strings(cs)=str
         cs=cs+1
      else
         write(10,*)'Problem in checktype, must abort'
         exit
      end if
   stmp=''
   end if
end do
return
end subroutine cutline 
