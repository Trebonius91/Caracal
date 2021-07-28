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
!     subroutine getf: Find out if an atom is fixed
!
!     part of QMDFF
!
subroutine getf(line,fix)
implicit none
character(len=*)::line
CHARACTER(len=3)::e
logical::fix
integer::i,j,k,l,n
e='   '
k=1
fix=.false.
do J=1,len(line)
   if (k.gt.2)exit
   N=ICHAR(line(J:J))
!
!     break if space after elem. symbol
!
   if (len_trim(e).ge.1 .and. n.eq.ichar(' ')) exit
!
!     break if tab after elem. symbol
!
   if (len_trim(e).ge.1 .and. n.eq.9) exit 
   if (n.ge.ichar('a') .and. n.le.ichar('z') ) then
      e(k:k)=line(j:j)
      k=k+1
   end if
end do
if (j.eq.len(line))then
   return
else
   do i=j,len(line)
      n=ichar(line(i:i))
      if (n.eq.ichar('f')) fix=.true.
   end do
end if

end subroutine getf
