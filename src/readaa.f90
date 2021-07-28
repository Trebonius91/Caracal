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
!     function readaa: analyze output file line for Turbomole
!
!     part of QMDFF
!
function readaa(a,istart,iend,iend2)
implicit real(kind=8) (a-h,o-z)
real(kind=8)::readaa
character(len=*)::a

NINE=ICHAR('9')
IZERO=ICHAR('0')
MINUS=ICHAR('-')
IDOT=ICHAR('.')
ND=ICHAR('D')
NE=ICHAR('E')
IBL=ICHAR(' ')

iend=0
iend2=0
idig=0
c1=0
c2=0
one=1.d0
x = 1.d0
nl=len(a)
do j=istart,nl-1
   n=ichar(a(j:j))
   m=ichar(a(j+1:j+1))
   if (n.le.nine.and.n.ge.izero .or.n.eq.idot) goto 20
   if (n.eq.minus.and.(m.le.nine.and.m.ge.izero &
   & .or. m.eq.idot)) goto 20
end do
readaa=0.d0
return
20 continue
iend=j
do i=j,nl
   n=ichar(a(i:i))
   if (n.le.nine.and.n.ge.izero) then
      idig=idig+1
      if (idig.gt.10) goto 60
      c1=c1*10+n-izero
   else if (n.eq.minus.and.i.eq.j) then
      one=-1.d0
   else if (n.eq.idot) then
      goto 40
   else
      goto 60
   end if
end do
40 continue
idig=0
do ii=i+1,nl
   n=ichar(a(ii:ii))
   if(n.le.nine.and.n.ge.izero) then
      idig=idig+1
      if (idig.gt.10) goto 60
      c2=c2*10+n-izero
      x = x /10
   else if (n.eq.minus.and.ii.eq.i) then
      x=-x
   else
      goto 60
   end if
end do
!                                                                               
! put the pieces together                                                       
! 
60 continue
readaa= one * ( c1 + c2 * x)
do 55 j=iend,nl
   n=ichar(a(j:j))
   iend2=j
   if (n.eq.ibl) return
   55 if(n.eq.nd .or. n.eq.ne) goto 57
return
57 c1=0.0d0
one=1.0d0

do 31 i=j+1,nl
   n=ichar(a(i:i))
   iend2=i
   if(n.eq.ibl)goto 70
   if(n.le.nine.and.n.ge.izero) c1=c1*10.0d0+n-izero
   if(n.eq.minus)one=-1.0d0
31 continue
61 continue
70 readaa=readaa*10**(one*c1)

return
end function readaa
