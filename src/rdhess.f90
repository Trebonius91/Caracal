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
!     subroutine rdhess: read QMDFF reference TURBOMOLE hessian
!
!     part of QMDFF
!
subroutine rdhess(nat3,h,fname)
use qmdff
implicit none
integer::nat3
real(kind=8)::h(nat3,nat3)
character(len=*)::fname
integer::iunit,i,j,mincol,maxcol,idum,jdum
character(len=5)::adum
character(len=80)::a80
logical::driver

driver=.false.

write(10,*)
write(10,*) 'reading <',trim(fname),'>'
if(index(fname,'driv').ne.0) driver=.true.

iunit=11
open(unit=iunit,file=fname)
!
!     read in the hessian from file
!
do
   read(iunit,'(a)')a80
   if(index(a80,'hessian').ne.0)then
      do  i=1,nat3
         maxcol = 0
         200 continue
         mincol = maxcol + 1
         maxcol = min(maxcol+5,nat3)
         if (driver) then
            read(iunit,*)(h(i,j),j=mincol,maxcol)
         else
            read(iunit,'(a5,5f15.10)')adum,(h(i,j),j=mincol,maxcol)
         end if
         if (maxcol.lt.nat3) goto 200
      end do
      close(iunit,status='keep')
      exit
   end if
end do
!
!     scale it by global defined value
!
300 continue  
h = h * scalh

return
end subroutine rdhess
