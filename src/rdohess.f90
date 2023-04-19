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
!     subroutine rdohess: read orca QMDFF reference hessian
!
!     part of QMDFF
!
subroutine rdohess(nat3,h,fname)
use qmdff
implicit none
integer::nat3
real(kind=8)::h(nat3,nat3)
character(len=*)::fname
integer::iunit,i,j,mincol,maxcol,idum,m,k
character(len=5)::adum
character(len=80)::a80
logical::snforca

snforca=.true.
!write(10,*)
!write(10,*) 'reading <',trim(fname),'>'
iunit=12
open(unit=iunit,file=fname)
!
!     define number of 6 column lines
!
m=nat3/5
if (mod(nat3,5).gt.0) m=m+1

!
!     read in the hessian from file
!

100 continue  
read(iunit,'(a)',end=300)a80
if (index(a80,'$act_energy').ne.0) snforca=.false.
if (index(a80,'$hessian').ne.0) then
   read(iunit,'(a)',end=300) a80
   maxcol = 0
   do k=1,m
      read(iunit,'(a)')a80
      mincol = maxcol + 1
      maxcol = min(maxcol+5,nat3)
      if (snforca)then
         do i=1,nat3
            read(iunit,*) adum,(h(i,j),j=mincol,maxcol)
         end do
      else
         do i=1,nat3
            read(iunit,*) idum,(h(i,j),j=mincol,maxcol)
         end do
      end if
   end do
end if
goto 100

!
!     scale it by global defined value
!

300 continue  
close(iunit,status='keep')

h=h*scalh
return
end subroutine rdohess
