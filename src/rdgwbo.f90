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
!     subroutine rdgwbo: Open Gaussian output with Wiberg-Mayer bond orders 
!        for QMDFF generation and read them in 
!
!     part of QMDFF
!
subroutine rdgwbo(n,wbo,ok,fn)
implicit none
integer::n,m
logical::ok,ex
real(kind=8)::wbo(n,n)
character(len=128)::atmp
character(len=80)::fn,a80 
real(kind=8)::xx(20)
integer::nn,k,i1,i2,kk,i,j,ll
integer::mincol,maxcol
character(len=40)::adum
ok=.false.
inquire(file=fn,exist=ex)
if(.not.ex) then
   write(10,*) "ERROR: No Gaussian.log-file found! It must have the", &
     &   "same prefix like the Gaussian.chk file!"  
   return
end if
write(10,*) '========================================'
write(10,*) 'reading Wiberg/Mayer BO from ',trim(fn) 
write(10,*) '========================================'     
wbo = 0
open(unit=11,file=fn)
!
!     define number of 9 column lines
!
!
m=(3*n)/9
if (mod(3*n,9).gt.0) m=m+1

100 continue
read(11,'(a)',end=300)a80
if (index(a80,'Wiberg bond index matrix in the NAO basis:').ne.0) then
   maxcol = 0
   do k=1,m
      read(11,'(a)',end=300) a80
      read(11,'(a)')a80
      read(11,'(a)')a80
      mincol = maxcol + 1
      maxcol = min(maxcol+9,n)
      do i=1,n
         read(11,*) adum,adum,(wbo(i,j),j=mincol,maxcol)
      end do
!
!     Fix for Gaussian bug (too small bond orders): Read in only the first block
!
      goto 300
   end do
end if
goto 100
!
!     scale it by global defined value
!

300 continue
close(11)

if (n.eq.2.and.wbo(1,2).lt.0.5)then
   write(10,*) 'warning: setting small diatomic WBO to unity'
   wbo(1,2)=1
   wbo(2,1)=1
endif

end subroutine rdgwbo
