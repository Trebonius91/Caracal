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
!     subroutine rdfix: add user defined restrains for FF calculation
!
!     part of QMDFF
!

subroutine rdfix(nat)
use qmdff
implicit none
real(kind=8)::xyz(3,nat),xx(10),fc,rij(3)
integer::n,nn,i,j,nat,k,fixed(nat)
character(len=128)::line
character(len=40)::fname
character(len=2)::asym
logical::f,ex

fc = 0.05
fixed = 0
nrs = 0

fname='coord.fix'
inquire(file=fname,exist=ex)

if(.not.ex)return

write(*,'(''Reading '',A,'' ...'')') trim(fname)
open(unit=3,file=fname)
n=0
k=0
do
   read(3,'(a)',end=200)line
   if(index(line,'$user').ne.0)goto 300
   if(index(line,'$red' ).ne.0)goto 300
   if(index(line,'$end' ).ne.0)goto 300
   if(index(line,'$fix').ne.0) then
      call readl(line,xx,nn)
      write(*,*) 'force constant set to',fc
      fc=xx(1)
   end if
   call readl(line,xx,nn)
   if (nn.ne.3) cycle ! go back to begin of loop
   n=n+1
   xyz(1,n)=xx(1)
   xyz(2,n)=xx(2)
   xyz(3,n)=xx(3)
   call getf(line,f)
   if (f) then
      fixed(n)=1
      k=k+1
   end if
end do
200  continue

close(3)
stop 'Error in rd!'
300  close(3)

fc=fc/(k-1)

k=0
do i=1,n
   do j=1,i-1
      if (fixed(i).eq.1.and.fixed(j).eq.1) then
         k=k+1
         if (k.gt.10000) stop 'too many restrains'
         rij=xyz(:,j)-xyz(:,i)
         vrest(1,k)=sqrt(sum(rij*rij))
         vrest(2,k)=fc
         rest(1,k)=i
         rest(2,k)=j
      end if
   end do
end do
nrs=k
write(*,*) '# of harmonic restrains ',nrs
write(*,*) 'force constant used     ',fc

return

end subroutine rdfix
