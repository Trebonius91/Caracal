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
!     subroutine rdwbo: Open orca output with Wiberg-Mayer bond orders 
!        for QMDFF generation and read them in 
!
!     part of QMDFF
!
subroutine rdwbo(n,wbo,ok,fn)
implicit none
integer::n
logical::ok,ex
real(kind=8)::wbo(n,n)
character(len=128)::atmp
character(len=80)::fn  
real(kind=8)::xx(20)
integer::nn,k,i1,i2,kk,i,ll

ok=.false.
inquire(file=fn,exist=ex)
if(.not.ex) then
   write(10,*) "ERROR: No orca.out-file found! It must have the", &
     &   "same prefix like the orca.hess-file!"  
   return
end if
write(10,*) '========================================'
write(10,*) 'reading Wiberg/Mayer BO from ',trim(fn) 
write(10,*) '========================================'     
wbo = 0
open(unit=11,file=fn)
do
   read(11,'(a)',end=100) atmp     
   if(index(atmp,'Mayer bond orders larger').ne.0) then         
   20 read(11,'(a)',end=100) atmp     
   ok=.true.
   do ll=1,len(atmp)
      if(atmp(ll:ll).eq.',') atmp(ll:ll)=' '
   enddo
   call readl(atmp,xx,nn)
   k=nn/3
   kk=0
   do i=1,k
      kk=kk+1
      i1=idint(xx(kk))+1              
      kk=kk+1
      i2=idint(xx(kk))+1              
      kk=kk+1
      wbo(i1,i2)=xx(kk)
      wbo(i2,i1)=xx(kk)
   enddo
   if (nn.gt.0) goto 20
endif
end do
100 continue 
close(11)


if (n.eq.2.and.wbo(1,2).lt.0.5)then
   write(10,*) 'warning: setting small diatomic WBO to unity'
   wbo(1,2)=1
   wbo(2,1)=1
endif

end subroutine rdwbo
