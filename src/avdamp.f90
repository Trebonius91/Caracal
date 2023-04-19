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
!    subroutine avdamp: average damping of QMDFF torsion!
!      depending on atomic distances
!
!    part of QMDFF
!
subroutine avdamp(n,at,nb,xyz,ii,jj,damp)
implicit none
integer n,at(n),ii,jj,nb(20,n),m,kk,ll,ni,nj
real*8 dampki,dampij,dampjl,dum,rij,rik,rjl
real*8 xyz(3,n),damp

rij=(xyz(1,ii)-xyz(1,jj))**2 &
 & +(xyz(2,ii)-xyz(2,jj))**2 &
 & +(xyz(3,ii)-xyz(3,jj))**2
call abdamp(at(ii),at(jj),rij,dampij,dum)

damp=0
m   =0
do ni=1,nb(20,ii)
   kk=nb(ni,ii)
   rik=(xyz(1,kk)-xyz(1,ii))**2 &
    & +(xyz(2,kk)-xyz(2,ii))**2 &
    & +(xyz(3,kk)-xyz(3,ii))**2
   call abdamp(at(kk),at(ii),rik,dampki,dum)
   do nj=1,nb(20,jj)
      ll=nb(nj,jj)
      rjl=(xyz(1,ll)-xyz(1,jj))**2 &
       & +(xyz(2,ll)-xyz(2,jj))**2 &
       & +(xyz(3,ll)-xyz(3,jj))**2
      call abdamp(at(jj),at(ll),rjl,dampjl,dum)
      damp = damp + dampki*dampij*dampjl
      m=m+1
   enddo
enddo
damp=damp/m

return
end subroutine avdamp
