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
!     subroutine splitmol: Splitting of molecule parts in QMDFF reference
!
!     part of QMDFF
!

subroutine splitmol(n,ic,xyz,wbo,mlist,mol1,iat1,iat2)
use qmdff
implicit none
real(kind=8)::xyz(3,n)
real(kind=8)::wbo(n,n)
integer::n,ic(n)
integer::nn(21,10000)
integer::mlist(*)
integer::mol1,iat1,iat2
integer::i,j,k,m,ibond,imin,jmin
real(kind=8)::r,rcova,f,wbocut

f=1.2
wbocut=0.5

nn=0
!
!     determine covalent neighbours
!
do i=1,n
   k=0
   do j=1,n
      if (i.eq.j) cycle
!
! excluded bond between iat1-iat2
!
      if (i.eq.iat1.and.j.eq.iat2) cycle
      if (i.eq.iat2.and.j.eq.iat1) cycle
      r=sqrt((xyz(1,i)-xyz(1,j))**2+ &
         &   (xyz(2,i)-xyz(2,j))**2+ &
         &   (xyz(3,i)-xyz(3,j))**2)*0.52917726
      rcova=rad(ic(i))+rad(ic(j))
      if (r.lt.f*rcova.and.wbo(j,i).gt.wbocut.and.k.lt.20) then
         k=k+1
         nn(k,i)=j
         nn(21,i)=k
      end if
   end do
end do
!
!     do it
!
mlist(1:n) = 0
mlist(1)=1
mol1=1
!
!     repeat the check at least n times
!
do m=1,n
   do 100 i=1,n
      if(mlist(i).eq.1)cycle
!
!     atom i is not in mol1 and has some partners
!
      do j=1,nn(21,i)
         ibond=nn(j,i)
         do k=1,n
!
!     one partner is in mol1, i is therfore also in mol1
!
            if (mlist(k).eq.1.and.ibond.eq.k) then
               mol1=mol1+1
               mlist(i)=1
               goto 100
            end if
         end do
      end do
   100 continue
end do

return
end subroutine splitmol
