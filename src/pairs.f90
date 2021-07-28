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
!     subroutine pairs: identify a pair of bonded atoms in QMDFF ref. structure 
!     loop over all atom combinations and proove it
!
!     part of QMDFF
!
subroutine pairs(n,nn,list,pair,tag)
implicit none
integer::n,nn(n),list(1000,n),tag
integer::i,j,k,ni,nj,ii,jj
integer::pair(n*(n+1)/2)
logical::dai,daj

k=0
do i=1,n
   ni=nn(i)
   do j=1,i
      k=k+1
      if (i.eq.j) cycle
      nj=nn(j)
      dai=.false.
      daj=.false.
      do ii=1,ni
         if(list(ii,i).eq.j)daj=.true.
      end do
      do jj=1,nj
         if(list(jj,j).eq.i)dai=.true.
      end do
      if (dai.and.daj.and.pair(k).eq.0) then
         pair(k)=tag
      end if
   end do
end do

return
end subroutine pairs

