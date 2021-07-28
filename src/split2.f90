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
!     subroutine split2: find breakpoints between fragments
!      in QMDFF reference structure 
!
!     part of QMDFF
!
subroutine split2(n,at,nn,nlist,mlist,nm)
use qmdff
implicit none
integer n,nlist(*),mlist(*),nm,nn(20,*),at(*)

integer i,j,k,m
integer mol1,ibond
!
!     do it: define a first molecule 

mol1=1
!
!     repeat the check at least n times: are all 
!     other atoms also part of that molecule?
do m=1,n
   do 100 i=1,n
      if(mlist(i).eq.1) cycle
      if(mlist(i).eq.-1)cycle
      if(metal(at(i)).eq.1)cycle
!
!     atom i is not in mol1 and has some partners
!
      do j=1,nn(20,i)
         ibond=nn(j,i)
         if(metal(at(ibond)).eq.1)cycle
         do k=1,n
!
!     one partner is in mol1, i is therfore also in mol1
!
            if (mlist(k).eq.1.and.ibond.eq.k)then
               mol1=mol1+1
               mlist(i)=1
               goto 100
            end if
         end do
      end do
   100 continue
end do
!
!     List the atoms of the current molecule
!
mol1=0
do i=1,n
   if(mlist(i).eq.1)then
      mol1=mol1+1
      nlist(mol1)=i
      mlist(i)=-1
   endif
enddo
nm=mol1

return
end subroutine split2
