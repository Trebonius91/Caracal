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
!     subroutine fragmentation: define molecular fragments: number 
!        and indices of atoms
!
!     part of QMDFF
!
subroutine fragmentation(n,at,nb,molist)
use qmdff
implicit none
integer::n,nb(20,n),nbmod(20,n),at(n),molist(n)
integer::mlist(n),nlist(n),nmol,i,j,k

nmol=0
mlist=0
mlist(1)=1

nbmod=nb
do i=1,n
   if (metal(at(i)).eq.1) nbmod(20,i)=0
enddo

do
   call split2(n,at,nbmod,nlist,mlist,j)
   if (j.lt.1) exit
   do i=1,n
      if (mlist(i).eq.0) then
          mlist(i)=1
          exit
      end if
   end do
   nmol=nmol+1
   do k=1,j
      molist(nlist(k))=nmol
   end do
end do

write(10,'(I3,A)') nmol,' molecular fragments found'
write(10,*)
write(10,*) '------------------------------------------'
write(10,*) '  #        belongs to fragment'
write(10,*) '------------------------------------------'
do i=1,n
   write(10,'(I3,A,I2)') i,"               ",molist(i)
end do
write(10,*) '------------------------------------------'
write(10,*) 
return
end subroutine fragmentation
