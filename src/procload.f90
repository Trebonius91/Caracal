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
!     subroutine procload: determine number of OMP processes
!      currently not needed, because no OMP is used
!
!     part of QMDFF
!

subroutine procload(nproc,ntask,ioff,bal,ise)
implicit none
integer nproc,ise(2,*),ntask,ioff
integer :: ni(nproc)
integer i,ntot,nmiss,istart
real*8 wi,ws,bal
if (nproc.gt.ntask)then
   ise(1:2,1:nproc)=0
   do i=1,ntask
      ise(1,i)=i
      ise(2,i)=i
   enddo
   return
end if
ws=0
do i=1,nproc
   wi=1.+bal*float(i-1)**1.
   ws=ws+wi
end do
ntot=0
do i=1,nproc
   wi=1.+bal*float(i-1)**1.
   wi=wi/ws
   ni(i)=idint(wi*ntask)
   ntot=ntot+ni(i)
end do
100 continue 
nmiss=ntask-ntot
if(nmiss.lt.1) goto 200
do i=1,nmiss
   if(nmiss.le.nproc)ni(i)=ni(i)+1
end do
ntot=sum(ni(1:nproc))
goto 100
200  continue
istart=1
do i=1,nproc
   ise(1,i)=istart+ioff
   ise(2,i) =istart+ni(i)-1+ioff
   istart=istart+ni(i)
end do

return
end subroutine procload

