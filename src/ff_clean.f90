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
!     subroutine ff_clean: all bonded terms are neglected whose force constants 
!           have too low values
!
!     part of QMDFF
!
subroutine ff_clean(n,xyz,at,refit)
use qmdff
implicit none  
integer::n,at(n)
real(kind=8)::xyz(3,n)
logical::refit
integer::iat,i,j,k,l,m,new,ii(6,ndim),minus,plus
real(kind=8),allocatable :: vv(:,:)
real(kind=8)::thr

allocate(vv(2+3*ntterm,ndim))
refit=.false.
!
!     under this value the force constant an the term will be erased
!
thr=1.d-4  
!
!     the bond terms
!
new=0
write(*,*)
do m=1,nbond
   if (vbond(2,m).gt.thr) then
      new=new+1
      vv(1,new)=vbond(1,m)
      vv(2,new)=vbond(2,m)
      vv(3,new)=vbond(3,m)
      ii(1,new)=bond (1,m)
      ii(2,new)=bond (2,m)
   else
      if (abs(vbond(2,m)).gt.0.0) &
       &  write(10,*) 'discarded bond ',bond(1,m),bond(2,m),vbond(2,m)
   end if
end do
write(10,*) 'bonds  old,new',nbond,new
nbond=new
vbond(1:3,1:new)=vv(1:3,1:new)
bond(1:2,1:new)=ii(1:2,1:new)
!
!     the 1-2-3 angle terms (bend angle)
!
new=0
do m=1,nangl
   if (vangl(2,m).gt.thr) then
      new=new+1
      vv(1,new)=vangl(1,m)
      vv(2,new)=vangl(2,m)
      ii(1,new)=angl (1,m)
      ii(2,new)=angl (2,m)
      ii(3,new)=angl (3,m)
   else
      if (abs(vangl(2,m)).gt.0.1) &
      &  write(10,*)'discarded angle ',angl(1:3,m),vangl(1:2,m)  
   end if
end do
write(10,*) 'angles old,new',nangl,new
nangl=new
vangl(1:2,1:new)=vv(1:2,1:new)
angl(1:3,1:new)=ii(1:3,1:new)

!
!     the torsional angle terms
!     here only positive force constants above the throughput 
!     are considered
!
new=0
minus=0
plus =0
do m=1,ntors
!     take only positive FCs            
   if (vtors(2,m).gt.1.and.tors(4,m).gt.0) then
!      write(*,'(''warning: large FC for torsion'',6i4,f10.2)') &
!            &   tors(1:6,m),vtors(2,m)
      plus=plus+1
   end if
   if (vtors(2,m).gt.thr.and.vtors(2,m).lt.99) then
      new=new+1
      do j=1,3*ntterm+2
         vv(j,new)=vtors(j,m)
      end do
      do j=1,6
         ii(j,new)=tors (j,m)
      end do
   end if
   if (vtors(2,m).lt.0) then
      minus=minus+1
!      write(*,'(''warning: FC<0 for torsion'',6i4)')tors(1:6,m)
   end if
end do
write(10,*) 'tors   old,new',ntors,new
write(10,*) '# tors FCs < 0',minus    
write(10,*) '# tors FCs > 1',plus     
ntors=new
vtors(1:2+3*ntterm,1:new)=vv(1:2+3*ntterm,1:new)
tors(1:6,1:new)=ii(1:6,1:new)

deallocate(vv)

return
end subroutine ff_clean
