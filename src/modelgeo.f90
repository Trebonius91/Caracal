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
!     subroutine modelgeo: create model geometry 
!        for EHT torsional scan
!
!     part of QMDFF
!
subroutine modelgeo(n,at,nb,xyz,q,wbo,coord,qm, &
            &    atm,ii,jj,map,mm,mmm, &
            &    cring,ringsize,list,err)
use qmdff
implicit none
integer::n,at(n),atm(*),ii,jj,nb(20,n),mm,mmm,map(*)
integer::cring(8,n),ringsize(n),nbc(50,2),list(*)
integer::ni,nii,nj,njj,mni,mnj,m,nba(50),nbb(50),i,j,iii,l,k
real(kind=8)::xyz(3,n),q(n),wbo(n,n)
real(kind=8)::coord(3,*),qm(*)
real(kind=8)::dx(3),dy(3),r,wbom(100,100)
logical::atomthere,inring,err

err=.false.

map(1)=ii
map(2)=jj
atm(1)=at(ii)
atm(2)=at(jj)
coord(1:3,1)=xyz(1:3,ii)
coord(1:3,2)=xyz(1:3,jj)

nba(50)=nb(20,ii)
nbb(50)=nb(20,jj)
do i=1,nb(20,ii)
   nba(i)=nb(i,ii)
end do
do i=1,nb(20,jj)
   nbb(i)=nb(i,jj)
end do
!
!     for the case that rings are contained in the system
!
if (ringsize(ii).gt.0) then
   do i=1,ringsize(ii)
      iii=cring(i,ii)
      if (atomthere(iii,nba(50),nba)) cycle
      nba(50)=nba(50)+1
      nba(nba(50))=iii
   end do
end if

if (ringsize(jj).gt.0) then
   do i=1,ringsize(jj)
      iii=cring(i,jj)
      if (atomthere(iii,nbb(50),nbb)) cycle
      nbb(50)=nbb(50)+1
      nbb(nbb(50))=iii
   end do
end if

mmm=0
mm=2

do m=1,nba(50)
   ni=nba(m)
   if(atomthere(ni,mm,map))cycle
!
!      this introduces artifacts for ADCBI, does not seem to be necessary
!      if(metal(at(ni)).eq.1)cycle
!
   mm=mm+1
!
!      No larger that 99 atoms!!
!
   if(mm.gt.99) stop 'model system too large'
   mmm=mmm+1
   atm(mm)=at(ni)
   map(mm)=ni
   coord(1:3,mm)=xyz(1:3,ni)
end do

do m=1,nbb(50)
   nj=nbb(m)
   if (atomthere(nj,mm,map)) cycle
!
!      this introduces artifacts for ADCBI, does not seem to be necessary
!        if(metal(at(nj)).eq.1)cycle
!
   mm=mm+1
   if (mm.gt.99) stop 'model system too large'
   map(mm)=nj
   atm(mm)=at(nj)
   coord(1:3,mm)=xyz(1:3,nj)
end do
k=mm

do m=3,k 
   ni=map(m)
!     this removes groups on the rotating atoms
!     instead dangle atoms are replaced by H for simplification..        
!        if(.not.inring(n,ni,ii,jj,cring,ringsize))cycle
   do l=1,nb(20,ni)
      nj=nb(l,ni)
      if (atomthere(nj,mm,map)) cycle
      mm=mm+1
      if (mm.gt.99) stop 'model system too large'
      map(mm)=nj
      atm(mm)=at(nj)
      coord(1:3,mm)=xyz(1:3,nj)
   end do
end do
k=mm

do i=1,k
   qm(i)=q(map(i))
   do j=1,k
      wbom(j,i)=wbo(map(i),map(j))
   end do
end do
!   
!     saturate dangling bonds: determine neighbor atoms   
!      
call nneighbor(mm,coord,atm,wbom,nbc)

do m=1,k
   if(nbc(m,1).eq.0)then
      err=.true.  
      return
   endif
enddo

do m=1,k 
   if(nbc(m,1).ne.nb(20,map(m)).and.atm(m).ne.1)then
      nj=nbc(m,2)
!
!     put the H at a standard distance   
!         
      dx(1:3)=coord(1:3,m)-coord(1:3,nj)
      r=0.5291677*sqrt(sum(dx*dx))/(rad(1)+rad(atm(nj)))
      coord(1:3,m)=coord(1:3,nj)+dx(1:3)/r
      atm  (    m)=1               
      qm   (    m)=0.1             
   end if
end do
!
!     part of left or right side?   
!   
 10   list(1)=1
list(2)=2
do m=3,k
   dx(1:3)=coord(1:3,m)-coord(1:3,1)
   dy(1:3)=coord(1:3,m)-coord(1:3,2)
   if (sum(dx*dx).lt.sum(dy*dy)) then
      list(m)=1
   else
      list(m)=2
   end if
end do

return
end subroutine modelgeo
