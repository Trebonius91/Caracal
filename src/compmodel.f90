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
!     subroutine compmodel: calculate a torsional profile with EHT 
!         used in FF setup: all atoms outside the 4 torsion atoms are 
!         replaced with H atoms! (directly bonded)
!
!     part of QMDFF
!

subroutine compmodel(echo,plot,ii,jj,n,xyz,at,q,nb, &
              & r0ab,zab,r0094,sr42,c66ab,wbo, &
              & nfit,ok,barr,mad,s,o,p,iwt,wt, &
              & cring,ringsize)
use qmdff
implicit none
integer::n,at(n),nb(20,n),ii,jj,nfit,iwt(2,10)
integer::cring(8,n),ringsize(n)
real(kind=8)::xyz(3,n),q(n),barr,s(*),o(*),p(*),mad,wt(10),wbo(n,n)
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),c66ab(n,n)
logical::ok,echo,plot

integer::nbf,nel,it,mm,m,ni,nj,mmm,kk,ll
integer::iii,jjj,kkk,lll,i,ndata
real(kind=8)::e,e0,e1,tangle,angle,c0,minbarr,r,damp,gscal
real(kind=8)::phi,phi0,c1,c2,c3,et,maxbarr,dx(3),valijkl

parameter (ndata=20)
real(kind=8)::fft1(ndata),fitted(ndata)

real(kind=8)::coord(3,100),qm(100)
real(kind=8)::dda(3),ddb(3),ddc(3),ddd(3)
integer::atm(100),map(100),imap(n),list(100),pair(100*101/2)
logical::okbas,pr,err

!   n1 --> nrot1,  n2--> nrot2

barr=0
mad=1
wt=0
iwt=0
gscal=0

tangle=(360.0/dble(ndata-1))*pi/180.
!
!     setup of model system 
!
call modelgeo(n,at,nb,xyz,q,wbo,coord,qm,atm,ii,jj,map,mm,mmm, &
        & cring,ringsize,list,err)
!
!     the model structure is strange: go back to main
!
if (err) then
   ok=.false.
   return
end if

imap=0
do i=1,mm
   imap(map(i))=i
end do
!
!     covalent topology for model system         
!     
call covbond(mm,atm,coord,pair)
!     call prmati(6,pair,mm,0,'pairs')
!
!     average damping of torsion    
!  
call avdamp(n,at,nb,xyz,ii,jj,damp)

!
!     rotate around ii-jj for 360 deg. , initialize basis   
!            
call basis0(mm,atm,nel,nbf)
!
!     if the cut-out is a radical, make cation out of it   
!   
if (mod(nel,2).ne.0) then
   nel=nel-1
   write(*,*) 'warning: open-shell model structure in EHT torsion'
   write(*,*) ' calculation. replaced by cation for that purpose'
end if
!
!     go into EHT calculation
!
call basis (mm,atm,nbf,okbas)           
if (.not.okbas) return
pr=.false.
if (ii.eq.nrot1.and.jj.eq.nrot2) then
   call outgeo(mm,atm,coord,0.0d0)
endif
call eht(mm,atm,coord,qm,nel,nbf,e0,gscal,list,pair)
!
!     determine weights
!
m=0
do ni=1,nb(20,ii)
   kk=nb(ni,ii)
   if (kk.eq.ii) cycle
   if (kk.eq.jj) cycle
   do nj=1,nb(20,jj)
      ll=nb(nj,jj)
      if (ll.eq.jj) cycle
      if (ll.eq.ii) cycle
      m=m+1
      iwt(1,m)=kk
      iwt(2,m)=ll
!           wt(m)=sqrt(dble(z(at(kk))*z(at(ll))))
!     slightly better without weigthing
!
      wt(m)=1
   end do
end do
iwt(1,10)=m
e=sum(wt(1:m))
wt(1:m)=dble(m)*wt(1:m)/e

fft1(1)=0
!
!     define highest rotation barrier possible
!
maxbarr=-1000.
minbarr= 1000.
do m=2,ndata     
   call rotfrag(mm,1,2,coord,tangle,list)
   pr=.false.
   call eht(mm,atm,coord,qm,nel,nbf,e,gscal,list,pair)
   if (ii.eq.nrot1.and.jj.eq.nrot2) then
      call outgeo(mm,atm,coord,(e-e0)*627.51)
   end if
   fft1(m)=(e-e0)*627.51
   if (fft1(m).gt.maxbarr) maxbarr=fft1(m)
   if (fft1(m).lt.minbarr) minbarr=fft1(m)
end do
barr=maxbarr-minbarr
!
!     fit the QMDFF-potential to EHT reference
!
call pfit(ndata,fft1,tangle,damp,s,o,p,mad,fitted,nfit)

p(1:nfit)=p(1:nfit)/627.51
!
!     reasonable barrier and fit ok?  
!    
if (barr.gt.0.and.mad.lt.1.0) ok=.true.
!
! second check on comp/fitted barrier     
! 
maxbarr=-1000.
minbarr= 1000.
do m=1,ndata     
   if (fitted(m).gt.maxbarr) maxbarr=fitted(m)
   if (fitted(m).lt.minbarr) minbarr=fitted(m)
end do
c0=maxbarr-minbarr
!
!     if unusual behavior occured: print infos and stop program
!
if(barr.gt.1.and.abs(c0-barr)/barr.gt.1) ok=.false.

if (ii.eq.nrot1.and.jj.eq.nrot2) then
   do m=1,ndata
      write(101,*) (m-1)*tangle*180/pi,fft1(m)
   end do
   write(101,*)
   do m=1,ndata
      write(101,*) (m-1)*tangle*180/pi,fitted(m)
   end do  
   write(*,*) "ERROR in rotation barrier calculation! stop program---"
   stop
end if

return
end subroutine compmodel
