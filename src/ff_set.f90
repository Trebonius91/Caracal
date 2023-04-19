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
!     subroutine ff_set: define a QMDFF amd initialize the forceconstants
!
subroutine ff_set(n,xyz,xyzvdw,at,nb,vdw,q,chb,cxb,pair, &
   &           r0ab,zab,r0094,sr42,cc6ab,echo,wbo,hyb, &
   &           cring,ringsize,mlist,mtors)
use qmdff
implicit none  
integer::at(n),n,nb(20,n),pair(n*(n+1)/2),hyb(n)
integer::cring(8,n),ringsize(n),mlist(n),mtors
real(kind=8)::xyz(3,n),q(n),chb(94),cxb(94),wbo(n,n)
real(kind=8)::xyzvdw(3,n)
logical::vdw,ok,echo,pr,ring
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),cc6ab(n,n)

integer::iat,i,j,k,l,neig,ii,jj,kk,ll,m,iperm,iipp(3,3),ni,nj
integer::mol1,mm,iz1,iz2,nh,oo,nn,tag2(ndim),jk,jl,nxb
integer::lina,ij,tag(n*(n+1)/2),ii1,ii2,i1,i2,i3,pqn,irow
integer::iwt(2,10,n*(n+1)/2),iiwt(2,10),mem1(n),mem2(n),jrow
real(kind=8)::dum,dum1,dum2,dum3,rco,ri,rj,barr,xn,tpot,error,bav,aaa
real(kind=8)::s(ntterm),p(ntterm),o(ntterm),pp(ntterm,n*(n+1)/2),rav,bbb
real(kind=8)::angle,aa1,aa2,phi,rmax1,rmax2,rmin2,c,bb,r,fitmad,valijk,bthr
real(kind=8)::athr2,sthr1,sthr2,valijkl,wt(10,n*(n+1)/2),wwt(10),r2,omega
logical::checkpair,checktrip,checkfour,l1,l2,l3,cdbi,linear
logical::checkinvstr,doeht,samering,metaltorsion,cdbj,skiplin13

!
!     symmetrize to linear. if angle is within this value linear                       
!     this was 1.0 deg in the published version but has to be
!     increased for polyethiine units, SG Wed Nov 26 13:49:50 CET 2014
!
bthr=5.0d0
!
!     inlcude linear 1,3-stretch terms (if .F.)
!
skiplin13=.false.

write(10,*) 
write(10,*)
write(10,*)'===================='
write(10,*)'FF setup: define all'
write(10,*)'QMDFF energy terms'
write(10,*)'===================='
write(10,*)
write(10,*) 'General parameters:'
write(10,*)
write(10,*)' parameter bthr: below which degree an angle is taken as liner?'
write(10,*) bthr 
write(10,*)
write(10,*)' parameter skiplin13: if set F, all linear 1,3-strech terms are included!' 
write(10,*) skiplin13
write(10,*)
iipp(1,1)=1
iipp(1,2)=2
iipp(1,3)=3
iipp(2,1)=2
iipp(2,2)=1
iipp(2,3)=3
iipp(3,1)=3
iipp(3,2)=2
iipp(3,3)=1
vdw=.false.
!
!      additional torsions around metal centers?      
!
metaltorsion=.false.
if (mtors.ne.0) metaltorsion=.true.
 
write(10,*)' Are additional torsions around metal centers included? (True/False:)'
write(10,*) metaltorsion 
write(10,*)
nangl=0
ntors=0
nhb  =0
!
!     do we have a vdW complex?     
! 
call splitmol(n,at,xyz,wbo,mlist,mol1,0,0)
if(mol1.ne.n) then
   write(10,*) ' found vdW complex: more than one molecule!'    
   write(10,*)  
!          
!     make a geometry with separated monomers!     
!     shift all coordinates around 1000!! 
!   
   vdw=.true.
   do i=1,n
      if(mlist(i).eq.1)then
         xyzvdw(1:3,i)=xyz(1:3,i)
      else
         xyzvdw(1:3,i)=xyz(1:3,i)+10000.
      endif
   enddo
else
   mlist=0
endif
!      
!     stretch. set in ff_bond.f!!!!!
!
do k=1,nbond
   i=bond(1,k)
   j=bond(2,k)
   r=sqrt((xyz(1,i)-xyz(1,j))**2 &
     &   +(xyz(2,i)-xyz(2,j))**2 &
     &   +(xyz(3,i)-xyz(3,j))**2)
   vbond(1,k)=r           
   vbond(2,k)=1.0/r/vbond(3,k)**2
end do
nbond12=nbond
!
!     add Hbonds
!
nhb=0   
k  =0
do i=1,n
   do j=1,i
      k=k+1
      if(i.eq.j)cycle
      if(pair(k).lt.3.and.pair(k).ne.0) cycle
      iz1=at(i)
      iz2=at(j)
      c=chb(iz1)*chb(iz2)
      if (c.gt.1.d-6) then
         do nh=1,n
            if(at(nh).ne.1) cycle
            ri=sqrt((xyz(1,i)-xyz(1,nh))**2 &
              &    +(xyz(2,i)-xyz(2,nh))**2 &
              &    +(xyz(3,i)-xyz(3,nh))**2)
            rj=sqrt((xyz(1,j)-xyz(1,nh))**2 &
              &    +(xyz(2,j)-xyz(2,nh))**2 &
              &    +(xyz(3,j)-xyz(3,nh))**2)
            dum1=1.3*(rad(at(i))+rad(1))/0.52917726
            dum2=1.3*(rad(at(j))+rad(1))/0.52917726
            if (ri.lt.dum1.or.rj.lt.dum2)then
               nhb=nhb+1
               hb(1,nhb)=i
               hb(2,nhb)=j
               hb(3,nhb)=nh
               call hbpara(10.0d0,5.0d0,q(i),c)
               vhb(1,nhb)=c*chb(iz1)
               call hbpara(10.0d0,5.0d0,q(j),c)
               vhb(2,nhb)=c*chb(iz2)
            end if
         end do
         100 continue        
      end if
   end do
end do
write(10,'(I3,A)') nhb ,' HB terms (hydrogen bonds)'
!
!     add halogen bonds
!
nxb=0   
k  =0
do i=1,n
   do j=1,i
      k=k+1
      if (i.eq.j)cycle
      if (pair(k).lt.3.and.pair(k).ne.0) cycle
      iz1=at(i)
      iz2=at(j)
      cdbi=iz1.eq.7.or.iz1.eq.8
      cdbj=iz2.eq.7.or.iz2.eq.8
      if (iz1.eq.1.or.iz2.eq.1) cycle
      do nh=1,n
         if (nh.eq.i) cycle
         if (nh.eq.j) cycle
         if (at(nh).eq.17.or.at(nh).eq.35.or. &
           &  at(nh).eq.53.or.at(nh).eq.85) then
            ri=sqrt((xyz(1,i)-xyz(1,nh))**2 &
                &  +(xyz(2,i)-xyz(2,nh))**2 &
                &  +(xyz(3,i)-xyz(3,nh))**2)
            rj=sqrt((xyz(1,j)-xyz(1,nh))**2 &
                &  +(xyz(2,j)-xyz(2,nh))**2 &
                &  +(xyz(3,j)-xyz(3,nh))**2)
            dum1=1.2*(rad(at(i))+rad(at(nh)))/0.52917726
            dum2=1.2*(rad(at(j))+rad(at(nh)))/0.52917726
!
!     Y(I=1)-X..D(J=2)           
!     the charges on X and Y determine the strength
! 
            if (ri.lt.dum1.and.cdbj) then
               nxb=nxb+1
               hb(1,nxb+nhb)=i
               hb(2,nxb+nhb)=j
               hb(3,nxb+nhb)=nh
               call hbpara(-6.5d0,1.0d0,q(nh),c)
               vhb(1,nxb+nhb)=c*cxb(at(nh))
            end if
!
!     Y(j=1)-X..D(I=2)       
!     another configuration
!    
            if (rj.lt.dum2.and.cdbi) then
               nxb=nxb+1
               hb(1,nxb+nhb)=j
               hb(2,nxb+nhb)=i
               hb(3,nxb+nhb)=nh
               call hbpara(-6.5d0,1.0d0,q(nh),c)
               vhb(1,nxb+nhb)=c*cxb(at(nh))
            end if
         end if
      end do
   end do
end do
write(10,'(I3,A)') nxb ,' XB terms (halogen bonds)'
write(10,*)
nhb=nhb+nxb

!
!     bend (bond angle between A-B-C)
!
do i=1,n
   neig=nb(20,i)
   ii=i
!
!     no bending for highly coordinated (metal) atoms 
!     without this skip, cp rings do not rotate
!     highest possible coordination is 6!
!
   if(neig.gt.6) cycle
   do j=1,neig
      do k=1,j-1 
         jj=nb(j,i)
         kk=nb(k,i)
         if (j.eq.k) cycle
         call bangle (XYZ,JJ,II,KK,ANGLE)          
         r=angle*180.0d0/pi
!
!     small angles are not significant and redundant     
!          
         athr2=80.0d0
         l1=metal(at(ii)).eq.1.or. &
            &  metal(at(jj)).eq.1.or. &
            &  metal(at(kk)).eq.1
         if (l1) athr2=60.0d0
         if (r.lt.athr2) cycle
!
!     avoid close to linear situations in highly 
!     coordinated systems
!
         l1= nb(20,jj).gt.3
         l2= nb(20,kk).gt.3
         if (r.gt.160.0d0.and.l1.and.l2) cycle
         do m=1,nangl
            if(checktrip(ii,jj,kk,angl(1,m),angl(2,m),angl(3,m))) &
                 & goto 199
         end do
         nangl=nangl+1
         if (nangl.gt.ndim) stop 'too many angles'
         angl(1,nangl)=ii     
         angl(2,nangl)=jj       
         angl(3,nangl)=kk           
         vangl(1,nangl)=r*pi/180.       
!           
!     symmetrize to linear. without this, close to 
!     linear situations become unstable
!
         if(180.0d0-r.lt.bthr) vangl(1,nangl)=pi  
         vangl(2,nangl)=0.02               
         199  continue
      end do
   end do
end do

!
!     compute for each bond the 
!     TB torsion potential 
!
write(10,*) "Check for rings in reference system..."
do i=1,n
   call getring(n,nb,i,cring(1,i),ringsize(i))
   if (ringsize(i).gt.0) then
      write(10,'('' atom'',i4,'' is in ring '',8i4)') &
         &   i,cring(1:ringsize(i),i)
   end if
end do
write(10,*)
write(10,*) 'computing torsion potentials by Extended Hückel theory'
write(10,*) '   for rotatable bonds ...'
!
!     define the torsion potential.
!
do i=1,ntterm
   o(i)=i
   s(i)=pi
end do

!     else this part is still parity violating, see ffpot_both*
!        do i=1,ntterm/2
!        o(i)=i
!        s(i)=pi
!        enddo
!        do i=ntterm/2+1,ntterm
!        o(i)=i-ntterm/2
!        s(i)=pi/2
!        enddo
!    endif

k=0
do i=1,n
   do j=1,i
      k=k+1
      pp(1:ntterm,k)=0
      ok=.false.
!
!     check cases for which we don't need to compute the potential   
!         
      if (i.eq.j.or.pair(k).ne.1) cycle
      if (nb(20,j).eq.1.or.nb(20,j).gt.4) cycle
      if (nb(20,i).eq.1.or.nb(20,i).gt.4) cycle
      r=sqrt((xyz(1,i)-xyz(1,j))**2 &
         &  +(xyz(2,i)-xyz(2,j))**2 &
         &  +(xyz(3,i)-xyz(3,j))**2)
!
!     check if TB/EHT should be done          
!  
      l1=doeht(n,at,i,j,r,hyb,cring,ringsize,wbo)
      if(l1) then   
         p(1:ntterm)=0.001
!
!     do the TB/EHT pot curve and the fit
!
         call compmodel(echo,.false.,i,j,n,xyz,at,q,nb, &
           &  r0ab,zab,r0094,sr42,c6ab,wbo,ntterm,ok,barr,fitmad, &
           &  s,o,p,iiwt,wwt,cring,ringsize)
!
!     this avoids the TB/EHT potential for conjugated systems like C60     
!        
         if (barr.gt.25) ok=.false.
!
!     this avoids artificial model systems
!
         if (barr.lt.0.05) ok=.false.
         write(10,'('' pair '',2i4,2x &
              &  ''barrier (kcal) :'',f6.2,2x, &
              &  '' %MAD(fit):'',f5.1, &
              &  '' fit pot used?'',l)') &
              &  i,j,barr,fitmad*100.,ok
      end if
      if (ok) then
!
!     normalize by the # of torsions for ij     
!         
         mm=0
         do m=1,ntterm
            pp(m,k)=p(m)/dble((nb(20,i)-1)*(nb(20,j)-1))
         end do
         iwt(1:2,1:10,k)=iiwt(1:2,1:10)
         wt(1:10,k)    =     wwt(1:10)
      end if
   end do
end do
!
!     torsional terms
!
tag=0
do i=1,n
   ii=i
!
!     no torsions if not coordinated or highly coordinated        
!
   if (nb(20,ii).eq.1.or.nb(20,ii).gt.5) cycle
   do j=1,nb(20,ii)
      jj=nb(j,ii)
      ij=lina(jj,ii)
!
!     already done this ij?   
!      
      if (tag(ij).eq.1) cycle
      if (nb(20,jj).eq.1.or.nb(20,jj).gt.5) cycle
      do k=1,nb(20,jj)
         kk=nb(k,jj)    
         if (kk.eq.ii) cycle
         do l=1,nb(20,ii)
            ll=nb(l,ii)
            if (ll.eq.jj.or.ll.eq.kk) cycle

            r=valijkl(n,xyz,ll,ii,jj,kk)*180./pi
!
!     check for linear arrangements   
!             
            call bangle(XYZ,kk,ii,jj,aa1)          
            call bangle(XYZ,ll,ii,jj,aa2)         
            aa1=aa1*180./pi
            aa2=aa2*180./pi
            if (aa1.lt.athr.or.180-aa1.lt.athr) cycle
            if (aa2.lt.athr.or.180-aa2.lt.athr) cycle
!
!     exit if rotational number is 0 (eg 5-membered rings)
!
            if (pp(1,ij).lt.1.d-6) then
               call getrot(n,at,ii,jj,hyb,cring,ringsize,xn)
               if(xn.lt.1.d-6) cycle
            end if

            ntors=ntors+1
            if(ntors.gt.ndim) stop 'too many torsions'
            tors(1,ntors)=ll
            tors(2,ntors)=ii
            tors(3,ntors)=jj
            tors(4,ntors)=kk
            vtors(1,ntors)=r*pi/180.0d0
!
!     fitted pot case (done with EHT above)     
!   
            if (pp(1,ij).ne.0) then
               tors(5,ntors)=ntterm
!
!     no fit
!
               tors(6,ntors)=0    
               vtors(2,ntors)=1
!
!     use computed weights     
!      
               do m=1,iwt(1,10,ij)
                  if (ll.eq.iwt(1,m,ij).and.kk.eq.iwt(2,m,ij)) then
                     vtors(2,ntors)=wt(m,ij)
                  end if
                  if (ll.eq.iwt(2,m,ij).and.kk.eq.iwt(1,m,ij)) then
                     vtors(2,ntors)=wt(m,ij)
                  end if
               end do
               mm=3
               do m=1,ntterm
                  vtors(mm  ,ntors)=o (m)
                  vtors(mm+1,ntors)=s (m)
                  vtors(mm+2,ntors)=pp(m,ij)
                  mm=mm+3
               end do
            else
!
!     fitted to Hessian with one-term (sine: shift=pi) potential  
!               
               call getrot(n,at,ii,jj,hyb,cring,ringsize,xn)
               tors(5,ntors)=1           
               tors(6,ntors)=1     
               vtors(2,ntors)=0.001
               vtors(3,ntors)=xn    
               vtors(4,ntors)=pi 
               vtors(5,ntors)=1.0   
            end if
! 
!     mark torsion as done for this ij pair (avoids doubles)
!
            tag(ij)=1
!
!     next lijk       
!  
         end do
      end do
   end do
end do

!
!     special case of metal atoms   
!    
     
if(metaltorsion)then
   write(10,*)'torsion terms at metal centers with rot=',mtors,&
       &  'switched on'
   do i=1,n
      m=0
      ni=nb(20,i)
      if (ni.gt.4) then 
         do i1=1,ni
            do i2=1,i1-1
               kk=nb(i1,i)
               ll=nb(i2,i)
               if (nb(20,kk).eq.1) cycle
               if (nb(20,ll).eq.1) cycle
               call bangle(xyz,kk,i,ll,angle)          
               r=angle*180/pi
               if (r.gt.100.and.r.lt.160.and.m.lt.6) then
                  jk=nb(1,kk)
                  jl=nb(1,ll)
                  if (.not.checkfour(kk,i,ll,jl)) cycle
                  if (.not.checkfour(ll,i,kk,jk)) cycle
                  m=m+1
                  r=valijkl(n,xyz,kk,i,ll,jl)
                  ntors=ntors+1
                  tors(1,ntors)=kk
                  tors(2,ntors)=i 
                  tors(3,ntors)=ll
                  tors(4,ntors)=jl
                  tors(5,ntors)=1  
                  tors(6,ntors)=1     
                  vtors(1,ntors)=r  
                  vtors(2,ntors)=0.001
!   Grubbs        vtors(3,ntors)=4     
                  vtors(3,ntors)=mtors 
                  vtors(4,ntors)=pi    
                  vtors(5,ntors)=1     
                  r=valijkl(n,xyz,ll,i,kk,jk)
                  ntors=ntors+1
                  tors(1,ntors)=ll
                  tors(2,ntors)=i 
                  tors(3,ntors)=kk
                  tors(4,ntors)=jk
                  tors(5,ntors)=1  
                  tors(6,ntors)=1     
                  vtors(1,ntors)=r  
                  vtors(2,ntors)=0.001
                  vtors(3,ntors)=mtors 
                  vtors(4,ntors)=pi    
                  vtors(5,ntors)=1     
               end if
            end do
         end do
      end if
   end do
end if

!
!     inversion type torsion (central atom has 3 bonds)     
!
do i=1,n
   ii=i
   l1=nb(20,ii).eq.3
   l2=nb(20,ii).eq.4.and.nb(18,ii).eq.1
   l3=nb(20,ii).eq.4.and.metal(at(ii)).eq.1
   if (l1.or.l2.or.l3) then  
      jj=nb(1,ii)
      kk=nb(2,ii)
      ll=nb(3,ii)
      r=omega (n,xyz,ii,jj,kk,ll)
      ntors=ntors+1
      if (ntors.gt.ndim) stop 'too many torsions'
      tors(1,ntors)=ii
      tors(2,ntors)=jj
      tors(3,ntors)=kk
      tors(4,ntors)=ll
      tors(5,ntors)=1  
      tors(6,ntors)=2     
      vtors(1,ntors)=r  
      vtors(2,ntors)=0.1
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
!      the 'if' takes the chiraliy breaking term for planar cases
!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
      if (abs(r*180.0d0/pi).le.athr) then
         vtors(3,ntors)=1   
!  
!     symmetrize it to planarity 
!
         vtors(1,ntors)=0  
      else
!
!     this is the double minimum around 0    
!               
         vtors(3,ntors)=0   
      end if
      vtors(4,ntors)=pi    
      vtors(5,ntors)=1     
   end if
end do

!cccccccccccccccccccccccccccc
!     torsion terms for linear R-A-X-B-R
!     where the torsion is then taken for
!     X-A..B-Y (eg allene)
!cccccccccccccccccccccccccccc
tag2=0
do i=1,nangl
   if (abs(vangl(1,i)*180./pi-180.).le.athr.and.tag2(i).eq.0) then
       ii=angl(1,i)       
       jj=angl(2,i)       
       kk=angl(3,i)       
       do j=1,nb(20,jj)
          if (nb(j,jj).eq.ii) cycle
          do k=1,nb(20,kk)
             if (nb(k,kk).eq.ii) cycle

             r=valijkl(n,XYZ,nb(j,jj),jj,kk,nb(k,kk))   
! 
!     check for linear arrangements             
!   
             call bangle(xyz,nb(j,jj),kk,jj,aa1)          
             call bangle(xyz,nb(k,kk),kk,jj,aa2)         
             aa1=aa1*180./pi
             aa2=aa2*180./pi
             if (aa1.lt.athr.or.180-aa1.lt.athr) cycle
             if (aa2.lt.athr.or.180-aa2.lt.athr) cycle

             ntors=ntors+1
             if (ntors.gt.ndim) stop 'too many torsions'
             tors(1,ntors)=nb(j,jj)
             tors(2,ntors)=jj 
             tors(3,ntors)=kk
             tors(4,ntors)=nb(k,kk)
             tors(5,ntors)=1 
             tors(6,ntors)=1     
             vtors(1,ntors)=r
             vtors(2,ntors)=0.001
             vtors(3,ntors)=2  
             vtors(4,ntors)=pi
             vtors(5,ntors)=1 
          end do
       end do
   end if
end do

!cccccccccccccccccccccccccccc
!     torsion terms for linear A-X-B-R
!     where the torsion is then taken for
!     A...R-Y      
!cccccccccccccccccccccccccccc
do i=1,nangl
   if (abs(vangl(1,i)*180./pi-180.).le.athr.and.tag2(i).eq.0) then
      ii=angl(1,i)       
      jj=angl(2,i)       
      kk=angl(3,i)       
      ll=0
      if (nb(20,jj).eq.1) then
         ll=kk
         mm=jj
      end if
      if (nb(20,kk).eq.1) then
         ll=jj
         mm=kk
      end if
      if(ll.eq.0) goto 200
      do j=1,nb(20,ll)
         nn=nb(j,ll)
         if (nn.eq.ll) cycle
         if (nn.eq.ii) cycle
         do k=1,nb(20,nn)
            oo=nb(k,nn)
            if (oo.eq.ll) cycle
            if (oo.eq.nn) cycle
            r=valijkl(n,xyz,mm,ii,nn,oo)   
! 
!     check for linear arrangements  
!              
            call bangle(xyz,mm,ii,nn,aa1)          
            call bangle(xyz,oo,nn,ii,aa2)         
            aa1=aa1*180./pi
            aa2=aa2*180./pi
            if (aa1.lt.athr.or.180-aa1.lt.athr) cycle
            if (aa2.lt.athr.or.180-aa2.lt.athr) cycle

            tag2(i)=1
            ntors=ntors+1
            if (ntors.gt.ndim) stop 'too many torsions'
            tors(1,ntors)=mm         
            tors(2,ntors)=ii 
            tors(3,ntors)=nn 
            tors(4,ntors)=oo       
            tors(5,ntors)=1 
            tors(6,ntors)=1     
            vtors(1,ntors)=r
            vtors(2,ntors)=0.001
            vtors(3,ntors)=4 
            vtors(4,ntors)=pi
            vtors(5,ntors)=1 
         end do
      end do
   end if
   200  continue
end do

!ccccccccccccccccccccccccccccccccccccc
!     stretch for bend: only define parameters!
!ccccccccccccccccccccccccccccccccccccc
88  continue 
rmax1=-1


! for fit details see ff_bond.f
!     open(unit=66,file='~/.param')
!     read(66,*) aaa,bbb
!     read(66,*) ade13(1),ade13(5:9)
!     close(66)

aaa=2.81
bbb=0.53

l=nbond
do i=1,nangl
   ii=angl(1,i)       
   jj=angl(2,i)       
   kk=angl(3,i)       
   linear=abs(vangl(1,i)*180./pi-180.).le.athr
   if (linear) then
!
!     its not entirely clear if 1,3-str should be kept for
!     linear cases. In the published version they were
!     discarded and if included e.g. for c2h2 they are
!     automatically removed (FC<0). for other cases this is, 
!     however, different
!
      if (skiplin13) cycle
   end if
   l1=metal(at(ii)).eq.1.or. &
     &  metal(at(jj)).eq.1.or. &
     &  metal(at(kk)).eq.1
   r=sqrt((xyz(1,kk)-xyz(1,jj))**2 &
     &  +(xyz(2,kk)-xyz(2,jj))**2 &
     &  +(xyz(3,kk)-xyz(3,jj))**2)*0.52917726
   if (mlist(jj).eq.mlist(ii).and.r.gt.rmax1) rmax1=r         
   rco=rad(at(kk))+rad(at(jj))
!
!     already there       
!  
   do m=1,l
      if(jj.eq.bond(1,m).and.kk.eq.bond(2,m)) goto 299
      if(jj.eq.bond(2,m).and.kk.eq.bond(1,m)) goto 299
   end do
   l=l+1
   bond(1,l)=jj
   bond(2,l)=kk
!
!     both neighbors are in a ring  
!         
   l1=ringsize(kk).gt.0.and.ringsize(jj).gt.0
!
!     De not sensitive to this if be values (ff_bond) are adjusted
!
   if (l1) then
      vbond(3,l)=aaa+bbb*ade13(at(jj))*ade13(at(kk))+0.7 
   else
      vbond(3,l)=aaa+bbb*ade13(at(jj))*ade13(at(kk))
   end if
!
!     adjusted to c4h2  
!      
   if (linear) vbond(3,l)=vbond(3,l)+2.5
   vbond(1,l)=r/0.52917726
   vbond(2,l)=1.0/vbond(1,l)/vbond(3,l)**2
   299  continue
end do
nbond=l

write(10,*)
write(10,'('' largest distance considered as covalent stretch (r(bohr), r(Ang)):'')') 
write(10,'(2F8.2)')  rmax1/0.52917726,rmax1
write(10,*) 

!ccccccccccccccccccccccccccccccccccccc
!     symmetrize (find out which bond is rotable)   
!ccccccccccccccccccccccccccccccccccccc

!
!     larger values (higher chance for symmetrization) is better   
!   
sthr1=1.2
sthr2=0.2
!     first get atoms which become equivalent upon single-bond rotation      
iwt=0
k=0
do i=1,n
   iz1=at(i)
   do j=1,i
      iz2=at(j)
      k=k+1
!
!     exit cases    
!        
      if (i.eq.j.or.pair(k).ne.1)         cycle
      if (wbo(j,i).gt.sthr1     )         cycle
      if (nb(20,j).eq.1.or.nb(20,j).gt.4) cycle
      if (nb(20,i).eq.1.or.nb(20,i).gt.4) cycle
      if (nb(18,i).eq.1.or.nb(18,j).eq.1) cycle
      if (metal(iz1).eq.1.or.metal(iz2).eq.1) cycle
      if (samering(n,i,j,cring,ringsize)) cycle           
!     on i atom
      mem1=0
      do ni=1,nb(20,i)
         ii1=nb(ni,i)
         do nj=1,nb(20,i)
            ii2=nb(nj,i)
            if (ii1.eq.ii2) cycle
            if (at(ii1).ne.at(ii2)) cycle
            if (abs(wbo(i,ii1)-wbo(i,ii2)).gt.sthr2) cycle
            if (ii1.eq.i.or.ii2.eq.i) cycle
            if (ii1.eq.j.or.ii2.eq.j) cycle
            mem1(ii1)=mem1(ii1)+1
            mem1(ii2)=mem1(ii2)+1
         end do
      end do
      mm=0
      do kk=1,n
         if (mem1(kk).gt.0) then
            mm=mm+1
            iwt(1,mm,k)=kk
            iwt(1,10,k)=mm
         end if
      end do
!     on j atom
      mem2=0
      do ni=1,nb(20,j)
         ii1=nb(ni,j)
         do nj=1,nb(20,j)
            ii2=nb(nj,j)
            if (ii1.eq.ii2) cycle
            if (at(ii1).ne.at(ii2)) cycle
            if (abs(wbo(j,ii1)-wbo(j,ii2)).gt.sthr2) cycle
            if (ii1.eq.j.or.ii2.eq.j) cycle
            if (ii1.eq.i.or.ii2.eq.i) cycle
            mem2(ii1)=mem2(ii1)+1
            mem2(ii2)=mem2(ii2)+1
         end do
      end do
      mm=0
      do kk=1,n
         if (mem2(kk).gt.0) then
            mm=mm+1
            iwt(2,mm,k)=kk
            iwt(2,10,k)=mm
         end if
      end do
   end do
end do

write(10,*) 'local symmetry analysis (rotable bonds) ...'
k=0
do i=1,n
   do j=1,i
      k=k+1
      kk=iwt(1,10,k)
      if (kk.gt.0)then
         write(10,'('' rotation on atom'',i4,'' symmetrizes:'',9i4)') &
           &  i,iwt(1,1:kk,k)
!
!     bonds and 1,3 stretches        
!                           
         rav=0
         bav=0
         do mm=1,kk
            jj=iwt(1,mm,k)
            r=sqrt((xyz(1,i)-xyz(1,jj))**2 &
              &   +(xyz(2,i)-xyz(2,jj))**2 &
              &   +(xyz(3,i)-xyz(3,jj))**2)
            rav=rav+r
            r=sqrt((xyz(1,j)-xyz(1,jj))**2 &
              &   +(xyz(2,j)-xyz(2,jj))**2 &
              &   +(xyz(3,j)-xyz(3,jj))**2)
            bav=bav+r
         end do
         rav=rav/(mm-1)
         bav=bav/(mm-1)
         do mm=1,kk
            jj=iwt(1,mm,k)
            do l=1,nbond
               i1=bond(1,l)
               i2=bond(2,l)
               if (i1.eq.i.and.i2.eq.jj) vbond(1,l)=rav
               if (i2.eq.i.and.i1.eq.jj) vbond(1,l)=rav
               if (i1.eq.j.and.i2.eq.jj) vbond(1,l)=bav
               if (i2.eq.j.and.i1.eq.jj) vbond(1,l)=bav
            end do
         end do
!
!      angles  
!             
         rav=0
         do mm=1,kk
            jj=iwt(1,mm,k)
            call bangle(xyz,jj,i,j,angle)          
            rav=rav+angle
         end do
         rav=rav/(mm-1)
         do mm=1,kk
            jj=iwt(1,mm,k)
            do l=1,nangl
               i1=angl(1,l)
               i2=angl(2,l)
               i3=angl(3,l)
               if (i1.eq.i.and.i2.eq.j.and.i3.eq.jj)vangl(1,l)=rav
               if (i1.eq.i.and.i3.eq.j.and.i2.eq.jj)vangl(1,l)=rav
            end do
         end do
      end if
      kk=iwt(2,10,k)
      if (kk.gt.0)then
          write(10,'('' rotation on atom'',i4,'' symmetrizes:'',9i4)') &
            &  j,iwt(2,1:kk,k)
         rav=0
         bav=0
         do mm=1,kk
            jj=iwt(2,mm,k)
            r=sqrt((xyz(1,j)-xyz(1,jj))**2 &
              &   +(xyz(2,j)-xyz(2,jj))**2 &
              &   +(xyz(3,j)-xyz(3,jj))**2)
            rav=rav+r
            r=sqrt((xyz(1,i)-xyz(1,jj))**2 &
              &   +(xyz(2,i)-xyz(2,jj))**2 &
              &   +(xyz(3,i)-xyz(3,jj))**2)
            bav=bav+r
         end do
         rav=rav/(mm-1)
         bav=bav/(mm-1)
         do mm=1,kk
            jj=iwt(2,mm,k)
            do l=1,nbond
               i1=bond(1,l)
               i2=bond(2,l)
               if (i1.eq.j.and.i2.eq.jj) vbond(1,l)=rav
               if (i2.eq.j.and.i1.eq.jj) vbond(1,l)=rav
               if (i1.eq.i.and.i2.eq.jj) vbond(1,l)=bav
               if (i2.eq.i.and.i1.eq.jj) vbond(1,l)=bav
            end do
         end do
!
!      angles  
!             
         rav=0
         do mm=1,kk
            jj=iwt(2,mm,k)
            call bangle(xyz,jj,j,i,angle)          
            rav=rav+angle
         end do
         rav=rav/(mm-1)
         do mm=1,kk
            jj=iwt(2,mm,k)
            do l=1,nangl
               i1=angl(1,l)
               i2=angl(2,l)
               i3=angl(3,l)
               if (i1.eq.j.and.i2.eq.i.and.i3.eq.jj) vangl(1,l)=rav
               if (i1.eq.j.and.i3.eq.i.and.i2.eq.jj) vangl(1,l)=rav
            end do
         end do
      end if
   end do
end do

write(10,*)
write(10,*) 'force field setup finished!.'
write(10,*) 

return
end subroutine ff_set
