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
!     This subroutine ff_anal calculates again the QMDFF energy and distinguishes
!     between the different groups of energy terms
!     It can calculate the parts for single atoms!
!
!     part of QMDFF
!

subroutine ff_anal(n,at,xyz,xyz2,q,r0ab,zab,r0094,sr42,c66ab)
use qmdff
implicit none  
integer::n,at(n)
real(kind=8)::xyz(3,n),xyz2(3,n),q(n)
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),c66ab(n,n)
integer::i1,i2,iz1,iz2,k,lina,kk,i,j,m,ic,nt,it,l,it2,mm
real(kind=8)::r2,r,r4,r6,r06,R0,t6,t8,c6t6,c6t8,t27,drij,c1,c2,aai,aai2
real(kind=8)::dx,dy,dz,c6,x,alpha,oner,aa,damp,damp2,e0,de,rij,kij,et
real(kind=8)::c0,kijk,va(3),vb(3),vc(3),dt,phi,phi0,rp,rkl,rjk,rjl
real(kind=8)::ea,cosa,rab2,rcb2,vlen,theta,vdc(3),vcb(3),vab(3),vp(3)
real(kind=8)::dampkl,damp2kl,damp2jl,dampjk,damp2jk,dampij,damp2ij,dampjl
real(kind=8)::valijkl,esum,rn,dphi,phipi,ef,x1cos,x2cos
real(kind=8)::dphi1,dphi2,omega,e1,e2,e3,e4
real(kind=8)::epair (n*(n+1)/2,4),epart(ndim,3)
real(kind=8)::epair2(n*(n+1)/2,4),epart2(ndim,3)
real(kind=8)::ptors (ndim),pbend (ndim)
real(kind=8)::ptors2(ndim),pbend2(ndim)
real(kind=8)::eatom(n),eatom2(n),third,quart,half

third=0.33333333333333333333333333333333d0
quart=0.25000000000000000000000000000000d0
half=0.50000000000000000000000000000000d0

write(*,*)
write(*,*) '=============================='   
write(*,*) 'analyzing energy contributions' 
write(*,*) 'from QMDFF energy terms'            
write(*,*) '=============================='  

epart=0
epart2=0
epair=0
epair2=0
ptors=0
pbend=0
ptors2=0
pbend2=0
eatom=0
eatom2=0
esum=0
!
!     non covalent interactions
!

do k=1,nnci               
   i1=nci(1,k)
   i2=nci(2,k)
   kk=lina(i1,i2)
   dx=xyz(1,i1)-xyz(1,i2)          
   dy=xyz(2,i1)-xyz(2,i2)          
   dz=xyz(3,i1)-xyz(3,i2)          
   r2=dx*dx+dy*dy+dz*dz
   r =sqrt(r2)
   oner =1.0d0/r
   iz1=at(i1) 
   iz2=at(i2) 
   R0=r0094(iz1,iz2)        
   c6=c66ab(i2,i1)
   r4=r2*r2
   r6=r4*r2
   r06=R0**6
   t6=r6+r06
   t8=r6*r2+r06*R0*R0
   c6t6=c6/t6
   c6t8=c6/t8
   t27=sr42(iz1,iz2)*c6t8
   e0=-c6t6-t27
   epair(kk,1)=epair(kk,1)+e0*eps2(nci(3,k))
   eatom(i1)=eatom(i1)+e0*eps2(nci(3,k))*half
   eatom(i2)=eatom(i2)+e0*eps2(nci(3,k))*half
   e0=q(i1)*q(i2)*oner
   epair(kk,2)=epair(kk,2)+e0*eps1(nci(3,k))
   eatom(i1)=eatom(i1)+e0*eps1(nci(3,k))*half
   eatom(i2)=eatom(i2)+e0*eps1(nci(3,k))*half
   x    =zab(iz1,iz2)  
   alpha=r0ab(iz1,iz2)
   t27  =x*dexp(-alpha*r)
   e0   =t27*oner
   epair(kk,3)=epair(kk,3)+e0*eps2(nci(3,k))
   eatom(i1)=eatom(i1)+e0*eps2(nci(3,k))*half
   eatom(i2)=eatom(i2)+e0*eps2(nci(3,k))*half
end do
!
!     hydrogen and halogen bonds
!

do k=1,nhb   
   i1 =hb(1,k)
   i2 =hb(2,k)
   kk=lina(i1,i2)
   r  =sqrt((xyz(1,i1)-xyz(1,i2))**2 &
     &    +(xyz(2,i1)-xyz(2,i2))**2 &
     &    +(xyz(3,i1)-xyz(3,i2))**2)
   i  =hb(3,k)
   c1 =vhb(1,k)
   c2 =vhb(2,k)
   if (at(i).eq.1) then
      call eabh0(n,i1,i2,i,r,xyz,c1,c2,e0)
   else
      call eabx (n,i1,i2,i,xyz,c1,e0)
   end if
   epair(kk,4)=epair(kk,4)+e0  
   eatom(i1)=eatom(i1)+e0*third
   eatom(i2)=eatom(i2)+e0*third
   eatom(i )=eatom(i )+e0*third
end do
!
!     the same as above for the second partner of pairs
!

do k=1,nnci
   i1=nci(1,k)
   i2=nci(2,k)
   kk=lina(i1,i2)
   dx=xyz2(1,i1)-xyz2(1,i2)
   dy=xyz2(2,i1)-xyz2(2,i2)
   dz=xyz2(3,i1)-xyz2(3,i2)
   r2=dx*dx+dy*dy+dz*dz
   r =sqrt(r2)
   oner =1.0d0/r
   iz1=at(i1)
   iz2=at(i2)
   R0=r0094(iz1,iz2)
   c6=c66ab(i2,i1)
   r4=r2*r2
   r6=r4*r2
   r06=R0**6
   t6=r6+r06
   t8=r6*r2+r06*R0*R0
   c6t6=c6/t6
   c6t8=c6/t8
   t27=sr42(iz1,iz2)*c6t8
   e0=-c6t6-t27
   epair2(kk,1)=epair2(kk,1)+e0*eps2(nci(3,k))
   eatom2(i1)=eatom2(i1)+e0*eps2(nci(3,k))*half
   eatom2(i2)=eatom2(i2)+e0*eps2(nci(3,k))*half
   e0=q(i1)*q(i2)*oner
   epair2(kk,2)=epair2(kk,2)+e0*eps1(nci(3,k))
   eatom2(i1)=eatom2(i1)+e0*eps1(nci(3,k))*half
   eatom2(i2)=eatom2(i2)+e0*eps1(nci(3,k))*half
   x    =zab(iz1,iz2)
   alpha=r0ab(iz1,iz2)
   t27  =x*dexp(-alpha*r)
   e0   =t27*oner
   epair2(kk,3)=epair2(kk,3)+e0*eps2(nci(3,k))
   eatom2(i1)=eatom2(i1)+e0*eps2(nci(3,k))*half
   eatom2(i2)=eatom2(i2)+e0*eps2(nci(3,k))*half
enddo

do k=1,nhb
   i1 =hb(1,k)
   i2 =hb(2,k)
   kk=lina(i1,i2)
   r  =sqrt((xyz2(1,i1)-xyz2(1,i2))**2 &
       &   +(xyz2(2,i1)-xyz2(2,i2))**2 &
       &   +(xyz2(3,i1)-xyz2(3,i2))**2)
   i  =hb(3,k)
   c1 =vhb(1,k)
   c2 =vhb(2,k)
   if(at(i).eq.1)then
      call eabh0(n,i1,i2,i,r,xyz2,c1,c2,e0)
   else
      call eabx (n,i1,i2,i,xyz2,c1,e0)
   endif
   epair2(kk,4)=epair2(kk,4)+e0
   eatom2(i1)=eatom2(i1)+e0*third
   eatom2(i2)=eatom2(i2)+e0*third
   eatom2(i )=eatom2(i )+e0*third
enddo

!
!     call prmat(6,epair,n,0,'Enci+Ehb')
!
epair=epair*627.51
epair2=epair2*627.51

write(10,*) 'non-covalent and HB terms    disp/ES/rep/HB'
k=0
e0=0
e1=0
e2=0
e3=0
e4=0
do i=1,n
   do j=1,i
      k=k+1
      de=sum(epair2(k,1:4))-sum(epair(k,1:4))
      e0=e0+de
      e1=e1+epair2(k,1)-epair(k,1)
      e2=e2+epair2(k,2)-epair(k,2)
      e3=e3+epair2(k,3)-epair(k,3)
      e4=e4+epair2(k,4)-epair(k,4)
      if (abs(de).gt.1) then
         write(10,'(''atom pair:'',2i4, &
           &    '' delta E (kcal)='',F8.2,4x,4F7.2)')i,j,de, &
           &      epair2(k,1)-epair(k,1), &
           &      epair2(k,2)-epair(k,2), &
           &      epair2(k,3)-epair(k,3), &
           &      epair2(k,4)-epair(k,4) 
      end if
   end do
end do
write(10,'(''total disp delta E='',F12.6,F10.3)')e1/627.51,e1
write(10,'(''total ES   delta E='',F12.6,F10.3)')e2/627.51,e2
write(10,'(''total rep  delta E='',F12.6,F10.3)')e3/627.51,e3
write(10,'(''total HB   delta E='',F12.6,F10.3)')e4/627.51,e4
write(10,'(''total NCI  delta E='',F12.6,F10.3)')e0/627.51,e0
esum=esum+e0
      
!
!     strech (bondlenghts)
!
 
do m=1,nbond
   i=bond(1,m)
   j=bond(2,m)
   r2=   ((xyz(1,i)-xyz(1,j))**2 &
     &   +(xyz(2,i)-xyz(2,j))**2 &
     &   +(xyz(3,i)-xyz(3,j))**2)
   r =sqrt(r2)
   rij=vbond(1,m)
   kij=vbond(2,m)
   aai=vbond(3,m)
   aai2 =aai/2
   epart(m,1)=kij*(1.+(rij/r)**aai - 2.*(rij/r)**aai2)
   eatom(i)=eatom(i)+epart(m,1)*half
   eatom(j)=eatom(j)+epart(m,1)*half
end do
!
!     bends (1-2-3 angles)    
!  

do m=1,nangl
   j = angl(1,m)
   i = angl(2,m)
   k = angl(3,m)
   c0  =vangl(1,m)
   kijk=vangl(2,m)
   va(1:3) = xyz(1:3,i)
   vb(1:3) = xyz(1:3,j)
   vc(1:3) = xyz(1:3,k)
   call vsub(va,vb,vab,3)
   call vsub(vc,vb,vcb,3)
   rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
   rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
   call crprod(vcb,vab,vp)
   rp = vlen(vp)+1.d-14
   call impsc(vab,vcb,cosa)
   cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
   theta= dacos(cosa)
   call abdamp(at(i),at(j),rab2,dampij,damp2ij)
   call abdamp(at(k),at(j),rcb2,dampjk,damp2jk)
   damp=dampij*dampjk
   if(pi-c0.lt.0.1)then
      dt  = theta - c0  
      ea  = kijk * dt**2
   else
      ea=kijk*(cos(theta)-cos(c0))**2 
   end if
   epart(m,2) = ea * damp
   eatom(i)=eatom(i)+epart(m,2)*third
   eatom(j)=eatom(j)+epart(m,2)*third
   eatom(k)=eatom(k)+epart(m,2)*third
   pbend(m)=theta
end do
!
!     torsional angles
!   
do m=1,ntors
   i=tors(1,m)
   j=tors(2,m)
   k=tors(3,m)
   l=tors(4,m)
   nt=tors(5,m)
   phi0 =vtors(1,m)
   if (tors(6,m).ne.2) then
      vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
      damp=dampjk*dampij*dampkl
      phi=valijkl(n,xyz,i,j,k,l)
      et=0
      mm=3
      do it=1,nt
         rn=vtors(mm,m)
         dphi1=phi-phi0
         dphi2=phi+phi0-pi2
         c1=rn*dphi1+vtors(mm+1,m)
         c2=rn*dphi2+vtors(mm+1,m)
         x1cos=cos(c1)
         x2cos=cos(c2)
         phipi=phi-pi
         ef   =erf(phipi)
         e1   =vtors(mm+2,m)*(1.+x1cos)
         e2   =vtors(mm+2,m)*(1.+x2cos)
         et   =et+0.5*(1.-ef)*e1+(0.5+0.5*ef)*e2
         mm=mm+3
      end do
      epart(m,3)=et*vtors(2,m)*damp
      eatom(i)=eatom(i)+epart(m,3)*quart
      eatom(j)=eatom(j)+epart(m,3)*quart
      eatom(k)=eatom(k)+epart(m,3)*quart
      eatom(l)=eatom(l)+epart(m,3)*quart
   else
      vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
      damp=dampjk*dampij*dampkl
      phi=omega(n,xyz,i,j,k,l)
      rn=vtors(3,m)
      if (rn.gt.1.d-6) then
         dphi1=phi-phi0
         c1=dphi1+pi              
         x1cos=cos(c1)
         et=(1.+x1cos)*vtors(2,m)
      else
         et=vtors(2,m)*(cos(phi)-cos(phi0))**2 
      end if
      epart(m,3)=et*damp
      eatom(i)=eatom(i)+epart(m,3)*quart
      eatom(j)=eatom(j)+epart(m,3)*quart
      eatom(k)=eatom(k)+epart(m,3)*quart
      eatom(l)=eatom(l)+epart(m,3)*quart
   end if
   ptors(m)=phi
end do
!
!     do the same again for the bond and angle terms!
!
do m=1,nbond
   i=bond(1,m)
   j=bond(2,m)
   r2=   ((xyz2(1,i)-xyz2(1,j))**2 &
     &   +(xyz2(2,i)-xyz2(2,j))**2 &
     &   +(xyz2(3,i)-xyz2(3,j))**2)
   r =sqrt(r2)
   rij=vbond(1,m)
   kij=vbond(2,m)
   aai=vbond(3,m)
   aai2 =aai/2
   epart2(m,1)=kij*(1.+(rij/r)**aai - 2.*(rij/r)**aai2)
   eatom2(i)=eatom2(i)+epart2(m,1)*half
   eatom2(j)=eatom2(j)+epart2(m,1)*half
end do

do m=1,nangl
   j = angl(1,m)
   i = angl(2,m)
   k = angl(3,m)
   c0  =vangl(1,m)
   kijk=vangl(2,m)
   va(1:3) = xyz2(1:3,i)
   vb(1:3) = xyz2(1:3,j)
   vc(1:3) = xyz2(1:3,k)
   call vsub(va,vb,vab,3)
   call vsub(vc,vb,vcb,3)
   rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
   rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
   call crprod(vcb,vab,vp)
   rp = vlen(vp)+1.d-14
   call impsc(vab,vcb,cosa)
   cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
   theta= dacos(cosa)
   call abdamp(at(i),at(j),rab2,dampij,damp2ij)
   call abdamp(at(k),at(j),rcb2,dampjk,damp2jk)
   damp=dampij*dampjk
   if (pi-c0.lt.0.1) then
      dt  = theta - c0  
      ea  = kijk * dt**2
   else
      ea=kijk*(cos(theta)-cos(c0))**2 
   end if
   epart2(m,2) = ea * damp
   eatom2(i)=eatom2(i)+epart2(m,2)*third
   eatom2(j)=eatom2(j)+epart2(m,2)*third
   eatom2(k)=eatom2(k)+epart2(m,2)*third
   pbend2(m)=theta
end do

do m=1,ntors
   i=tors(1,m)
   j=tors(2,m)
   k=tors(3,m)
   l=tors(4,m)
   nt=tors(5,m)
   phi0 =vtors(1,m)
   if (tors(6,m).ne.2) then
      vab(1:3) = xyz2(1:3,i)-xyz2(1:3,j)
      vcb(1:3) = xyz2(1:3,j)-xyz2(1:3,k)
      vdc(1:3) = xyz2(1:3,k)-xyz2(1:3,l)
      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
      damp=dampjk*dampij*dampkl
      phi=valijkl(n,xyz2,i,j,k,l)
      et=0
      mm=3
      do it=1,nt
         rn=vtors(mm,m)
         dphi1=phi-phi0
         dphi2=phi+phi0-pi2
         c1=rn*dphi1+vtors(mm+1,m)
         c2=rn*dphi2+vtors(mm+1,m)
         x1cos=cos(c1)
         x2cos=cos(c2)
         phipi=phi-pi
         ef   =erf(phipi)
         e1   =vtors(mm+2,m)*(1.+x1cos)
         e2   =vtors(mm+2,m)*(1.+x2cos)
        et   =et+0.5*(1.-ef)*e1+(0.5+0.5*ef)*e2
         mm =mm+3
      end do
      epart2(m,3)=et*vtors(2,m)*damp
      eatom2(i)=eatom2(i)+epart2(m,3)*quart
      eatom2(j)=eatom2(j)+epart2(m,3)*quart
      eatom2(k)=eatom2(k)+epart2(m,3)*quart
      eatom2(l)=eatom2(l)+epart2(m,3)*quart
   else
      vab(1:3) = xyz2(1:3,i)-xyz2(1:3,j)
      vcb(1:3) = xyz2(1:3,j)-xyz2(1:3,k)
      vdc(1:3) = xyz2(1:3,k)-xyz2(1:3,l)
      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
      damp=dampjk*dampij*dampkl
      phi=omega(n,xyz2,i,j,k,l)
      rn=vtors(3,m)
      if(rn.gt.1.d-6)then
         dphi1=phi-phi0
         c1=dphi1+pi              
         x1cos=cos(c1)
         et=(1.+x1cos)*vtors(2,m)
      else
         et=vtors(2,m)*(cos(phi)-cos(phi0))**2 
      end if
      epart2(m,3)=et*damp
      eatom2(i)=eatom2(i)+epart2(m,3)*quart
      eatom2(j)=eatom2(j)+epart2(m,3)*quart
      eatom2(k)=eatom2(k)+epart2(m,3)*quart
      eatom2(l)=eatom2(l)+epart2(m,3)*quart
   end if
   ptors2(m)=phi
end do
!
!     write out the results for the full QMDFF
!
epart=epart*627.51
epart2=epart2*627.51

write(*,*) 'bond terms'
e0=0
do m=1,nbond
   i=bond(1,m)
   j=bond(2,m)
   de=epart2(m,1)-epart(m,1)
   e0=e0+de
   write(10,'(''bond :'',2i4, &
        &    '' delta E (kcal)='',F10.3)')i,j,de
end do
write(10,'(''total bond delta E='',F10.3)')e0
esum=esum+e0

pbend=pbend*180./pi
pbend2=pbend2*180./pi
write(*,*) 'angle terms'
e0=0
do m=1,nangl
   j = angl(1,m)
   i = angl(2,m)
   k = angl(3,m)
   de=epart2(m,2)-epart(m,2)
   e0=e0+de
   write(10,'(''angle :'',3i4, &
       &     '' delta E (kcal)='',F10.3,3F7.2)') &
       & j,i,k,de,vangl(1,m)*180./pi,pbend(m),pbend2(m)
end do
write(10,'(''total bend delta E='',F10.3)')e0
esum=esum+e0

write(10,*) 'torsion terms'
ptors=ptors*180./pi
ptors2=ptors2*180./pi
e0=0
do m=1,ntors
   i=tors(1,m)
   j=tors(2,m)
   k=tors(3,m)
   l=tors(4,m)
   de=epart2(m,3)-epart(m,3)
   e0=e0+de
   write(10,'(i3,'' torsion :'',6i4, &
    &   '' delta E (kcal)='',F10.4,'' angles:'',3F7.2)')m, &
    &    i,j,k,l,tors(5,m),idint(vtors(3,m)),de, & 
    &          vtors(1,m)*180./pi,ptors(m),ptors2(m)
end do
write(10,'(''total torsion delta E='',F10.3)')e0
esum=esum+e0

write(10,'(''total delta E='',F10.3)')esum

write(10,*) 'atomic strain'
write(10,*) 'atom   coord1  coord2    delta'
do i=1,n
   write(10,'(i4,3F9.2)')  &
   &  i,eatom(i)*627.51,eatom2(i)*627.51,(eatom2(i)-eatom(i))*627.51
end do

return
end subroutine ff_anal

