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
!     subroutine ff_e: calculate only energy of the qmdff for the reference
!                 (used for dissociation energy calculation)
!
!     part of QMDFF
!

subroutine ff_e(n,at,xyz,q,r0ab,zab,r0094,sr42,c66ab,e)
use qmdff
implicit none  
integer::n,at(n)
real(kind=8)::xyz(3,n),q(n),e
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),c66ab(n,n)
integer::i1,i2,iz1,iz2,k,lina,i,j,m,ic,nt,it,l,it2,mm
real(kind=8)::r2,r,r4,r6,r06,R0,t6,t8,c6t6,c6t8,t27,drij,c1,c2,aai,aai2
real(kind=8)::dx,dy,dz,c6,x,alpha,oner,aa,damp,damp2,e0,de,rij,kij,et
real(kind=8)::c0,kijk,va(3),vb(3),vc(3),dt,phi,phi0,rp,rkl,rjk,rjl
real(kind=8)::ea,cosa,rab2,rcb2,vlen,theta,vdc(3),vcb(3),vab(3),vp(3)
real(kind=8)::dampkl,damp2kl,damp2jl,dampjk,damp2jk,dampij,damp2ij,dampjl
real(kind=8)::valijkl,esum,rn,dphi,e1,e2,phipi,ef,x1cos,x2cos,ban
real(kind=8)::dphi1,dphi2,omega

e=0
!
!     non covalent interactions
!
do k=1,nnci               
   i1=nci(1,k)
   i2=nci(2,k)
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
   e=e+e0*eps2(nci(3,k))
   e0=q(i1)*q(i2)*oner
   e=e+e0*eps1(nci(3,k))
   x    =zab(iz1,iz2)  
   alpha=r0ab(iz1,iz2)
   t27  =x*dexp(-alpha*r)
   e0   =t27*oner
   e=e+e0*eps2(nci(3,k))
end do
!
!     hydrogen and halogen bonds
!
do k=1,nhb   
   i1 =hb(1,k)
   i2 =hb(2,k)
   r  =sqrt((xyz(1,i1)-xyz(1,i2))**2 &
      &    +(xyz(2,i1)-xyz(2,i2))**2 &
      &    +(xyz(3,i1)-xyz(3,i2))**2)
   i  =hb(3,k)
   c1 =vhb(1,k)
   c2 =vhb(2,k)
   call eabh0(n,i1,i2,i,r,xyz,c1,c2,e0)
   e=e+e0  
end do
!
!     strech (bondlenghts)
!
do m=1,nbond
   i=bond(1,m)
   j=bond(2,m)
   r2=((xyz(1,i)-xyz(1,j))**2 &
     &   +(xyz(2,i)-xyz(2,j))**2 &
     &   +(xyz(3,i)-xyz(3,j))**2)
   r =sqrt(r2)
   rij=vbond(1,m)
   kij=vbond(2,m)
   aai=vbond(3,m)
   aai2 =aai/2
   e=e+kij*(1.+(rij/r)**aai - 2.*(rij/r)**aai2)
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
   e=e+ea*damp
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
!
!     normal torsions
!
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
      e=e+et*vtors(2,m)*damp
!
!     torsional inversions
!
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
      e=e+et*damp
   end if
end do

return
end subroutine ff_e
