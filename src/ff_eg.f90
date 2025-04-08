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
!     subroutine ff_eg: energies and gradients of the bonded 
!            potential terms    
!
!       Virial tensor components are calculated as explained in:
!       J. Comput. Chem. 23; 667-672, 2002 (Computation of Pressure Components
!       due to Class II Force Fields)
!
!     part of QMDFF
!


subroutine ff_eg(n,at,xyz,e,g)
use qmdff
use debug ! additional printout if desired
use pbc_mod
implicit none  
integer::n,at(n)
real(kind=8)::xyz(3,n),e,g(3,n)
real(kind=8)::e_tmp
real(kind=8)::vp(3),dt,deddt,rp,cosa,rab2,rcb2,vlen,rmul1,rmul2
real(kind=8)::cosna(4),sinna(4),v(4),sa(4),ca(4),phix(4),dphi(4)
real(kind=8)::va(3),vb(3),vc(3),vd(3),vba(3),vcb(3),vdc(3),vab(3)
real(kind=8)::vca(3),vdb(3),vt(3),vu(3),vtu(3),dedt(3),dedu(3)
real(kind=8)::deda(3),dedb(3),dedc(3),dedd(3),vtmp1(3),vtmp2(3)
real(kind=8)::rt2,ru2,rtu,rcb,rtru,damp,damp2,rac2,racut,omega
real(kind=8)::aai,aai2,r2,rjk,dampjk,damp2jk,dampij,damp2ij,dampjl
real(kind=8)::rkl,rjl,dampkl,damp2kl,damp2jl,term1(3),term2(3),term3(3)
integer::ib,id
real(kind=8)::dda(3),ddb(3),ddc(3),ddd(3)
real(kind=8)::rbond(3)
!   The newly added gradient components
real(kind=8)::g_local_a(3),g_local_b(3),g_local_c(3),g_local_d(3) 

real(kind=8)::r,thab,thbc,ra(3),fac,eb,or,kij,rij,alpha,dij
real(kind=8)::dei(3),dej(3),dek(3),kijk,rb(3),ea,dedtheta,ban
real(kind=8)::c0,c1,c2,theta,sinth,costh,expo,valijkl
real(kind=8)::dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),step
real(kind=8)::et,phi0,rn,kijkl,dedphi,phi,cosrph,cosrph0,er,el
real(kind=8)::dphi1,dphi2,x1cos,x2cos,x1sin,x2sin,phipi,e1,e2,ef
integer::i,j,k,l,m,ic,ia,list(4),jj,mm,ii,kk,it,nt,it2

e=0
g=0
!
!     If the eval_cutoff option is activated, reset the within_cut array
!     and the e_partition array with the energy for each atom
!  
if (eval_cutoff) then
   within_cut=.false.
   e_partition=0.d0
end if
!
!     calculate energy and gradients: For all coordinate types, loop over 
!     all included coordinates and sum up the parts
!
!     strech terms (bondlengths)
!
do m=1,nbond
   i=bond(1,m)
   j=bond(2,m)
!
!     calculate bondlength
!
   do ic=1,3
      rbond(ic)=xyz(ic,i)-xyz(ic,j)
   end do
 
!
!     apply periodic boundaries, if needed
!
   if (periodic) then
      call box_image(rbond)
   end if
   r2=rbond(1)**2+rbond(2)**2+rbond(3)**2
   r =sqrt(r2)

   rij=vbond(1,m)
   kij=vbond(2,m)
   aai=vbond(3,m)
   aai2 =aai/2
!
!     the energy
!
   e_tmp=kij*(1.+(rij/r)**aai-2.*(rij/r)**aai2)
   e=e+e_tmp
!
!     For energy partition in eval_cutoff: add half energy to each atom
!
   if (eval_cutoff) then
      e_partition(i)=e_partition(i)+e_tmp/2.d0
      e_partition(j)=e_partition(j)+e_tmp/2.d0
   end if
   if (eval_cutoff) then
      within_cut(i,j) = .true.
      within_cut(j,i) = .true.
   end if

!
!     the gradient
!
   fac=aai*kij*(-(rij/r)**aai+(rij/r)**aai2)/r2
   do ic=1,3
      g_local_a(ic)=fac*rbond(ic)
   end do

   g(:,i)=g(:,i)+g_local_a
   g(:,j)=g(:,j)-g_local_a
!
!     the components of the virial tensor, if needed
!
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+rbond(1)*g_local_a(1)
      vir_ten(2,1)=vir_ten(2,1)+rbond(2)*g_local_a(1)
      vir_ten(3,1)=vir_ten(3,1)+rbond(3)*g_local_a(1)
      vir_ten(1,2)=vir_ten(1,2)+rbond(1)*g_local_a(2)
      vir_ten(2,2)=vir_ten(2,2)+rbond(2)*g_local_a(2)
      vir_ten(3,2)=vir_ten(3,2)+rbond(3)*g_local_a(2)
      vir_ten(1,3)=vir_ten(1,3)+rbond(1)*g_local_a(3)
      vir_ten(2,3)=vir_ten(2,3)+rbond(2)*g_local_a(3)
      vir_ten(3,3)=vir_ten(3,3)+rbond(3)*g_local_a(3)
   end if
!
!     If debugging is activated, store the component in the array
!
   if (do_debug) then
      comp_no=m
      act_part=kij*(1.+(rij/r)**aai-2.*(rij/r)**aai2)
      qmdff_parts(1,struc_no,comp_no)=act_part
      write(288,*) & 
       &       "bond",i,j," energy",kij*(1.+(rij/r)**aai-2.*(rij/r)**aai2)," eq. length:",r," act. length:",rij, &
       & " f.c.",kij," exp.a.",aai," dr/dxi:", g(1,i)," dr/dyi:", g(2,i)," dr/dzi:", g(3,i)," dr/dxj:", &
       &   g(1,j)," dr/dyj:", g(2,j)," dr/dzj:",g(3,j)
   end if
end do
!
!     bend terms (1-2-3 angles)   
!   
do m=1,nangl
   j=angl(1,m)
   i=angl(2,m)
   k=angl(3,m)
   c0=vangl(1,m)
   kijk=vangl(2,m)
   va(1:3)=xyz(1:3,i)
   vb(1:3)=xyz(1:3,j)
   vc(1:3)=xyz(1:3,k)
   
   vab=va-vb
   vcb=vc-vb
!
!     apply periodic boundaries, if needed
!
   if (periodic) then
      call box_image(vab)
      call box_image(vcb)
   end if


! 
!     vectors of both included bonds
!
   rab2=vab(1)*vab(1)+vab(2)*vab(2)+vab(3)*vab(3)
   rcb2=vcb(1)*vcb(1)+vcb(2)*vcb(2)+vcb(3)*vcb(3)
!
!     calculate bond angle
!
   call crprod(vcb,vab,vp)
   rp=vlen(vp)+1.d-14
   call impsc(vab,vcb,cosa)
   cosa=dble(min(1.0d0,max(-1.0d0,cosa)))
   theta=dacos(cosa)
!
!     damping for large distances
!
   call abdamp(at(i),at(j),rab2,dampij,damp2ij)
   call abdamp(at(k),at(j),rcb2,dampjk,damp2jk)
   damp=dampij*dampjk
!
!     if the angle is almost linear
!
   if (pi-c0.lt.1.d-6) then
      dt=theta-c0  
      ea=kijk*dt**2
      deddt=2.d0*kijk*dt 
   else
      ea=kijk*(cosa-cos(c0))**2 
      deddt=2.*kijk*sin(theta)*(cos(c0)-cosa)               
  
    end if
!
!     add the calculated energy
!
!     --> for debugging
!
!   if (do_debug) then
!      comp_no=comp_no+1
!      act_part=ea*damp
!      qmdff_parts(1,struc_no,comp_no)=act_part
!      write(parts_labels(1,comp_no),'(a,i4,i4,i4)') "angle",i,j,k
!   end if
   e=e+ea*damp
!
!     For energy partition in eval_cutoff: add third energy to each atom
!
   if (eval_cutoff) then
      e_partition(i)=e_partition(i)+ea*damp/3.d0
      e_partition(j)=e_partition(j)+ea*damp/3.d0
      e_partition(k)=e_partition(k)+ea*damp/3.d0
   end if
   if (eval_cutoff) then
      within_cut(i,j) = .true.
      within_cut(i,k) = .true.
      within_cut(j,i) = .true.
      within_cut(k,i) = .true.
      within_cut(j,k) = .true.
      within_cut(k,j) = .true.
   end if

!
!     calculate the gradient
!
   call crprod(vab,vp,deda)
   rmul1 = -deddt / (rab2*rp)
   deda=deda*rmul1
   call crprod(vcb,vp,dedc)
   rmul2 =  deddt / (rcb2*rp)
   dedc=dedc*rmul2
   dedb=deda+dedc   ! sum of a and c!
   term1(1:3)=ea*damp2ij*dampjk*vab(1:3)
   term2(1:3)=ea*damp2jk*dampij*vcb(1:3)
   
   g_local_a = deda(1:3)*damp+term1(1:3)
   g_local_b = -dedb(1:3)*damp-term1(1:3)-term2(1:3)
   g_local_c = dedc(1:3)*damp+term2(1:3)

   g(1:3,i)=g(1:3,i)+g_local_a
   g(1:3,j)=g(1:3,j)+g_local_b
   g(1:3,k)=g(1:3,k)+g_local_c
!
!     Calculate the virial tensor components, if needed!   
!
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)+vcb(1)*g_local_c(1)
      vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)+vcb(2)*g_local_c(1)
      vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)+vcb(3)*g_local_c(1)
      vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)+vcb(1)*g_local_c(2)
      vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)+vcb(2)*g_local_c(2)
      vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)+vcb(3)*g_local_c(2)
      vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)+vcb(1)*g_local_c(3)
      vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)+vcb(2)*g_local_c(3)
      vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)+vcb(3)*g_local_c(3)
   end if
!
!     If debugging is activated, store the component in the array
!

   if (do_debug) then
      comp_no=m
      act_part=kij*(1.+(rij/r)**aai-2.*(rij/r)**aai2)
      qmdff_parts(1,struc_no,comp_no)=act_part
      write(289,*) &
       &       "angle",i,j,k, " energy",ea*damp," act. angle:",c0, " eq. angle:",acos(cosa), &
       & " f.c.",kijk," damping.",damp," da/dxi:", g(1,i)," da/dyi:", g(2,i)," da/dzi:", g(3,i)," da/dxj:", &
       &   g(1,j)," da/dyj:", g(2,j)," da/dzj:",g(3,j), "da/dxk:", &
       &   g(1,k)," da/dyk:", g(2,k)," da/dzk:",g(3,k)
   end if



end do
!
!     torsion ter,s  
!      
do m=1,ntors
!
!     the four involved atoms
!
   i=tors(1,m)
   j=tors(2,m)
   k=tors(3,m)
   l=tors(4,m)
   nt=tors(5,m)
   phi0 =vtors(1,m)
!----------------------------------------------------
!     the normal torsions
!----------------------------------------------------
   if(tors(6,m).ne.2)then
      vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
!
!     apply periodic boundaries, if needed
!
      if (periodic) then
         call box_image(vab)
         call box_image(vcb)
         call box_image(vdc)
      end if

!
!     three bondlengths are needed
!
      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
!
!     calculate all damping functions
!
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
      damp= dampjk*dampij*dampkl
!
!     the torsional angle phi
!
      phi=valijkl(n,xyz,i,j,k,l)
      call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
      et=0
      dij=0
      mm=3!
!
!     loop over the number of torsional minima along curve
!
      if (do_debug) then
         write(290,'(a,i4,i4,i4,i4)') "torsion",i,j,k,l
      end if
      do it=1,nt
         rn=vtors(mm,m)
         dphi1=phi-phi0
         dphi2=phi+phi0-pi2
         c1=rn*dphi1+vtors(mm+1,m)
         c2=rn*dphi2+vtors(mm+1,m)
         x1cos=cos(c1)
         x2cos=cos(c2)
         x1sin=sin(c1)
         x2sin=sin(c2)
         phipi=phi-pi
!
!     deviation of the angle from linearity (=pi)
!
         ef=erf(phipi)
         e1=vtors(mm+2,m)*(1.+x1cos)
         e2=vtors(mm+2,m)*(1.+x2cos)
         et=et+0.5*(1.-ef)*e1+(0.5+0.5*ef)*e2
         expo=exp(-phipi**2)/spi
         dij=dij-expo*e1- 0.5*(1.-ef)*vtors(mm+2,m)*x1sin*rn+ &
           &   expo*e2-(0.5+0.5*ef)*vtors(mm+2,m)*x2sin*rn    
         if (do_debug) then
       !      write(*,*) it,mm,e1,e2,phi,phi0,dampjk,dampij,dampkl,damp,c1,c2,vtors(2,m)
            write(290,*) "Minimum No.",it,"of",nt,"damp(ij)",dampij,"damp(jk)",dampjk,&
                       & "damp(kl)",dampkl,"damp_tot",damp,"eq. angle",phi0,"act. angle",phi,"force const.",vtors(2,m)

         end if
         mm=mm+3
      end do
!
!     sum up the torsional energy
!
!     --> for debugging
!
!   if (do_debug) then
!      comp_no=comp_no+1
!      act_part=ea*damp
!      qmdff_parts(1,struc_no,comp_no)=act_part
!      write(parts_labels(1,comp_no),'(a,i4,i4,i4,i4)') "torsion",i,j,k,l
!   end if

      et=et*vtors(2,m)

!
!     the torsional gradient!
!
      dij=dij*vtors(2,m)*damp

      term1(1:3)=et*damp2ij*dampjk*dampkl*vab(1:3)
      term2(1:3)=et*damp2jk*dampij*dampkl*vcb(1:3)
      term3(1:3)=et*damp2kl*dampij*dampjk*vdc(1:3)

      g_local_a = dij*dda(1:3)+term1
      g_local_b = dij*ddb(1:3)-term1+term2
      g_local_c = dij*ddc(1:3)+term3-term2
      g_local_d = dij*ddd(1:3)-term3
  
      g(1:3,i)=g(1:3,i)+g_local_a
      g(1:3,j)=g(1:3,j)+g_local_b
      g(1:3,k)=g(1:3,k)+g_local_c
      g(1:3,l)=g(1:3,l)+g_local_d

!
!     Calculate the virial tensor components, if needed!   
!

      if (calc_vir) then
         vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)+vcb(1)*(g_local_c(1)+ &
                     & g_local_d(1))+vdc(1)*g_local_d(1)
         vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)+vcb(2)*(g_local_c(1)+ &
                     & g_local_d(1))+vdc(2)*g_local_d(1)
         vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)+vcb(3)*(g_local_c(1)+ &
                     & g_local_d(1))+vdc(3)*g_local_d(1)
         vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)+vcb(1)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(1)*g_local_d(2)
         vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)+vcb(2)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(2)*g_local_d(2)
         vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)+vcb(3)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(3)*g_local_d(2)
         vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)+vcb(1)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(1)*g_local_d(3)
         vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)+vcb(2)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(2)*g_local_d(3)
         vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)+vcb(3)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(3)*g_local_d(3)
      end if

!
!     If debugging is activated, store the component in the array
!

      if (do_debug) then
         comp_no=m
         write(290,*) &
          &       "total energy",et, &
          & " dt/dxi:", g(1,i)," dt/dyi:", g(2,i)," dt/dzi:", g(3,i)," dt/dxj:", &
          &   g(1,j)," dt/dyj:", g(2,j)," dt/dzj:",g(3,j),"dt/dxk:", &
          &   g(1,k)," dt/dyk:", g(2,k)," dt/dzk:",g(3,k),"dt/dxl:", &
          &   g(1,l)," dt/dyl:", g(2,l)," dt/dzl:",g(3,l)
      end if



!
!     scale energy by damping factor
!
!     --> for debugging
!
!      if (do_debug) then
!         comp_no=nbond+nangl+m
!         act_part=et*damp
!         qmdff_parts(1,struc_no,comp_no)=act_part
!      end if
      e=e+et*damp
!
!     For energy partition in eval_cutoff: add quarter energy to each atom
!
      if (eval_cutoff) then
         e_partition(i)=e_partition(i)+et*damp/4.d0
         e_partition(j)=e_partition(j)+et*damp/4.d0
         e_partition(k)=e_partition(k)+et*damp/4.d0
         e_partition(l)=e_partition(l)+et*damp/4.d0
      end if
      if (eval_cutoff) then
         within_cut(i,j) = .true.
         within_cut(i,k) = .true.
         within_cut(i,l) = .true.
         within_cut(j,i) = .true.
         within_cut(k,i) = .true.
         within_cut(l,i) = .true.
         within_cut(j,k) = .true.
         within_cut(j,l) = .true.
         within_cut(k,j) = .true.
         within_cut(l,j) = .true.
         within_cut(k,l) = .true.
         within_cut(l,k) = .true.
      end if



!----------------------------------------------------    
!     torsions for inversions!
!----------------------------------------------------
   else
!
!     do basically the same as above...
!
      vab(1:3) = xyz(1:3,j)-xyz(1:3,i)
      vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
      vdc(1:3) = xyz(1:3,j)-xyz(1:3,l)
!
!     apply periodic boundaries, if needed
!
      if (periodic) then
         call box_image(vab)
         call box_image(vcb)
         call box_image(vdc)
      end if

      rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      rjl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
      call abdamp(at(i),at(j),rij,dampij,damp2ij)
      call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
      call abdamp(at(j),at(l),rjl,dampjl,damp2jl)
      damp= dampjk*dampij*dampjl

      phi=omega(n,xyz,i,j,k,l)
      call domegadr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)
       
      rn=vtors(3,m)
      if (rn.gt.1.d-6)then
         dphi1=phi-phi0
         c1=dphi1+pi              
         x1cos=cos(c1)
         x1sin=sin(c1)
         et   =(1.+x1cos)*vtors(2,m)
         dij  =-x1sin*vtors(2,m)*damp
      else
         et =   vtors(2,m)*(cos(phi) -cos(phi0))**2 
         dij=2.*vtors(2,m)* sin(phi)*(cos(phi0)-cos(phi))*damp
      end if

      term1(1:3)=et*damp2ij*dampjk*dampjl*vab(1:3)
      term2(1:3)=et*damp2jk*dampij*dampjl*vcb(1:3)
      term3(1:3)=et*damp2jl*dampij*dampjk*vdc(1:3)

      g_local_a = dij*dda(1:3)-term1
      g_local_b = dij*ddb(1:3)+term1+term2+term3
      g_local_c = dij*ddc(1:3)-term2
      g_local_d = dij*ddd(1:3)-term3


      g(1:3,i)=g(1:3,i)+g_local_a
      g(1:3,j)=g(1:3,j)+g_local_b
      g(1:3,k)=g(1:3,k)+g_local_c
      g(1:3,l)=g(1:3,l)+g_local_d
      e=e+et*damp

!
!     For energy partition in eval_cutoff: add quarter energy to each atom
!
      if (eval_cutoff) then
         e_partition(i)=e_partition(i)+et*damp/4.d0
         e_partition(j)=e_partition(j)+et*damp/4.d0
         e_partition(k)=e_partition(k)+et*damp/4.d0
         e_partition(l)=e_partition(l)+et*damp/4.d0
      end if
      if (eval_cutoff) then
         within_cut(i,j) = .true.
         within_cut(i,k) = .true.
         within_cut(i,l) = .true.
         within_cut(j,i) = .true.
         within_cut(k,i) = .true.
         within_cut(l,i) = .true.
         within_cut(j,k) = .true.
         within_cut(j,l) = .true.
         within_cut(k,j) = .true.
         within_cut(l,j) = .true.
         within_cut(k,l) = .true.
         within_cut(l,k) = .true.
      end if


!
!     Calculate the virial tensor components, if needed!   
!     a=i   b=j   c=k  d=l
!     --> Maybe the out of plane variant needed with b as central atom??
!
      if (calc_vir) then
         vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)+vcb(1)*(g_local_c(1)+ &
                     & g_local_d(1))+vdc(1)*g_local_d(1)
         vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)+vcb(2)*(g_local_c(1)+ &
                    & g_local_d(1))+vdc(2)*g_local_d(1)
         vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)+vcb(3)*(g_local_c(1)+ &
                     & g_local_d(1))+vdc(3)*g_local_d(1)
         vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)+vcb(1)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(1)*g_local_d(2)
         vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)+vcb(2)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(2)*g_local_d(2)
         vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)+vcb(3)*(g_local_c(2)+ &
                     & g_local_d(2))+vdc(3)*g_local_d(2)
         vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)+vcb(1)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(1)*g_local_d(3)
         vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)+vcb(2)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(2)*g_local_d(3)
         vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)+vcb(3)*(g_local_c(3)+ &
                     & g_local_d(3))+vdc(3)*g_local_d(3)
      end if



!
!     --> for debugging
!
!      if (do_debug) then
!         comp_no=comp_no+1
!         act_part=et*damp
!         qmdff_parts(1,struc_no,comp_no)=act_part
!         write(parts_labels(1,comp_no),'(a,i4,i4,i4,i4)') "inversopm",i,j,k,l
!      end if

   end if
end do
!
!     here restraints can be added! 
!     maybe some metal surface could be simulated using restraints??? 
!
!      call egrestrain(n,xyz,g,er)
!      e = e + er
!

return
end subroutine ff_eg
