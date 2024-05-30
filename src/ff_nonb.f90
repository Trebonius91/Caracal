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
!     subroutine ff_nonb: calculate the nonbonded part of energy and gradient of the FF
!
!     part of QMDFF
!
subroutine ff_nonb(n,at,xyz,q,r0ab,zab,r0094,sr42,c66ab,enb,g)
use qmdff
use debug
use pbc_mod
implicit none
integer n,at(*)
real(kind=8)::xyz(3,n),enb,g(3,n),q(n)
real(kind=8)::r0ab(94,94),zab(94,94),r0094(94,94),sr42(94,94),c66ab(n,n)

integer::i1,i2,iz1,iz2,k,nk,i,j
real(kind=8)::r2,r,r4,r6,r06,R0,t6,t8,c6t6,c6t8,t27,drij,e
real(kind=8)::dx,dy,dz,c6,x,alpha,oner,aa,damp,damp2,e0
!  for smooth cutoff
real(kind=8)::switch,exp_switch
!  for the Ewald summation
real(kind=8)::k_ew_cut  ! adjustable scaling parameters
real(kind=8)::f12,erfterm
real(kind=8)::medium
real(kind=8)::de,dedx,dedy,dedz
real(kind=8)::r_ewald
real(kind=8)::e_ew  ! the Ewald sum Coulomb energy
real(kind=8)::de_ew(3,n)  ! the Ewald sum Coulomb gradient
real(kind=8)::e_rec ! the reciprocal Ewald Coulomb energy
real(kind=8)::de_rec(3,n)  ! the reciprocal Ewald sum Coulomb gradient
real(kind=8)::xd,yd,zd  ! cell dipole correction
real(kind=8)::e_local,r_buf,r_buf2
real(kind=8)::fs_par,term
real(kind=8)::bohr
real(kind=8)::e_brute ! energy of the brute force Ewald method
integer::nbox,xbox,ybox,zbox  ! for brute force Ewald
integer::nbox_sq,rbox,rbox_sq
real(kind=8),allocatable::pot_shell(:)  ! energy for shells
real(kind=8)::rbox_vec(3)
!   Dipole moment calculation
real(kind=8)::dipolemom(3),charge_com(3)
!   The newly added gradient components
real(kind=8)::g_local_a(3),g_local_b(3)
real(kind=8)::vab(3)  ! for periodic image calculation

bohr=0.52917721092d0
e=0.d0
!
!     if no nonbonded terms are existent
!
if (nnci.le.1 .and. nmols.eq.0) return
!
!     For nonperiodic systems, set the VDW cutoff always to 10 Angstroms
!
if (.not. periodic) then
   vdw_cut=10.d0/bohr
end if 
!
!     loop over dispersion/Lennard-Jones energy terms
!
do k=1,nnci
   i1=nci(1,k)
   i2=nci(2,k)
   nk=nci(3,k)
   vab=xyz(1:3,i1)-xyz(1:3,i2)
!
!     apply periodic boundaries, if needed
!
   if (periodic) then 
      call box_image(vab(1),vab(2),vab(3))
   end if
 
   r2=dot_product(vab,vab)

   r =sqrt(r2)

!
!     only proceed if the distance is below the cutoff
!
   if (periodic .and. r .gt. vdw_cut) then 
      cycle
   end if 
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
   e0=c6t6+t27
   e=e-e0*eps2(nk)
!   write(*,*) "e",e,r,k,i1,i2,e0,eps2(nk),nk

   drij=eps2(nk)*(c6t6*6.0d0*r4/t6+8.0d0*t27*r6/t8)
   
!
!     determine gradient vector entries
!
   g_local_a(1:3)=vab(1:3)*drij

   if (r .lt. 25) then
      x    =zab(iz1,iz2)
      alpha=r0ab(iz1,iz2)
      t27  =x*dexp(-alpha*r)
      e0   =t27*oner
      e    =e + e0*eps2(nk)
      drij=eps2(nk)*t27*(alpha*r+1.0d0)*oner/r2
      g_local_a(1:3)=g_local_a(1:3)-vab(1:3)*drij
   end if
   
   g(1:3,i1)=g(1:3,i1)+g_local_a(1:3)
   g(1:3,i2)=g(1:3,i2)-g_local_a(1:3)


!
!     the components of the virial tensor, if needed
!
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)
      vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)
      vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)
      vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)
      vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)
      vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)
      vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)
      vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)
      vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)
   end if


end do

!
!    Add a second component for the inter-molecule interaction for solvent QMDFFs
!
if (nmols .gt. 1) then
!
!    Loop over all atoms and calculate all Van der Waals and Coulomb terms to all atoms 
!    of other molecules
!
   do i=1,n-1   ! Loop to n or n-1?
      do j=i+1,n
         if (molnum(i) .ne. molnum(j)) then
            i1=i
            i2=j

            vab=xyz(1:3,i1)-xyz(1:3,i2)
!
!     apply periodic boundaries, if needed
!
            if (periodic) then
               call box_image(vab)
            end if

            r2=dot_product(vab,vab)

            r =sqrt(r2)
!
!     only proceed if the distance is below cutoff
!
            if (periodic .and. r .gt. vdw_cut) cycle

            oner =1.0d0/r

            iz1=at(i1)
            iz2=at(i2)
!
!     If noncovalent parameters were optimized separately, apply the factors
!
            if (ff_mod_noncov) then
               R0=r0094(iz1,iz2)*mn_par((i-1)*7+2)*mn_par((j-1)*7+2)
               c6=c66ab(i2,i1)*mn_par((i-1)*7+3)*mn_par((j-1)*7+3)
               r4=r2*r2
               r6=r4*r2
               r06=R0**6
               t6=r6+r06
               t8=r6*r2+r06*R0*R0
               c6t6=c6/t6
               c6t8=c6/t8
               t27=sr42(iz1,iz2)*c6t8*mn_par((i-1)*7+4)*mn_par((j-1)*7+4)
            else 
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
            end if
            e0=c6t6+t27
            e=e-e0
            drij=(c6t6*6.0d0*r4/t6+8.0d0*t27*r6/t8)

             
!
!     determine gradient vector entries
!
            g_local_a(1:3)=vab(1:3)*drij

!
!    Van der Waals energy: only, if distance is below the 25 bohr cutoff!
!
            if (r .lt. 25) then
               if (ff_mod_noncov) then
                  x    =zab(iz1,iz2)*mn_par((i-1)*7+5)*mn_par((j-1)*7+5)
                  alpha=r0ab(iz1,iz2)*mn_par((i-1)*7+6)*mn_par((j-1)*7+6)
               else 
                  x    =zab(iz1,iz2)
                  alpha=r0ab(iz1,iz2)
               end if
               t27  =x*dexp(-alpha*r)
               e0   =t27*oner
               e    =e + e0
               drij=t27*(alpha*r+1.0d0)*oner/r2
               g_local_a(1:3)=g_local_a(1:3)-vab(1:3)*drij
            end if
!            write(*,*) "eges",e

            g(1:3,i1)=g(1:3,i1)+g_local_a(1:3)
            g(1:3,i2)=g(1:3,i2)-g_local_a(1:3)


!
!     the components of the virial tensor, if needed
!
            if (calc_vir) then
               vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)
               vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)
               vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)
               vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)
               vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)
               vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)
               vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)
               vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)
               vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)
            end if
         end if
      end do
   end do
end if
!
!     The Coulomb energy: Do usual summation for nonperiodic system over 
!      QMDFF and (if existent) multi-molecule parts
!
ewald=.false.
if (.not. ewald) then
   do k=1,nnci
      i1=nci(1,k)
      i2=nci(2,k)
      nk=nci(3,k)

      vab=xyz(1:3,i1)-xyz(1:3,i2)
!
!     apply periodic boundaries, if needed
!
      if (periodic) then
         call box_image(vab)
      end if

      r2=dot_product(vab,vab)


      r =sqrt(r2)
!
!     For periodic systems: proceed only for distances below cutoff
!      for distances between cutoff and lower cutoff, perform a 
!      switching function approach
!
      switch=1.d0
      if (periodic) then
         if (r .gt. coul_cut) cycle
         if (.not. zahn) then
            if (r .gt. cut_low) then 
               switch=exp_switch(cut_low,coul_cut,r)           
            end if
         end if
      end if      
  

      oner =1.0d0/r
!
!     Use the Zahn smooth cutoff method, if ordered
!
      if (zahn) then
         e0=q(i1)*q(i2)*((erfc(zahn_a*r)*oner)-zahn_par*(r-coul_cut))
      else 
         e0=q(i1)*q(i2)*oner*eps1(nk)*switch
      end if
      e=e+e0
      drij=e0/r2

      g_local_a(1:3)=-vab(1:3)*drij

      g(1:3,i1)=g(1:3,i1)+g_local_a(1:3)
      g(1:3,i2)=g(1:3,i2)-g_local_a(1:3)

!
!     the components of the virial tensor, if needed
!
!
!     the components of the virial tensor, if needed
!
      if (calc_vir) then
         vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)
         vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)
         vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)
         vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)
         vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)
         vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)
         vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)
         vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)
         vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)
      end if
   end do
!
!    Add a second component for the inter-molecule interaction for solvent QMDFFs
!
!   e=0
   if (nmols .gt. 1) then

      do i=1,n-1
         do j=i+1,n
            if (molnum(i) .ne. molnum(j)) then
               i1=i
               i2=j


               vab=xyz(1:3,i1)-xyz(1:3,i2)
!
!     apply periodic boundaries, if needed
!
               if (periodic) then
                  call box_image(vab)
               end if

               r2=dot_product(vab,vab)

               r =sqrt(r2)
!
!     For periodic systems: proceed only for distances below cutoff
!
               switch=1.d0
               if (periodic) then
                  if (r .gt. coul_cut) cycle
                  if (.not. zahn) then
                     if (r .gt. cut_low) then
                        switch=exp_switch(cut_low,coul_cut,r)
                     end if
                  end if
               end if


               oner =1.0d0/r
       
!
!     Use the Zahn smooth cutoff method, if ordered
!
               if (zahn) then
                  e0=q(i1)*q(i2)*((erfc(zahn_a*r)*oner)-zahn_par*(r-coul_cut))
               else
                  e0=q(i1)*q(i2)*oner*switch
               end if

!
!     If local optimization of parameters was done 
! 
               if (ff_mod_noncov) then
                  e0=e0*mn_par(1)*mn_par(1)
               end if
               e=e+e0
               
               drij=e0/r2

               g_local_a(1:3)=-vab(1:3)*drij

               g(1:3,i1)=g(1:3,i1)+g_local_a(1:3)
               g(1:3,i2)=g(1:3,i2)-g_local_a(1:3)

!
!     the components of the virial tensor, if needed
!
!
!     the components of the virial tensor, if needed
!
               if (calc_vir) then
                  vir_ten(1,1)=vir_ten(1,1)+vab(1)*g_local_a(1)
                  vir_ten(2,1)=vir_ten(2,1)+vab(2)*g_local_a(1)
                  vir_ten(3,1)=vir_ten(3,1)+vab(3)*g_local_a(1)
                  vir_ten(1,2)=vir_ten(1,2)+vab(1)*g_local_a(2)
                  vir_ten(2,2)=vir_ten(2,2)+vab(2)*g_local_a(2)
                  vir_ten(3,2)=vir_ten(3,2)+vab(3)*g_local_a(2)
                  vir_ten(1,3)=vir_ten(1,3)+vab(1)*g_local_a(3)
                  vir_ten(2,3)=vir_ten(2,3)+vab(2)*g_local_a(3)
                  vir_ten(3,3)=vir_ten(3,3)+vab(3)*g_local_a(3)
               end if
            end if
         end do
      end do
   end if
!
!     The Coulomb energy: for a periodic system, apply the particle mesh
!      Ewald summation technique!
!     Implemented is the Smooth Particle Mesh Ewald Method according to 
!      J. Chem. Phys. 103, 8577 (1995)
!
else if (ewald) then
!
!     Recalculate the different adjustable parameter, mainly for scaling 
!     and accuracy
!
!   if (n .lt. 200) then
!      a_ewald=spi*(trtf*200/volbox**2)**(1d0/6d0) 
!   else 
!      a_ewald=spi*(trtf*n/volbox**2)**(1d0/6d0)  
!   end if

!   r_ew_cut=sqrt_p/a_ewald
!   k_ew_cut=2d0*sqrt_p*a_ewald
   coul_cut=r_ew_cut**2   
 
   if (r_ew_cut .gt. min(boxlen_x2,boxlen_y2,boxlen_z2)) then
      
      write(*,*) "The Ewald cutoff distance is larger than the half of the " 
      write(*,*) "periodic box! This would cause problems!"
      write(*,'(a,f14.7,a,f14.7)') " Half of periodic box: ",boxlen_x2*0.52917721092d0, &
              & " Ewald cutoff: ",r_ew_cut*0.52917721092d0
      write(*,*) "Increase the system size or decrease the Ewald accuracy!"
      call fatal
   end if

!
!     Reset the values of Ewald energy and gradient
! 
   e_ew=0.d0
   do i=1,n
      de_ew(1,i)=0.d0
      de_ew(2,i)=0.d0
      de_ew(3,i)=0.d0
   end do
 

!
!     Auxiliary parameter 
!   
   fs_par=-a_ewald/spi
   write(*,*) "fspar",fs_par
!
!     Calculate the Ewald self-energy term over all charges/atoms
!     --> constant term, no derivative!
!
   do i=1,n
      e_ew=e_ew - fs_par*q(i)**2
    end do

   write(*,*) "self enrergy",e_ew   
!
!     Compute the cell dipole boundary correction term
!     --> not needed for usual Ewald sums..
!
!   xd = 0.0d0
!   yd = 0.0d0
!   zd = 0.0d0
!   do i = 1, n
!      xd = xd + q(i)*xyz(1,i)
!      yd = yd + q(i)*xyz(2,i)
!      zd = zd + q(i)*xyz(3,i)
!   end do
!   term = (2.0d0/3.0d0) * (pi/volbox)
!   e_ew = e_ew + term * (xd*xd+yd*yd+zd*zd)


   do i = 1, n
      de = 2.0d0 * term * q(i)
      dedx = de * xd
      dedy = de * yd
      dedz = de * zd
      de_ew(1,i) = de_ew(1,i) + dedx
      de_ew(2,i) = de_ew(2,i) + dedy
      de_ew(3,i) = de_ew(3,i) + dedz
   end do

 !  write(*,*) "energy dop",e_ew
!
!     Compute the reciprocal space Ewald energy and gradient
!     In a separate subroutine for the full system of interest!
!
   call ewald_recip(n,xyz,q,e_rec,de_rec)
   write(*,*) "after rec",e_rec

!
!     Compute the real space Ewald energy and gradient for the 
!     local intramolecule interactions
!
!   e_ew=0
   do k=1,nnci
      i1=nci(1,k)
      i2=nci(2,k)
      nk=nci(3,k)
      dx=xyz(1,i1)-xyz(1,i2)
      dy=xyz(2,i1)-xyz(2,i2)
      dz=xyz(3,i1)-xyz(3,i2)
     
      call box_image(dx,dy,dz)

      r2=dx*dx+dy*dy+dz*dz
!
!     Proceed only if the squared distance is smaller than the cutoff
!
      if (r2 .gt. coul_cut) cycle

      r =sqrt(r2)
      oner =1.0d0/r
      r_buf=r

      r_buf2=r_buf*r_buf
      f12=q(i1)*q(i2)
      r_ewald=r*a_ewald

      erfterm=erfc(r_ewald)

      e_local=(f12/r_buf)*erfterm
      de=-f12*(erfterm/r_buf2+(2d0*a_ewald/spi)*exp(-r_ewald**2)/r)
      
      de=de/r
      dedx=de*dx
      dedy=de*dy
      dedz=de*dz

      e_ew=e_ew+e_local
      de_ew(1,i1)=de_ew(1,i1)+dedx
      de_ew(2,i1)=de_ew(2,i1)+dedy
      de_ew(3,i1)=de_ew(3,i1)+dedz
      de_ew(1,i2)=de_ew(1,i2)-dedx
      de_ew(2,i2)=de_ew(2,i2)-dedy
      de_ew(3,i2)=de_ew(3,i2)-dedz
   end do

!
!     Compute the real space Ewald energy for the intermolecule 
!     interactions, if a solvent system is present
!     see eq. (2.3) in paper
!
!    e_ew=0.d0
    if (nmols .gt. 1) then

      do i=1,n
         do j=i+1,n
            if (molnum(i) .ne. molnum(j)) then
               i1=i
               i2=j
               dx=xyz(1,i1)-xyz(1,i2)
               dy=xyz(2,i1)-xyz(2,i2)
               dz=xyz(3,i1)-xyz(3,i2)
               if (periodic) then
                  call box_image(dx,dy,dz)
               end if

               r2=dx*dx+dy*dy+dz*dz
!
!     Proceed only if the squared distance is smaller than the cutoff
!
               if (r2 .gt. coul_cut) cycle

               r =sqrt(r2)
               oner =1.0d0/r

               e_ew=e_ew+q(i1)*q(i2)*erfc(a_ewald*r)*oner

    !           r_buf=r+ebuffer

    !           r_buf2=r_buf*r_buf
    !           f12=q(i1)*q(i2)
!               write(*,*) "f12",f12
    !           r_ewald=r*a_ewald
!               write(*,*) "REW",r_ewald
    !           erfterm=erfc(r_ewald)

    !           e_local=(f12/r_buf)*(erfterm)
             !  write(*,*) "ewald",k,e_local
    !           de=-f12*((erfterm)/r_buf2+(2d0*a_ewald/spi)*exp(-r_ewald**2)/r)

      !         de=de/r
      !         dedx=de*dx
      !         dedy=de*dy
      !         dedz=de*dz

      !         e_ew=e_ew+e_local
!               write(*,*) "eew",e_ew/bohr
         !      stop "GupgpU"
      !         de_ew(1,i1)=de_ew(1,i1)+dedx
      !         de_ew(2,i1)=de_ew(2,i1)+dedy
      !         de_ew(3,i1)=de_ew(3,i1)+dedz
      !         de_ew(1,i2)=de_ew(1,i2)-dedx
      !         de_ew(2,i2)=de_ew(2,i2)-dedy
      !         de_ew(3,i2)=de_ew(3,i2)-dedz
         !   else 
         !      i1=i
         !      i2=j
         !      dx=xyz(1,i1)-xyz(1,i2)
         !      dy=xyz(2,i1)-xyz(2,i2)
         !      dz=xyz(3,i1)-xyz(3,i2)
         !      if (periodic) then
         !         call box_image(dx,dy,dz)
         !      end if

         !      r2=dx*dx+dy*dy+dz*dz
!
!     Proceed only if the squared distance is smaller than the cutoff
!
              ! if (r2 .gt. coul_cut) cycle

         !      r =sqrt(r2)
         !      oner =1.0d0/r
             !  e_ew=e_ew-0.5*q(i1)*q(i2)*erf(a_ewald*r)*oner

            end if
         end do
      end do
   end if
   write(*,*) "energy Reciproc",e_rec
   e_ew=e_ew-e_rec
   write(*,*) "energy Ewald",e_ew!/bohr

   if (ewald_brute) then
      write(*,'(a,f18.9,a)') " The SPME-energy of the system is ",e_ew," Hartrees."
   end if      
!   stop "in Ewald!"

end if
!
!     Calculate the charge center of mass for dipole moment calculation
!
if (sum(q(1:n)) .gt. 0.1d0) then
   charge_com=0.d0
   do i=1,n
      do j=1,3
         charge_com(j)=charge_com(j)+q(i)*xyz(j,i)
      end do
   end do
   charge_com=charge_com/sum(q(1:n))
else
   charge_com=0.d0
end if
 
!
!    Calculate the dipole moment relative to the center of charge
!
dipolemom=0.d0
!dipolemom=sum(charges*(xyz_local-charge_com),dim=2)
do i=1,n
   dipolemom=dipolemom-q(i)*(xyz(:,i)-charge_com(:))
end do
!
!     For intensities of numerical frequencies: store dipole vector in global array
!
if (calc_freq_int) then
   dip_list(:,int_incr) = dipolemom
end if



!
!    Benchmark for Smooth particle mesh Ewald method: calculate the 
!    brute force direct Ewald sum!
!
ewald_brute=.false.
if (ewald_brute) then
   e_brute=0.d0
   write(*,*)
   write(*,*) "The brute force Ewald method will be applied for the Coulomb"
   write(*,*) " energy. This might take some time..."


!
!    Part 1: the usual Coulomb energy with all particles in the box
!   
!    1A: The intramolecular energy contribution
!
   do k=1,nnci
      i1=nci(1,k)
      i2=nci(2,k)
      nk=nci(3,k)
      dx=xyz(1,i1)-xyz(1,i2)
      dy=xyz(2,i1)-xyz(2,i2)
      dz=xyz(3,i1)-xyz(3,i2)

      r2=dx*dx+dy*dy+dz*dz
      r =sqrt(r2)
      oner =1.0d0/r

      e0=q(i1)*q(i2)*oner*eps1(nk)

      e_brute=e_brute+e0
   end do
!
!   1B: The intermolecular energy contribution (if existent)
!
   if (nmols .gt. 1) then

      do i=1,n
         do j=i+1,n
            if (molnum(i) .ne. molnum(j)) then
               i1=i
               i2=j
               dx=xyz(1,i1)-xyz(1,i2)
               dy=xyz(2,i1)-xyz(2,i2)
               dz=xyz(3,i1)-xyz(3,i2)

               r2=dx*dx+dy*dy+dz*dz
               r =sqrt(r2)
               oner =1.0d0/r

               e0=q(i1)*q(i2)*oner

               e_brute=e_brute+e0
            end if
         end do
      end do
   end if
!
!    Part 2: The interaction to all cell replicas in the environment
!
   nbox=3
   nbox_sq=(nbox+1)**2
   allocate(pot_shell(nbox+1))
   pot_shell = 0.0d0
   do xbox = -nbox, nbox
      do ybox = -nbox, nbox
         do zbox = -nbox, nbox
            rbox_sq = xbox**2 + ybox**2 + zbox**2
            if ( rbox_sq > nbox_sq ) cycle ! Skip if outside maximum sphere of boxes
            rbox_vec = REAL ( [xbox,ybox,zbox] )*boxlen_x
            do i = 1, n
               do j = 1, n
                  if (rbox_sq .eq. 0 ) cycle ! Skip only for central box
                  ! Bare Coulomb term, including box vector, no periodic box correction
                  dx=xyz(1,i)-xyz(1,j)-rbox_vec(1)
                  dy=xyz(2,i)-xyz(2,j)-rbox_vec(2)
                  dz=xyz(3,i)-xyz(3,j)-rbox_vec(3)
             
                  r2=dx*dx+dy*dy+dz*dz
                  r =sqrt(r2)
                  oner =1.0d0/r

                  e0=q(i)*q(j)*oner
                  pot_shell(rbox_sq+1) = pot_shell(rbox_sq+1) + e0
               end do
            end do
         end do
      end do
   end do
   pot_shell=pot_shell/2.d0
  ! Convert to cumulative sum
   
   pot_shell(1) = e_brute
   do rbox_sq = 1, nbox_sq
      pot_shell(rbox_sq+1) = pot_shell(rbox_sq+1) + pot_shell(rbox_sq)
   end do

!
!     Write out the energies for the increasing shells to show convergence
!
   write(*,*) "The results for the full Ewald sum are:"
   write (*, fmt='(a)' ) ' Shell(box replicas) Energy (Hartree)'
   do rbox = 1, nbox
      write (*, fmt='(i5,f18.9)' ) rbox-1, pot_shell(rbox)
   end do

   write(*,*) "The calculation will be stopped here."
   stop
end if
!do i=1,n
!   write(*,*) "grad",g(:,i)
!end do

enb = enb + e
return
end subroutine ff_nonb
