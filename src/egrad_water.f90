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
!     Subroutine egrad_water: Calculates energy and gradient of a box 
!      of water molecules with the flexible SPC water model
!
!     part of EVB
!

subroutine egrad_water(xyz_act,g_act,e_act)
use evb_mod
use general
use qmdff
implicit none
!     Loop indices
integer::i,j
!     The current coordinates (bohr)
real(kind=8),intent(in)::xyz_act(3,natoms)
!     The calculated gradient (hartree/bohr)
real(kind=8),intent(out)::g_act(3,natoms)
!     The calculated energy (hartree)
real(kind=8),intent(out)::e_act
integer::i_O,i_H1,i_H2
!     The flexible SPC parameters 
real(kind=8)::r_zero,theta_0,d_e,a_par,k_theta,k_rtheta
real(kind=8)::k_rr,e_H,e_O,sigma_OO,eps_OO,r_0HH
real(kind=8)::OH1_vec(3),OH2_vec(3),HH_vec(3)  ! the current distance vectors
real(kind=8)::dist_vec(3)
real(kind=8)::r_OH1,r_OH2,r_HH,rij  ! the current distances 
real(kind=8)::dr_OH1,dr_OH2,dr_HH ! the current distance deviations from minimum
real(kind=8)::oner,sigr6
real(kind=8)::exp1,exp2,de1,de2,e0
real(kind=8)::g_vec(3)
real(kind=8)::g_local_a(3),g_local_b(3)  ! for virial tensor
!
!     Initialize the local parameters 
!
r_zero=water_pars(1)
r_0HH=water_pars(2)
d_e=water_pars(3)
a_par=water_pars(4)
k_theta=water_pars(5)
k_rtheta=water_pars(6)
k_rr=water_pars(7)
e_H=water_pars(8)
e_O=water_pars(9)
sigma_OO=water_pars(10)
eps_OO=water_pars(11)

!
!     First, calculate the intramolecular energy part 
!
e_act=0.d0
g_act=0.d0
do i=1,nwater
   i_O=(i-1)*3+1
   i_H1=(i-1)*3+2
   i_H2=(i-1)*3+3
!
!     Calculate all needed distances between the atoms 
!
   OH1_vec=xyz_act(:,i_O)-xyz_act(:,i_H1)
   if (periodic) then
      call box_image(OH1_vec)
   end if
   OH2_vec=xyz_act(:,i_O)-xyz_act(:,i_H2)
   if (periodic) then
      call box_image(OH2_vec)
   end if
   HH_vec=xyz_act(:,i_H1)-xyz_act(:,i_H2)
   if (periodic) then
      call box_image(HH_vec)
   end if
   r_OH1=sqrt(dot_product(OH1_vec,OH1_vec))
   dr_OH1=r_OH1-r_zero
   r_OH2=sqrt(dot_product(OH2_vec,OH2_vec))
   dr_OH2=r_OH2-r_zero
   r_HH=sqrt(dot_product(HH_vec,HH_vec))
   dr_HH=r_HH-r_0HH
 
!
!     Calculate the OH bond stretches (Morse potentials)
!
   exp1=exp(a_par*dr_OH1)
   exp2=exp(a_par*dr_OH2)

   e_act=e_act+d_e*(1.d0-exp1)**2
   e_act=e_act+d_e*(1.d0-exp2)**2
!
!     Calculate gradients of OH bond stretches
!
   de1=-2d0*a_par*d_e*exp1*(1.d0-exp1)/r_OH1

   g_local_a(:)=de1*OH1_vec(:)
   g_act(:,i_H1)=g_act(:,i_H1)-g_local_a(:)
   g_act(:,i_O)=g_act(:,i_O)+g_local_a(:)

!
!     the components of the virial tensor, if needed
! 
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+OH1_vec(1)*g_local_a(1)
      vir_ten(2,1)=vir_ten(2,1)+OH1_vec(2)*g_local_a(1)
      vir_ten(3,1)=vir_ten(3,1)+OH1_vec(3)*g_local_a(1)
      vir_ten(1,2)=vir_ten(1,2)+OH1_vec(1)*g_local_a(2)
      vir_ten(2,2)=vir_ten(2,2)+OH1_vec(2)*g_local_a(2)
      vir_ten(3,2)=vir_ten(3,2)+OH1_vec(3)*g_local_a(2)
      vir_ten(1,3)=vir_ten(1,3)+OH1_vec(1)*g_local_a(3)
      vir_ten(2,3)=vir_ten(2,3)+OH1_vec(2)*g_local_a(3)
      vir_ten(3,3)=vir_ten(3,3)+OH1_vec(3)*g_local_a(3)
   end if

   de2=-2d0*a_par*d_e*exp2*(1.d0-exp2)/r_OH2
  
   g_local_a(:)=de2*OH2_vec(:) 
   g_act(:,i_H2)=g_act(:,i_H2)-g_local_a(:)
   g_act(:,i_O)=g_act(:,i_O)+g_local_a(:)

!
!     the components of the virial tensor, if needed
! 
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+OH2_vec(1)*g_local_a(1)
      vir_ten(2,1)=vir_ten(2,1)+OH2_vec(2)*g_local_a(1)
      vir_ten(3,1)=vir_ten(3,1)+OH2_vec(3)*g_local_a(1)
      vir_ten(1,2)=vir_ten(1,2)+OH2_vec(1)*g_local_a(2)
      vir_ten(2,2)=vir_ten(2,2)+OH2_vec(2)*g_local_a(2)
      vir_ten(3,2)=vir_ten(3,2)+OH2_vec(3)*g_local_a(2)
      vir_ten(1,3)=vir_ten(1,3)+OH2_vec(1)*g_local_a(3)
      vir_ten(2,3)=vir_ten(2,3)+OH2_vec(2)*g_local_a(3)
      vir_ten(3,3)=vir_ten(3,3)+OH2_vec(3)*g_local_a(3)
   end if

!
!     Calculate the HH bond stretch (Harmonic potential)
! 
   
   e_act=e_act+0.5d0*k_theta*dr_HH**2

!
!     Gradients of the HH bond stretch
!
   g_local_a(:)=k_theta*(HH_vec(:))*dr_HH/r_HH

   g_act(:,i_H1)=g_act(:,i_H1)+g_local_a(:)
   g_act(:,i_H2)=g_act(:,i_H2)-g_local_a(:)

!
!     the components of the virial tensor, if needed
! 

   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+HH_vec(1)*g_local_a(1)
      vir_ten(2,1)=vir_ten(2,1)+HH_vec(2)*g_local_a(1)
      vir_ten(3,1)=vir_ten(3,1)+HH_vec(3)*g_local_a(1)
      vir_ten(1,2)=vir_ten(1,2)+HH_vec(1)*g_local_a(2)
      vir_ten(2,2)=vir_ten(2,2)+HH_vec(2)*g_local_a(2)
      vir_ten(3,2)=vir_ten(3,2)+HH_vec(3)*g_local_a(2)
      vir_ten(1,3)=vir_ten(1,3)+HH_vec(1)*g_local_a(3)
      vir_ten(2,3)=vir_ten(2,3)+HH_vec(2)*g_local_a(3)
      vir_ten(3,3)=vir_ten(3,3)+HH_vec(3)*g_local_a(3)
   end if

!
!     Calculate the bond-angle coupling term
! 
   e_act=e_act+k_rtheta*dr_HH*(dr_OH1+dr_OH2)
!
!     Gradients of the bond-angle coupling term
!
   g_local_a(:)=dr_HH/r_OH1*OH1_vec(:)*k_rtheta
   g_local_b(:)=dr_HH/r_OH2*OH2_vec(:)*k_rtheta

   g_act(:,i_H1)=g_act(:,i_H1)+k_rtheta*2d0*dr_OH1/r_HH*HH_vec(:)- &
                    & g_local_a
   g_act(:,i_H2)=g_act(:,i_H2)-k_rtheta*2d0*dr_OH2/r_HH*HH_vec(:)- &
                    & g_local_b
   g_act(:,i_O)=g_act(:,i_O)+g_local_a + g_local_b

!
!     the components of the virial tensor, if needed
! 
   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+OH1_vec(1)*g_local_a(1)+OH2_vec(1)*g_local_b(1)
      vir_ten(2,1)=vir_ten(2,1)+OH1_vec(2)*g_local_a(1)+OH2_vec(2)*g_local_b(1)
      vir_ten(3,1)=vir_ten(3,1)+OH1_vec(3)*g_local_a(1)+OH2_vec(3)*g_local_b(1)
      vir_ten(1,2)=vir_ten(1,2)+OH1_vec(1)*g_local_a(2)+OH2_vec(1)*g_local_b(2)
      vir_ten(2,2)=vir_ten(2,2)+OH1_vec(2)*g_local_a(2)+OH2_vec(2)*g_local_b(2)
      vir_ten(3,2)=vir_ten(3,2)+OH1_vec(3)*g_local_a(2)+OH2_vec(3)*g_local_b(2)
      vir_ten(1,3)=vir_ten(1,3)+OH1_vec(1)*g_local_a(3)+OH2_vec(1)*g_local_b(3)
      vir_ten(2,3)=vir_ten(2,3)+OH1_vec(2)*g_local_a(3)+OH2_vec(2)*g_local_b(3)
      vir_ten(3,3)=vir_ten(3,3)+OH1_vec(3)*g_local_a(3)+OH2_vec(3)*g_local_b(3)
   end if

!
!     Calculate the bond-bond coupling term
!
   e_act=e_act+k_rr*dr_OH1*dr_OH2
!
!     Gradients of the bond-bond coupling term
!
   g_local_a(:)=k_rr*(dr_OH2/r_OH1*OH1_vec(:))
   g_local_b(:)=k_rr*(dr_OH1/r_OH2*OH2_vec(:))


   g_act(:,i_H1)=g_act(:,i_H1)-g_local_a(:)
   g_act(:,i_H2)=g_act(:,i_H2)-g_local_b(:)
   g_act(:,i_O)=g_act(:,i_O)+g_local_b(:)+ g_local_a(:)  

   if (calc_vir) then
      vir_ten(1,1)=vir_ten(1,1)+OH1_vec(1)*g_local_a(1)+OH2_vec(1)*g_local_b(1)
      vir_ten(2,1)=vir_ten(2,1)+OH1_vec(2)*g_local_a(1)+OH2_vec(2)*g_local_b(1)
      vir_ten(3,1)=vir_ten(3,1)+OH1_vec(3)*g_local_a(1)+OH2_vec(3)*g_local_b(1)
      vir_ten(1,2)=vir_ten(1,2)+OH1_vec(1)*g_local_a(2)+OH2_vec(1)*g_local_b(2)
      vir_ten(2,2)=vir_ten(2,2)+OH1_vec(2)*g_local_a(2)+OH2_vec(2)*g_local_b(2)
      vir_ten(3,2)=vir_ten(3,2)+OH1_vec(3)*g_local_a(2)+OH2_vec(3)*g_local_b(2)
      vir_ten(1,3)=vir_ten(1,3)+OH1_vec(1)*g_local_a(3)+OH2_vec(1)*g_local_b(3)
      vir_ten(2,3)=vir_ten(2,3)+OH1_vec(2)*g_local_a(3)+OH2_vec(2)*g_local_b(3)
      vir_ten(3,3)=vir_ten(3,3)+OH1_vec(3)*g_local_a(3)+OH2_vec(3)*g_local_b(3)
   end if


end do
!
!      Second, calculate the intermolecular energy part
!
do i=1,natoms
   do j=i+1,natoms
!
!      Only between atoms of different water molecules 
!      
      if (water_act(i) .ne. water_act(j)) then
         dist_vec=xyz_act(:,i)-xyz_act(:,j)
         if (periodic) then
            call box_image(dist_vec)
         end if
         rij=sqrt(dot_product(dist_vec,dist_vec))         
!
!      Only below the cutoff
!         
         if (rij .lt. coul_cut .or. .not. periodic) then
!
!      The Coulomb interaction
! 
            oner=1.d0/rij
!
!      Use the Zahn smooth cutoff method, if ordered
!
            if (zahn) then
               e0=q(i)*q(j)*((erfc(zahn_a*rij)*oner)-zahn_par*(rij-coul_cut))
            else 
               e0=q(i)*q(j)*oner
            end if

            e_act=e_act+e0
!
!      The Coulomb interaction gradient
!
            g_vec=e0*oner*oner*dist_vec(:)
            g_act(:,i)=g_act(:,i)-g_vec
            g_act(:,j)=g_act(:,j)+g_vec
!
!     the components of the virial tensor, if needed
! 
            if (calc_vir) then
               vir_ten(1,1)=vir_ten(1,1)-dist_vec(1)*g_vec(1)
               vir_ten(2,1)=vir_ten(2,1)-dist_vec(2)*g_vec(1)
               vir_ten(3,1)=vir_ten(3,1)-dist_vec(3)*g_vec(1)
               vir_ten(1,2)=vir_ten(1,2)-dist_vec(1)*g_vec(2)
               vir_ten(2,2)=vir_ten(2,2)-dist_vec(2)*g_vec(2)
               vir_ten(3,2)=vir_ten(3,2)-dist_vec(3)*g_vec(2)
               vir_ten(1,3)=vir_ten(1,3)-dist_vec(1)*g_vec(3)
               vir_ten(2,3)=vir_ten(2,3)-dist_vec(2)*g_vec(3)
               vir_ten(3,3)=vir_ten(3,3)-dist_vec(3)*g_vec(3)
             end if

!
!      The dispersion (LJ 12-6) interaction: only between O and H!
!
            if (name(i) .eq. "O" .and. name(i) .eq. "O") then
               sigr6=(sigma_OO*oner)**6
               e_act=e_act+4d0*eps_OO*(sigr6**2-sigr6)
!
!      The dispersion gradient
!
               g_vec=24.d0*eps_OO*oner*oner*dist_vec(:)* &
                             & sigr6*(2.d0*sigr6-1.d0)
               g_act(:,i)=g_act(:,i)-g_vec
               g_act(:,j)=g_act(:,j)+g_vec
!
!     the components of the virial tensor, if needed
! 
               if (calc_vir) then
                  vir_ten(1,1)=vir_ten(1,1)+dist_vec(1)*g_vec(1)
                  vir_ten(2,1)=vir_ten(2,1)+dist_vec(2)*g_vec(1)
                  vir_ten(3,1)=vir_ten(3,1)+dist_vec(3)*g_vec(1)
                  vir_ten(1,2)=vir_ten(1,2)+dist_vec(1)*g_vec(2)
                  vir_ten(2,2)=vir_ten(2,2)+dist_vec(2)*g_vec(2)
                  vir_ten(3,2)=vir_ten(3,2)+dist_vec(3)*g_vec(2)
                  vir_ten(1,3)=vir_ten(1,3)+dist_vec(1)*g_vec(3)
                  vir_ten(2,3)=vir_ten(2,3)+dist_vec(2)*g_vec(3)
                  vir_ten(3,3)=vir_ten(3,3)+dist_vec(3)*g_vec(3)
               end if
            end if
         end if         
      end if
   end do
end do
!stop "gup"
return
end subroutine egrad_water
