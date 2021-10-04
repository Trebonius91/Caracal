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
real(kind=8)::oner,sigr6
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
   r_OH1=abs(sqrt(dot_product(OH1_vec,OH1_vec))-r_zero)
   r_OH2=abs(sqrt(dot_product(OH2_vec,OH2_vec))-r_zero)
   r_HH=abs(sqrt(dot_product(HH_vec,HH_vec))-r_0HH)
!
!     Calculate the OH bond stretches (Morse potentials)
!
   e_act=e_act+d_e*(1.d0-exp(-a_par*r_OH1))**2
   e_act=e_act+d_e*(1.d0-exp(-a_par*r_OH2))**2
!
!     Calculate the HH bond stretch (Harmonic potential)
!
   
   e_act=e_act+0.5d0*k_theta*r_HH**2
!
!     Calculate the bond-angle coupling term
! 
   e_act=e_act+k_rtheta*r_HH*(r_OH1+r_OH2)
!
!     Calculate the bond-bond coupling term
!
   e_act=e_act+k_rr*r_OH1*r_OH2
end do
!
!      Second, calculate the intermolecular energy part
!
!goto 234
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
!         if (rij .lt. coul_cut) then
!
!      The Coulomb interaction
! 
            oner=1.d0/rij
 !           write(*,*) q(i),q(j)
            e_act=e_act+q(i)*q(j)*oner
!
!      The dispersion (LJ 12-6) interaction: only between O and H!
!
            if (name(i) .eq. "O" .and. name(i) .eq. "O") then
               sigr6=(sigma_OO*oner)**6
               e_act=e_act+4d0*eps_OO*(sigr6**2-sigr6)
            end if

 !        end if         

      end if
   end do
end do
!stop "gup"
!234 continue
return
end subroutine egrad_water
