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
!     subroutine mdinit: initializes the velocities and accelerations
!     for a molecular dynamics trajectory, including restarts
!     literature reference:
!     K. A. Feenstra, B. Hess and H. J. C. Berendsen, "Improving
!     Efficiency of Large Time-Scale Molecular Dynamics Simulations
!     of Hydrogen-Rich Systems", Journal of Computational Chemistry,
!     8, 786-798 (1999)
!
!     part of EVB
!
subroutine mdinit(derivs,xi_ideal,dxi_act,bias_mode,rank)
use general
use qmdff
use evb_mod
use pbc_mod

implicit none
integer::i,j  ! loop index
integer::rank ! dummy MPI parameter
integer::backup  ! backup for andersen_step
real(kind=8),dimension(3,natoms,nbeads)::derivs
real(kind=8),dimension(3,natoms)::derivs_1d,q_1b
real(kind=8),dimension(nat6)::act_int,int_ideal
real(kind=8)::dxi_act(3,natoms)
real(kind=8)::xi_ideal,xi_real
real(kind=8)::epot  ! the actual potential energy
real(kind=8)::centroid(3,natoms)  !  the centroid with all the 
                                  !  com's of all beads
integer::bias_mode  ! usual sampling traj.: 1, constraint: 2
real(kind=8)::qterm,k_B,ekt  ! for NHC initialization
parameter(k_B=0.316679D-5) ! the boltzmann constant in au/kelvin
!
!     Calculate new box dimensions, if the NPT ensemble is used
!
if (periodic) then
   volbox=boxlen_x*boxlen_y*boxlen_z
end if

!
!     initalize the random number generator for the thermostat
!
!rank=0
!if (calc_modus .eq. 1) stop "Jgdi"
!call random_init_local(rank)
!
!     get the potential energy and atomic forces
!     for each bead at once: define its current structure and collect 
!     all results in a global derivs array thereafter
!
do i=1,nbeads
   q_1b=q_i(:,:,i)
   call gradient (q_1b,epot,derivs_1d,i,rank)
    
   derivs(:,:,i)=derivs_1d
end do
!if (calc_modus .eq. 1) stop "Jgdi2"
!
!     Calculate the centroid (center of masses) for the 
!     ring polymer and use its positions for application 
!     of the force constant
!
call get_centroid(centroid)
!
!     Add the bias potential for umbrella samplings!
!
if (bias_mode .eq. 1) then
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,1)
else if (bias_mode .eq. 2) then
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,0)
end if

!if (calc_modus .eq. 1) stop "Jgdi3"
!
!     Reset the momentum to a pseudo random distribution
!     --> apply the andersen thermostat
!
!     Beta-mode:  do the same for the Nose-Hoover thermostat!
!
if (.not. nve) then
   backup = andersen_step
   andersen_step = 1
   if (thermostat .eq. 0 .or. thermostat .eq. 1) then
      if (read_vel) then
         do i=1,natoms
            p_i(:,i,1)=-vel_start(i,:)*mass(i)
         end do
      else 
         call andersen
      end if
   end if
   andersen_step = backup
else 
   p_i=0
end if
!
!     Set the initial values for the Nose-Hoover-Chain thermostat
!
if (thermostat .eq. 2) then
   call andersen
!   p_i=0.d0
   nfree=3.d0*natoms
   ekt = k_B * kelvin
   qterm = ekt * nose_q * nose_q
   do j = 1, maxnose
      qnh(j) = qterm
      vnh(j)=0.d0
      gnh(j)=0.d0
   end do
   qnh(1) = dble(nfree) * qnh(1)
!
!    If the Nose-Hoover-Chain barostat it used, set its initial values as well
!
   qterm = ekt * nose_tau * nose_tau
   qbar = dble(nfree+1) * qterm

end if
!
!     If the Berendsen barostat is used
!
if (barostat .eq. 1) then
   taupres=2.d0
   compress=0.000046d0*prescon/1000.d0
!
!     If the Nose-Hoover chain barostat is used
!
else if (barostat .eq. 2) then
   qterm = ekt * nose_tau * nose_tau
   qbar = dble(nfree+1) * qterm
end if

if (fix_atoms) then
   do i=1,fix_num
      p_i(:,fix_list(i),:)=0.d0
   end do
end if


return
end subroutine mdinit
