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
!     subroutine mdinit_bias: initializes the velocities and momenta
!      before each new umbrella sampling trajectory
!
!     part of EVB
!
subroutine mdinit_bias(xi_ideal,dxi_act,derivs,thermos,bias_mode)
use general
use evb_mod
implicit none
integer::i  ! loop index
integer::thermos ! if the thermostat shall be applied 
                 ! (0: no,  1: yes)
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
!
!     initalize the random number generator for the thermostat
!
!if (calc_modus .eq. 1) stop "Jgdi"
!
!     get the potential energy and atomic forces
!     for each bead at once: define its current structure and collect 
!     all results in a global derivs array thereafter
!
do i=1,nbeads
   q_1b=q_i(:,:,i)
   call gradient (q_1b,epot,derivs_1d,i)
   derivs(:,:,i)=derivs_1d
end do
!if (calc_modus .eq. 1) stop "Jgdi2"
!
!     Calculate the centroid (center of masses) for the 
!     ring polymer and use its positions for application 
!     of the force constant
!
call get_centroid(centroid)
!if (calc_modus .eq. 1) stop "Jgdi3"
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
!     --> apply the andersen thermostat or the GLE for partly resampling
!
backup = andersen_step
andersen_step = 1
if ((thermos .eq. 1) .and. (thermostat .eq. 0)) then 
   call andersen
else if (thermostat .eq. 1) then 
   stop "No GLE implemented!"
end if
andersen_step = backup
!
!     set momentum to zero for fixed atoms 
!
if (fix_atoms) then
   do i=1,fix_num
      p_i(:,fix_list(i),:)=0.d0
   end do
end if
return
end subroutine mdinit_bias
