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
!     subroutine verlet_bias: Changed version of verlet.f90 (dynamic), here
!        with incorporation of an additional bias potential for umbrella
!        samplings within rpmd.f90
!
!     part of EVB
!
subroutine verlet_bias (istep,dt,xi_ideal,xi_real,dxi_act,energy,derivs,round,constrain)
use general 
use evb_mod
use debug
!  ambigious reference of z-coord to qmdff z(94) array!!
implicit none
integer::i,j,k,istep   ! loop indices and actual step number
real(kind=8)::dt,dt_2
real(kind=8)::etot,epot
real(kind=8)::eksum
real(kind=8)::energy   ! the actual energy of the ring polymer
real(kind=8)::temp,pres
real(kind=8),dimension(natoms)::charge  !dummy array for charges
!   for umbrella sampling
real(kind=8),dimension(3,natoms)::derivs_1d,q_1b
real(kind=8),dimension(3,natoms,nbeads)::derivs
real(kind=8),dimension(nat6)::act_int,int_ideal
real(kind=8)::xi_ideal,xi_real
real(kind=8)::dxi_act(3,natoms)
! for the free ring polymer step:
! array with normal mode force interactions between RPMD beads:
real(kind=8)::poly(4,nbeads)  
real(kind=8)::beta_n    ! the inverse temperature per bead
real(kind=8)::twown   ! inverse inverse temperature
real(kind=8)::pi_n    ! half circle section for each bead
real(kind=8)::wk,wt,wm  ! circle section scaled parameters 
real(kind=8)::cos_wt,sin_wt   ! sinus and cosinus value of the circle section
real(kind=8)::pi   ! the pi
real(kind=8)::infinity,dbl_prec_var  ! test if one of the coordinates is infinity
real(kind=8)::p_new ! the actual new momentum component
real(kind=8)::centroid(3,natoms)  !the centroid with all the com's of all beads
integer::round  ! the actual window to sample (for calculation of avg/variance)
integer::constrain  ! if the trajectory shall be constrained to the actual xi value
                    ! (only for recrossing calculations)
                    !  0 : usual umbrella sampling
                    !  1 : constrain to dividing surface
                    !  2 : child trajectory: no restrain at all
                    !  3 : for pre sampling with tri/tetramolecular reactions
integer::const_good  ! 0:success, 1:failure (in constrain_q)
integer::num,modul

parameter(pi=3.1415926535897932384626433832795029d0)
infinity = HUGE(dbl_prec_var)   ! set the larges possible real value

!
!     update the momentum (half time step)
!
!write(*,*) "pi_STAT1",p_i
!write(*,*) "dt",dt
p_i=p_i-0.5d0*dt*derivs
!
!     set momentum to zero for fixed atoms 
!
if (fix_atoms) then
   do i=1,fix_num
      p_i(:,fix_list(i),:)=0.d0
   end do
end if

!
!     update the positions: for one bead, do the usual verlet procedure
!
if (nbeads .eq. 1) then
   do i=1,3
      do j=1,natoms
         q_i(i,j,1)=q_i(i,j,1)+p_i(i,j,1)*dt/mass(j)
      end do
   end do
else 
!
!     For the ring polymer: calculate the harmonic free ring polymer 
!     interactions between the beats: do it in normal mode space
!     tranform with Fast Fourier transformation to the normal mode 
!     space and back thereafter to the cartesian space
!
!     Transform to normal mode space
!     --> What is done there, exactly??
!
   do i = 1, 3
      do j = 1, Natoms
         call rfft(p_i(i,j,:), Nbeads)
         call rfft(q_i(i,j,:), Nbeads)
      end do
   end do
!   write(*,*) "beta",beta
!   beta=789.43666151836840
!   mass=1837.1525556457066
!   write(*,*) "pi",pi
!   write(*,*) "mass",mass(1:natoms)
   do j = 1, Natoms

      poly(1,1) = 1.0d0
      poly(2,1) = 0.0d0
      poly(3,1) = dt / mass(j)
      poly(4,1) = 1.0d0

      if (Nbeads .gt. 1) then
         beta_n = beta / nbeads
         twown = 2.0d0 / beta_n
         pi_n = pi / Nbeads
         do k = 1, nbeads / 2
            wk = twown * dsin(k * pi_n)
            wt = wk * dt
            wm = wk * mass(j)
            cos_wt = dcos(wt)
            sin_wt = dsin(wt)
            poly(1,k+1) = cos_wt
            poly(2,k+1) = -wm*sin_wt
            poly(3,k+1) = sin_wt/wm
            poly(4,k+1) = cos_wt
         end do
         do k = 1, (Nbeads - 1) / 2
            poly(1,Nbeads-k+1) = poly(1,k+1)
            poly(2,Nbeads-k+1) = poly(2,k+1)
            poly(3,Nbeads-k+1) = poly(3,k+1)
            poly(4,Nbeads-k+1) = poly(4,k+1)
         end do
      end if

      do k = 1, Nbeads
         do i = 1, 3
            p_new = p_i(i,j,k) * poly(1,k) + q_i(i,j,k) * poly(2,k)
            q_i(i,j,k) = p_i(i,j,k) * poly(3,k) + q_i(i,j,k) * poly(4,k)
            p_i(i,j,k) = p_new
         end do
      end do 
   end do
!   stop "hpihip"
!
!     Transform back to Cartesian space
!
   do i = 1, 3
      do j = 1, Natoms
         call irfft(p_i(i,j,:), Nbeads)
         call irfft(q_i(i,j,:), Nbeads)
      end do
   end do

end if
!
!     Calculate the centroid (center of masses) for the 
!     ring polymer and use its positions for application 
!     of the force constant
!
call get_centroid(centroid)

!
!     constrain the system to the xi value if desired
!
const_good=0.d0
if (constrain .eq. 1) then
   call constrain_q(centroid,xi_ideal,dxi_act,const_good,dt)
end if
!return
!
!     get the potential energy and atomic forces
!     for each bead at once: define its current structure and collect 
!     all results in a global derivs array thereafter
!
!write(*,*) "act_struc",q_i
!write(*,*) "derivs,before,",derivs
if (const_good .eq. 0) then
   energy=0.d0
else 
   write(*,*) "The SHAKE algorithm for constraint exeeded the maximum number of"
   write(*,*) "iterations! Restarting the trajectory! (energy penalty given)"
   energy=100000.d0
end if
xi_test=xi_real  ! TEST TEST
do i=1,nbeads 
   q_1b=q_i(:,:,i)              
   call gradient (q_1b,epot,derivs_1d,i)
   derivs(:,:,i)=derivs_1d
   energy=energy+epot
end do
!
!     For Mechanochemistry calculations: Add the additional forces here 
!     to the gradient of the reference method!
!

if (add_force) then
   do i=1,nbeads
      derivs(:,force1_at,i)=derivs(:,force1_at,i)+force1_v(:)*force1_k
      derivs(:,force2_at,i)=derivs(:,force2_at,i)+force2_v(:)*force2_k
   end do
end if


!
!     Add the bias potential for umbrella samplings!
!
if (constrain .eq. 0) then
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,0)
else if (constrain .eq. 1) then
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,1)
else if (constrain .eq. 2) then
!   write(*,*) "derivs_before",derivs,xi_real
!   write(*,*) "potential!",derivs
!   write(*,*) "centrooid!",centroid
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,1)
!   write(*,*) "xiiiiss",xi_real
!   write(*,*) "dxii",dxi_act
!   stop "Bogougvo"
!   write(*,*) "derivs_after",derivs,xi_real
else if (constrain .eq. 3) then
   call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,0)
end if
!
!     update the momentum (half time step)
!
p_i=p_i-0.5d0*dt*derivs
!write(*,*) "pi_STAT3",p_i
!write(*,*);write(*,*)

!
!     constrain the momentum (zero for dxi) if desired
!
if (constrain .eq. 1) then
   call constrain_p(dxi_act)
end if


!
!     Apply andersen thermostat (if chosen) to apply an random hit and 
!     change the momentum (only every andersen_step steps!)
!     --> Replace momenta with a fresh sampling from a Gaussian
!     distribution at the temperature of interest
!
if (constrain .ne. 2) then
   if (thermostat .eq. 0) then
      if (mod(istep,andersen_step) .eq. 0) then
         call andersen
      end if
   end if
end if
!
!     If the status is activated, update the average and the variance for this 
!     trajectors for the upfollowing umbrella-integration
!     ---> now shifted to main routine!
!
!if (constrain .eq. 0) then
!   if (writestat) then
!      average(round)=average(round)+xi_real
!      variance(round)=variance(round)+xi_real*xi_real
 !     if (abs(xi_ideal -1.0100) .lt. 0.0001) write(*,*) "average1,average2",average(round),variance(round)
 !     write(*,*) round,"average",average(round),"variance",variance(round)
!   end if
!end if

!
!     For debug mode: write all structures to one big file..
!     --> for all beads!
!
if (do_debug) then
   write(debug_unit,*) natoms*nbeads
   write(debug_unit,*)
   do i=1,nbeads
      do j=1,natoms
         write(debug_unit,*) name(j),q_i(:,j,i)*bohr
      end do
   end do
end if
!write(*,*) "q_i",q_i
!write(*,*) "p_i",p_i
!write(*,*) "derivs",derivs
return
end subroutine verlet_bias
