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
!     subroutine recross_serial: Perform a complete recrossing factor 
!            calculation for the RPMD rate constant task
!            NO MPI will be used here! (for machines with MPI bugs)    
!
!     part of EVB
!
subroutine recross_serial(rank,dt,xi_ideal,energy_ts,kappa)
use general
use evb_mod
implicit none
integer::i,j,k,l,m,o
integer::rank  ! the actual MPI processor ID
real(kind=8)::pi   ! pi..
real(kind=8)::dt  ! the timestep in atomic units
real(kind=8)::xi_ideal ! the xi value of the PMF max, set as constraint
real(kind=8)::xi_real ! the real xi value, here =xi_ideal
real(kind=8)::dxi_act(3,natoms)  ! the derivative of xi
real(kind=8)::d2xi_act(3,natoms,3,natoms) ! the hessian of xi
real(kind=8)::derivs(3,natoms,nbeads)
real(kind=8)::q_start(3,natoms,nbeads)  ! the initial ts structure for reset 
real(kind=8)::q_save(3,natoms,nbeads)  ! the start structure for a bunch of childs
real(kind=8)::p_save(3,natoms,nbeads)  ! the start momentum for a child pair
real(kind=8)::act_coord(3,natoms) ! coordinates for xi determination
real(kind=8)::q_1b(3,natoms)   ! coordinate vector for one bead
real(kind=8),dimension(3,natoms)::derivs_1d  ! for gradient calculation
real(kind=8)::kappa_num(child_evol)  ! the time dependent recrossing factor
real(kind=8)::Vpot ! the potential energy of the used PES
real(kind=8)::fs,vs  ! recrossing factor fractions
real(kind=8),allocatable::k_save(:)  ! the stored force constant
real(kind=8)::kappa_denom  ! denominator for the recrossing factor calculation
real(kind=8)::num_total(child_evol),denom_total  ! the final kappa calculation
real(kind=8)::kappa  ! the final result!
real(kind=8)::energy_act ! the actual energy
real(kind=8)::energy_ts  ! energy of the transition state
real(kind=8)::dt_info ! time in ps for information
real(kind=8)::t_total,t_actual ! print out of elapsed parent times
integer::andersen_save ! frequency for andersen thermostat (stored)
integer::err_act,err_count,struc_reset,traj_error ! error handling variables
!
!     number of errors for trajectory checking
!
err_count=0
err_act=0
pi=3.1415926535897932384d0    ! the pi 
dt_info=0.001*dt*2.41888428E-2 ! time in ps for information write out
t_total=child_times*child_interv*0.001*dt*2.41888428E-2  ! the total
                         ! time in ps that the parent trajectory will evolve
t_actual=0.d0  ! the already vanished parent trajectory time
!
!     store the start structure 
!
q_start=q_i
!
!     Allocate array for saved force constants
!
allocate(k_save(size(k_force)))

!
!     first, equilibrate the parent trajectory: umbrella trajectory with
!     infinite strong potential
!     start with the TS structure (should be a good initial guess)
!
!     Initialize the dynamics
!
andersen_step=0
if (andersen_step .eq. 0) then
   andersen_step=int(dsqrt(dble(recr_equi)))
end if
!
!     first reset point
!
112 continue
write(15,'(a,f11.4,a,f11.4,a)') " Evolving parent trajectory at Xi= ",xi_ideal, &
           & " for",real(recr_equi)*dt_info," ps."
write(*,'(a,f11.4,a,f11.4,a)') " Evolving parent trajectory at Xi= ",xi_ideal, &
           & " for",real(recr_equi)*dt_info," ps."
call mdinit(derivs,xi_ideal,dxi_act,2,1)
!
!     do dynamics with the verlet subroutine, constrain is activated
!
do i=1,recr_equi
   call verlet(i,dt,derivs,energy_act,0d0,0d0,xi_ideal,xi_real,dxi_act,j,1,.false.,1)
!
!     check the trajectory, if failed, go to beginning of trajectory
!     --> check after some equilibration time first!
!
   if (recross_check .and. i .gt. 1000) then
      call rpmd_check(i,child_interv,1,0,0,energy_act,xi_ideal,xi_ideal,energy_ts,&
               & err_act,struc_reset,err_count,traj_error)
      if (traj_error .eq. 1) then
!
!     take the TS structure as good starting condition
!
         q_i=q_start 
         goto 112
      end if
   end if
end do
!
!     Take the final structure of the first trajectory as reset point
!
q_start=q_i
!
!     Set the time dependent recrossing factor to 0!
! 
num_total=0.d0
denom_total=0.d0
!
!     After the parent trajectory is equilibrated, sample it child_interv
!     MD steps and start child_point child trajectories after it 
!     this will be repeated until all child trajectories are started
!     (this are child_times big cycles
!  
do i=1,child_times 
   q_save=q_i   ! store the current position to recap it in each child
   k_save=k_force  ! store the force constant for the next parent time
   k_force=0.d0  ! set the force to zero for all childs!
   andersen_save=andersen_step  ! save interval for thermostat
   andersen_step=0  ! no thermostat for the childs!
!
!     print info message before each spawn wafe of child trajectories
!  
   write(15,'(a,i6,a,f11.4,a)') " Spawning",child_point," child trajectories."
   write(15,'(a,f11.4,a,f11.4,a)') "Already",t_actual," ps of in total",t_total,&
            & " ps are elapsed." 
   write(*,'(a,i6,a,f11.4,a)') " Spawning",child_point," child trajectories."
   write(*,'(a,f11.4,a,f11.4,a)') "Already",t_actual," ps of in total",t_total,&
            & " ps are elapsed."
   
   do j=1,child_point/2
!
!     start pairs of two childs with momenta of different signs!
!
      call andersen
      p_save=p_i
      do k=1,2
         kappa_num=0.d0
         kappa_denom=0.d0
         if (k .eq. 1) then
            p_i=p_save
         else if (k .eq. 2) then
            p_i=-p_save
         end if
!         act_coord=q_save
         q_i=q_save
         call get_centroid(act_coord)
         call calc_xi(act_coord,xi_ideal,xi_real,dxi_act,d2xi_act,2)
!         write(*,*) "xi_act",xi_ideal,xi_real
!
!     calculate gradient for the whole ring polymer
!
         do l=1,nbeads
            q_1b=q_i(:,:,l)
            call gradient (q_1b,vpot,derivs_1d,l)
            derivs(:,:,l)=derivs_1d
         end do

!
!     calculate the recrossing velocity for the recrossing factor
!        
         vs = 0.0d0
         do l = 1, 3
            do m = 1, Natoms
               do o=1,nbeads
                  vs = vs + dxi_act(l,m) * p_i(l,m,o) / mass(m)
               end do
            end do
         end do
         vs = vs/nbeads  !later, correct by beads
!
!     calculate the recrossing flux for the recrossing factor
!
         fs = 0.0d0
         do l = 1, 3
            do m = 1, Natoms
                fs = fs + dxi_act(l,m) * dxi_act(l,m) / mass(m)
            end do
         end do
         fs = sqrt(fs / (2.0d0 * pi * beta))
!
!     increase the denominator of kappa if the velocity is positive
! 
         if (vs .gt. 0) then
            kappa_denom=kappa_denom+vs/fs
         end if
!
!     now do the actual evaluation of the child trajectory
!     no bias potential and no temperature shall be applied!
!
         do l=1,child_evol
            call verlet(l,dt,derivs,energy_act,0d0,0d0,xi_ideal,xi_real,dxi_act,m,2,.false.,1)
            if (xi_real .gt.0) then
               kappa_num(l)=kappa_num(l)+vs/fs
            end if
         end do
!
!     Avoid negative transmission coefficients.
!     Has no visible effect on benchmark reactions
!
!         if (kappa_num(child_evol) .ge. 0.d0) then
            num_total=num_total+kappa_num
!         end if
         denom_total=denom_total+kappa_denom
!         stop "gofozfoz"
      end do
   end do
!   stop "Gofzf"
!
!     use the configuration stored before the sampling of the child 
!     trajectories as if they never appeared!
!

   k_force=k_save
   andersen_step=andersen_save
   andersen_step=0

   if (andersen_step .eq. 0) then
     andersen_step=int(dsqrt(dble(child_interv)))
   end if

!   act_coord=q_save
   q_i=q_save
!
!     second reset point
!
   113 continue
!
!     Also redo the initialization of the dynamics
! 
   call mdinit(derivs,xi_ideal,dxi_act,2,1)

!
!     Sample the parent trajectory to generate a new starting configuration
!     for the childs
!
   do j=1,child_interv
      call verlet(i,dt,derivs,energy_act,0d0,0d0,xi_ideal,xi_real,dxi_act,j,1,.false.,1)
!
!     check the trajectory, if failed, go to beginning of trajectory
!

      if (recross_check) then
         call rpmd_check(j,child_interv,1,0,0,energy_act,xi_ideal,xi_ideal,energy_ts,&
                  & err_act,struc_reset,err_count,traj_error)
         if (traj_error .eq. 1) then
!
!     take the TS structure as good starting condition
!
            q_i=q_save
            goto 113
         end if
      end if

   end do
!
!     increment the informational time that has elapsed so far
!
   t_actual=t_actual+child_interv*0.001*dt*2.41888428E-2
   write(15,'(a,f12.8)') "The actual recrossing coefficient is: ",num_total(child_evol)/denom_total  
   write(*,'(a,f12.8)') "The actual recrossing coefficient is: ",num_total(child_evol)/denom_total
end do
!
!     Write out the success message with the final recrossing factor
!
kappa=num_total(child_evol)/denom_total
write(15,*)
write(15,*) "The recrossing calculation was successful!"
write(15,'(a,f12.8)') "The final recrossing factor is:",kappa
write(15,*)
write(*,*)
write(*,*) "The recrossing calculation was successful!"
write(*,'(a,f12.8)') "The final recrossing factor is:",kappa
write(*,*)
!
!     write out the time dependent recrossing factor!
!
open(unit=56,file="recrossing_time.dat",status="unknown")
write(56,*) "# This is the time dependent recrossing factor calculated with EVB-QMDFF!"
write(56,*) "#" 
write(56,*) "#  t(fs)           kappa(t)      kappa_num   kappa_denom"
write(56,*) "#"
do i=1,child_evol
   write(56,*) i*dt*2.41888428E-2,num_total(i)/denom_total,num_total(i),denom_total
end do
close(56)
!
!    write out the final recrossing factor as restart for later calculations
!
open(unit=33,file="recross_finished",status="unknown")
write(33,*) kappa
close(33)

return
end subroutine recross_serial
