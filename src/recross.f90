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
!     subroutine recross: Perform a complete recrossing factor 
!            calculation for the RPMD rate constant task
!
!     Disclaimer (parts (Recrossing factor calculation) were taken from):
!    
!     part of EVB
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
!
!   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
!                         William H. Green (whgreen@mit.edu)
!                         Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
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

subroutine recross(rank,psize,dt,xi_ideal,energy_ts,kappa)
use general
use evb_mod
implicit none
!
!     include MPI library
!
include 'mpif.h'
integer::i,j,k,l,m,o  ! loop index
integer::rank  ! the actual MPI processor ID
real(kind=8)::pi   ! pi..
real(kind=8)::dt  ! the timestep in atomic units
real(kind=8)::xi_ideal ! the xi value of the PMF max, set as constraint
real(kind=8)::xi_real ! the real xi value, here =xi_ideal
real(kind=8)::dxi_act(3,natoms)  ! the derivative of xi
real(kind=8)::d2xi_act(3,natoms,3,natoms) ! the hessian of xi
real(kind=8)::derivs(3,natoms,nbeads)   ! the gradient vector for all beads
real(kind=8)::q_start(3,natoms,nbeads)  ! the initial ts structure for reset 
real(kind=8)::q_save(3,natoms,nbeads)  ! the start structure for a bunch of childs
real(kind=8)::p_save(3,natoms,nbeads)  ! the start momentum for a child pair
real(kind=8)::q_1b(3,natoms)   ! coordinate vector for one bead
real(kind=8)::act_coord(3,natoms) ! coordinates for xi determination
real(kind=8),dimension(3,natoms)::derivs_1d  ! for gradient calculation
real(kind=8)::kappa_num(child_evol)  ! the time dependent recrossing factor
real(kind=8)::Vpot ! the potential energy of the used PES
real(kind=8)::en_act  ! the actual energy
real(kind=8)::fs,vs  ! recrossing factor fractions
real(kind=8),allocatable::k_save(:)  ! the stored force constant
real(kind=8)::kappa_denom  ! denominator for the recrossing factor calculation
real(kind=8)::num_total(child_evol),denom_total  ! the final kappa calculation
real(kind=8)::num_local(child_evol),denom_local ! the kappa calculation of the workers
real(kind=8)::kappa  ! the final result!
real(kind=8)::dt_info ! time in ps for information
real(kind=8)::t_total,t_actual ! print out of elapsed parent times
real(kind=8)::test123  ! for MPI debug
real(kind=8)::energy_act ! the actual energy
real(kind=8)::energy_ts  ! energy of the transition state
integer::andersen_save ! frequency for andersen thermostat (stored)
integer::recross_status  ! if the recrossing continues a previous calculation
integer::readstat  ! for restart file io
integer::steps
!     for MPI parallelization
integer::ierr,ID
integer::psize,source,status(MPI_STATUS_SIZE)
integer::count,tag_mpi,dest
real(kind=8)::message(child_evol+2)
integer::maxwork,numwork
integer::loop_large   ! number of loops over all slaves
integer::loop_rest    ! number of loops over remaining tasks
integer::schedule   ! info message for process
integer::err_act,err_count,struc_reset,traj_error ! error handling variables

pi=3.1415926535897932384d0    ! the pi 
dt_info=0.001*dt*2.41888428E-2 ! time in ps for information write out
t_total=child_times*child_interv*0.001*dt*2.41888428E-2  ! the total
                         ! time in ps that the parent trajectory will evolve
t_actual=0.d0  ! the already vanished parent trajectory time
!
!     Allocate array for saved force constants
!
allocate(k_save(size(k_force)))
!
!     Set the time dependent recrossing factor to 0!
! 
num_total=0.d0
denom_total=0.d0

!
!     Check if there are signs for an already started calculation
!
recross_status=0
open(unit=56,file="recross_status",status="old",iostat=readstat)
if (readstat .eq. 0) then
   read(56,*,iostat=readstat) recross_status
   if (readstat .ne. 0) then
      recross_status=0
   end if
end if
close(56)

if (recross_status .ne. 0) then
   if (rank .eq. 0) then
      write(*,*) "The recrossing calculation has been started already!"
      write(*,'(a,i6,a,i6,a)') "  ",recross_status," out of ",child_times," bunches of "
      write(*,*) "child trajectories were already calculated."
      if (recross_status .ge. child_times) then
         write(*,*) "All child trajectories were already calculated!"
      else        
         write(*,'(a,i6)') " We will start directly with bunch No.",recross_status+1
      end if
   end if
end if

!
!     If a previous calculation has been detected, read in the num_total and denom_total
!     values from the output files!
!
if (recross_status .ne. 0) then
   open(unit=78,file="recross_num_tmp.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      if (rank .eq. 0) then
         write(*,*) "The file 'recross_num_tmp.dat' is not there! We will start the recrossing"
         write(*,*) " calculation directly from beginning."
      end if
      recross_status=0
   else
      do i=1,child_evol
         read(78,*,iostat=readstat) num_total(i)
         if (readstat .ne. 0) then
            if (rank .eq. 0) then
               write(*,*) "The file 'recross_num_tmp.dat' is not complete! We will start the"
               write(*,*) " recrossing calculation directly from beginning."
            end if
            recross_status=0
            exit
         end if
      end do
   end if
   close(78)
   open(unit=78,file="recross_denom_tmp.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      if (rank .eq. 0) then
         write(*,*) "The file 'recross_denom_tmp.dat' is not there! We will start the recrossing"
         write(*,*) " calculation directly from beginning."
      end if
      recross_status=0
   else
      read(78,*,iostat=readstat) denom_total
      if (readstat .ne. 0) then
         if (recross_status .lt. child_times) then     
            if (rank .eq. 0) then
               write(*,*) "The file 'recross_denom_tmp.dat' is not complete! We will start the"
               write(*,*) " recrossing calculation directly from beginning."
            end if
            recross_status=0
         end if
      end if
   end if
   close(78)

   open(unit=78,file="recross_parent_pos.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      if (rank .eq. 0) then
         write(*,*) "The file 'recross_parent_pos.dat' is not there! We will start the recrossing"
         write(*,*) " calculation directly from beginning."
      end if
      recross_status=0
   else
      outer: do i=1,nbeads
         do j=1,natoms
            read(78,*,iostat=readstat) q_i(:,j,i)
            if (readstat .ne. 0) then
               if (rank .eq. 0) then
                  write(*,*) "The file 'recross_parent_pos.dat' is not complete! We will start"
                  write(*,*) " the recrossing calculation directly from beginning."
                  exit outer
               end if
               recross_status=0
            end if
         end do
      end do outer
   end if
   close(78)
end if
!
!     store the start structure 
!
q_start=q_i
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
err_count=0
!call mdinit(derivs,xi_ideal,dxi_act,2,rank)
!if (recross_status .eq. 0) then
if (rank .eq. 0) then
   call mdinit(derivs,xi_ideal,dxi_act,2,rank)
!
!     first reset point
!
   112 continue

   write(15,'(a,f11.4,a,f11.4,a)') " Evolving parent trajectory at Xi= ",xi_ideal, &
           & " for",real(recr_equi)*dt_info," ps."
   write(*,'(a,f11.4,a,f11.4,a)') " Evolving parent trajectory at Xi= ",xi_ideal, &
           & " for",real(recr_equi)*dt_info," ps."
   call mdinit(derivs,xi_ideal,dxi_act,2,rank)
!
!     do dynamics with the verlet subroutine, constrain is activated
!
   if (recross_status .ne. 0) then
      steps=100
   else 
      steps=recr_equi
   end if
   do i=1,steps
      call verlet(i,dt,derivs,en_act,0d0,0d0,xi_ideal,xi_real,dxi_act,j,1,.false.,rank)
!
!     check the trajectory, if failed, go to beginning of trajectory
!     --> check after some equilibration time first!
!
      if (recross_check .and. i .gt. 1000) then
         call rpmd_check(i,child_interv,4,0,0,en_act,xi_ideal,xi_ideal,energy_ts,&
                  & err_act,struc_reset,err_count,traj_error)
!
!     If too many errors occured, terminate the recrossing calculation, set a
!     default recrossing coefficient of 0.5 and go back to main program
!
         if (err_count .gt. err_max) then
            write(*,*) "ERROR! The recrossing calculation cannot be finished due to"
            write(*,*) " a large amount of errors!"
            write(*,*) "A default recrossing value will be used instead (kappa=0.5)"
            kappa=0.5d0
            return
         end if

         if (traj_error .eq. 1) then
!
!     take the TS structure as good starting condition
!
            q_i=q_start 
            goto 112
         end if
      end if
   end do
end if
!     define default tag value
tag_mpi=0
!     define default count value
count=1
!     define the number of worksteps 
maxwork=child_times
!
!     Now perform the recrossing factor calculation in parallel
!     The parent trajectory will be calculated by one process and each
!     batch of child trajectories will be distributed between the slaves
!     All is done via the MPI master worker setup. The recrossing 
!     factor is transferred between the processes 
!     Before starting, broadcast the positions and momenta of the systems 
!     particles
!
call mpi_bcast(q_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(p_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(xi_ideal,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

call mpi_barrier(mpi_comm_world,ierr)
!     Take the final structure of the first trajectory as reset point
!
q_start=q_i
!
!     Set the time dependent recrossing factor to 0!
! 
num_total=0.d0
denom_total=0.d0

!
!     The master process (rank=0) it sends the actual child sampling tasks
!     to the worker that will do the samplings for the child 
!     trajectories
!

loop_large=int((child_point/2)/(psize-1))
loop_rest=mod(child_point/2,psize-1)
if (recross_status .ge. child_times) goto 166

if (rank .eq. 0) then

  

   t_actual=recross_status*child_interv*0.001*dt*2.41888428E-2
 
   do k=recross_status+1,child_times
      

!
!     reset the recross parameters for this child run
!
      q_save=q_i   ! store the current position to recap it in each child
      k_save=k_force  ! store the force constant for the next parent time
      k_force=0.d0  ! set the force to zero for all childs!
      andersen_save=andersen_step  ! save interval for thermostat
      andersen_step=0  ! no thermostat for the childs!

      write(15,'(a,i6,a,f11.4,a)') " Spawning",child_point," child trajectories."
      write(15,'(a,f11.4,a,f11.4,a)') "Already",t_actual," ps of in total",t_total,&
            & " ps are elapsed."
      write(*,'(a,i6,a,f11.4,a)') " Spawning",child_point," child trajectories."
      write(*,'(a,f11.4,a,f11.4,a)') "Already",t_actual," ps of in total",t_total,&
            & " ps are elapsed."

     
!
!     at first, give each slave a task until all full rounds are finished
!
      do i=1,loop_large+1
         do j=1,psize-1
            schedule=i
            if (i .eq. loop_large+1) then
               if (k .eq. child_times .and. j .gt. loop_rest) then
                  schedule=-1
               else
                  cycle
               end if
            end if
            dest=j
!
!     Send number of trajectory and starting structure for children to slave process
!
            call mpi_send(schedule, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
            call mpi_send(q_i,natoms*3*nbeads, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
         end do
!
!     recieve results from all slaves for the i'th round
!
         do j=1,psize-1
            if (i .lt. loop_large+1) then
               message=0.d0
               call mpi_recv(message, child_evol+2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,tag_mpi, &
                  & MPI_COMM_WORLD,status,ierr)
               denom_total=denom_total+message(2)
               num_total=num_total+(message(3:child_evol+2))
            end if
         end do

      end do

!
!     Now, do the remaining jobs (modulo)
!
      do j=1,loop_rest
         dest=j
!
!     Send number of trajectory and starting structure for children to slave process
!  
         call mpi_send(i, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
         call mpi_send(q_i,natoms*3*nbeads, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
      end do
      do j=1,loop_rest
         call mpi_recv(message, child_evol+2, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,tag_mpi, &
               & MPI_COMM_WORLD,status,ierr)
         denom_total=denom_total+message(2)
         num_total=num_total+(message(3:child_evol+2))
      end do
!
!     Write the current numerator and denominator values as well as the parent structure
!      to updated files for restarts
!
      open(unit=67,file="recross_num_tmp.dat",status="replace")
      do j=1,child_evol
         write(67,*) num_total(j)
      end do
      close(67)
      open(unit=67,file="recross_denom_tmp.dat",status="replace")
      write(67,*) denom_total
      close(67)
      open(unit=67,file="recross_parent_pos.dat",status="replace")
      do j=1,nbeads
         do l=1,natoms
            read(67,*,iostat=readstat) q_i(:,l,j)
         end do
      end do
      close(67)
!
!     Update the recross status file for possible restarts
!
      open(unit=67,file="recross_status",status="replace")
      write(67,*) k 
      close(67)

      write(15,'(a,f12.8)') "The actual recrossing coefficient is: ",&
                   & num_total(child_evol)/denom_total
      write(*,'(a,f12.8)') "The actual recrossing coefficient is: ",&
                   & num_total(child_evol)/denom_total

      k_force=k_save
      andersen_step=andersen_save
      andersen_step=0
      q_i=q_save
!
!     second reset point
!

      if (andersen_step .eq. 0) then
        andersen_step=int(dsqrt(dble(child_interv)))
      end if

      113 continue
!
!     Also redo the initialization of the dynamics
!   
      call mdinit(derivs,xi_ideal,dxi_act,2,rank)

!
!     Sample the parent trajectory to generate a new starting configuration
!     for the children

      do j=1,child_interv
         call verlet(i,dt,derivs,en_act,0d0,0d0,xi_ideal,xi_real,dxi_act,j,1,.false.,rank)
!
!     check the trajectory, if failed, go to beginning of trajectory
!

         if (recross_check) then
            call rpmd_check(j,child_interv,4,0,0,en_act,xi_ideal,xi_ideal,energy_ts,&
                  & err_act,struc_reset,err_count,traj_error)
!
!     If too many errors occured, terminate the recrossing calculation, set a
!     default recrossing coefficient of 0.5 and go back to main program
!
            if (err_count .gt. err_max) then
               write(*,*) "ERROR! The recrossing calculation cannot be finished due to"
               write(*,*) " a large amount of errors!"
               write(*,*) "A default recrossing value will be used instead (kappa=0.5)"
               kappa=0.5d0
               return
            end if
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
   end do
!
!     Switch off the remaining processors
!
   do j=1,loop_rest
      dest=j
      call mpi_send(-1, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
      call mpi_send(q_i,natoms*3*nbeads, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
   end do
else 
   source=0
   message=rank
   do 
      call mpi_recv(numwork, count, MPI_DOUBLE_PRECISION, 0,tag_mpi,MPI_COMM_WORLD,status,ierr)
      call mpi_recv(q_i,natoms*3*nbeads, MPI_DOUBLE_PRECISION, 0,tag_mpi,MPI_COMM_WORLD,status,ierr)

!
!     If no more samplings are to do, exit with the current worker
!
      if (numwork .eq. -1) then
         exit
      end if
!
!     reset the recross parameters for this child run
!
      num_local=0.d0
      denom_local=0.d0
      t_actual=0.d0
      q_save=q_i   ! store the current position to recap it in each child
      k_save=k_force  ! store the force constant for the next parent time
      k_force=0.d0  ! set the force to zero for all childs!
      andersen_save=andersen_step  ! save interval for thermostat
      andersen_step=0  ! no thermostat for the childs!
!
!     start pairs of two childs with momenta of different signs!
!
!         write(*,*) "before andersen!",p_i
      call andersen
!         write(*,*) "after andersen!",p_i
      p_save=p_i
      do k=1,2
         kappa_num=0.d0
         kappa_denom=0.d0
         if (k .eq. 1) then
            p_i=p_save
         else if (k .eq. 2) then
            p_i=-p_save
         end if
         q_i=q_save
         call get_centroid(act_coord)  !extract centroid from q_i
         call calc_xi(act_coord,xi_ideal,xi_real,dxi_act,d2xi_act,2)
!
!     calculate gradient for the whole ring polymer
!
         do l=1,nbeads
            q_1b=q_i(:,:,l)
            call gradient (q_1b,vpot,derivs_1d,l,rank)
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
         vs = vs/nbeads
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
! 
         do l=1,child_evol
     !       write(*,*) "centroid3!",q_i
            call verlet(l,dt,derivs,en_act,0d0,0d0,xi_ideal,xi_real,dxi_act,m,2,.false.,rank)
            if (xi_real .gt.0) then
               kappa_num(l)=kappa_num(l)+vs/fs
            end if
         end do
!
!     Avoid negative transmission coefficients.
!     Has no visible effect on benchmark reactions
!   
!     Seems to lead to too large recrossing coefficients! Better deactivate it!
!
!         if (kappa_num(child_evol) .ge. 0.d0) then
            num_local=num_local+kappa_num
!         end if
         denom_local=denom_local+kappa_denom
     !       stop "HUohouhuo"
      end do
!
!     use the configuration stored before the sampling of the child 
!     trajectories as if they never appeared!
!

      k_force=k_save
      andersen_step=andersen_save
      q_i=q_save
      andersen_step=0
      message(1)=0.d0
      message(2)=denom_local
      message(3:child_evol+1)=num_local
      call mpi_send(message, child_evol+2, MPI_DOUBLE_PRECISION, source,tag_mpi,MPI_COMM_WORLD,ierr)
   end do

end if
166 continue
call mpi_barrier(mpi_comm_world,ierr)
!
!     Write out the success message with the final recrossing factor
!     Only for the master process
!
kappa=num_total(child_evol)/denom_total
if (rank .eq. 0) then
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
   write(56,*) "#  t(fs)           kappa(t)    "
   write(56,*) "#"
   do i=1,child_evol
      write(56,*) i*dt*2.41888428E-2,num_total(i)/denom_total
   end do
   close(56)
!
!    write out the final recrossing factor as restart for later calculations
!
   open(unit=33,file="recross_finished",status="unknown")
   write(33,*) kappa
   close(33)
end if
call mpi_barrier(mpi_comm_world,ierr)
return
end subroutine recross
