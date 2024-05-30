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
!     subroutine rpmd_check: check if an error occured during a
!       trajectory ran into the program rpmd.x
!       There are four kinds of trajectories: 
!       - structure generation (type 1)
!       - umbrella equilibration (type 2)
!       - umbrella sampling (type 3)
!       - recrossing parent (type 4)
!       For each of these types, perform the suited operations and 
!       reset appropriate results or the structures 
!
!       part of EVB
!
subroutine rpmd_check(istep,allstep,traj_type,umbr_run,num_umbrella,act_energy,xi_ideal,&
                    & xi_real,ts_energy,err_act,struc_reset,err_count,traj_error)
use general
use evb_mod 
use debug
implicit none
integer::i,j,k  ! loop indices
integer::traj_type  ! which type of trajectory shall be checked (see above)
integer::umbr_run   ! if type = 2 or 3, the number of umbrella window
integer::num_umbrella  ! total number of umbrella windows
integer::istep ! current step of verlet trajectory
integer::allstep ! total number of steps for current trajectory
real(kind=8)::act_energy  ! the actual energy of the system 
real(kind=8)::ts_energy   ! energy of the systems TS
real(kind=8)::infinity,dbl_prec_var  ! test if one of the coordinates is infinity
real(kind=8)::xi_real,xi_ideal   ! compare Xi values, if they are too different...
integer::traj_error ! result: 0 if OK, 1 if error
integer::err_count  ! total number of errors, will be incremented if a problem arises
integer::err_act    ! number of errors for this umbrella position
integer::struc_reset  ! number of newly chosen start structures 
integer::struc_index  ! if the index of the start structure will be changed


infinity = HUGE(dbl_prec_var)   ! set the largest possible real value
traj_error=0
!
!     PART A: Infinities and NaN values (complete fail of calculation)
!
!     Check structure: If an error occured in the trajectory and one of the 
!     entries is either NaN or Infinity
!
do i=1,3
   do j=1,natoms
      do k=1,nbeads
         if (q_i(i,j,k) .ne. q_i(i,j,k) .or. q_i(i,j,k) .gt. infinity) then
            traj_error=1
            write(*,*) "ERROR! The actual structure is not acceptable!"
            write(*,*) "The value of atom",j,", bead",k,", coordinate",j," is",q_i(i,j,k),"."
            goto 11
         end if
      end do
   end do
end do
11 continue
!
!     Check resulted energy, if its NaN or infinity
!
if (act_energy .ne. act_energy .or. act_energy .gt. infinity) then
   traj_error=1
   write(*,*) "ERROR! The energy of the actual structure is not acceptable! (E=",act_energy,")"
end if
!
!     PART B: High but finite energies (exploded systems)
!   
!     Add a tolerance value to the TS energy, if the actual energy is still higher,
!     throw the error
!
if (act_energy .gt. (ts_energy+energy_tol)*nbeads) then
   traj_error=1
   write(*,*) "ERROR! The energy of the actual structure (",act_energy," Eh) is much"
   write(*,*) "higher than the maximum allowed energy! (",(ts_energy+energy_tol)*nbeads,"Eh!"
   write(*,*) "Add the keyword RPMD_EN_TOL to change the tolerance energy per bead (kJ/mol)"
   write(*,*) "Actual value:",energy_tol*hartree*joule,"kJ/mol"
end if

!
!     PART C: The position of the structure with respect to the reactive coordinate Xi
!     os too much elongated!
!
if (abs(xi_real-xi_ideal) .gt. xi_tol) then
   traj_error=1 
   write(*,*) "ERROR! The actual xi value (",xi_real,") of the structure is too much away"
   write(*,*) "from the ideal value (",xi_ideal,")!"
end if

!
!     If an error occured, start appropriate countermeasures
!
if (traj_error.eq.1) then
!
!     For loose check option, go directly back to the dynamics routine without doing 
!     anything, except for start structure generation
!
   if (loose_check .and. traj_type .ne. 1) then
      write(*,*) "Loose_check activated! Error will not be counted!"
      return
   end if

   err_count=err_count+1
   write(*,*) "The error occured after",istep,"verlet steps. Total number of errors:",err_count
   write(*,*) "Allowed number of errors:",err_max
   !   TEST write structure that generated error to xyz file
   write(289,*) natoms
   write(289,*) 
   do i=1,natoms
      write(289,*) name(i), q_i(:,i,1)*bohr
   end do
   !   TESR  write energy, coupling components of all internal coordinates 
!   write(389,*) "test",nat6,"energy      coupling"
!   do i=1,nat6
!      write(389,*) i,error_vec_ens(i),error_vec_v12(i)
!   end do
!
!     If more than the allowed total number of errors occured, terminate program
!     If this occurs during the recrossing calculation, do not terminate the program
!     but set the recrossing to a standard value instead!
!
   if (err_count .gt. err_max) then
      if (traj_type .eq. 4) then
         write(*,*) "Too many errors occured during the recrossing!"
         return
      else 
         write(*,*) "Too many errors occured! The program will be stopped. Check your settings!"
         call fatal
      end if
   end if
!
!     Distinguish between different types of trajectories
!
   if (traj_type .eq. 1) then
      write(*,*) "Failed start structure generation trajectory. It will be restarted..."
      return
   else if (traj_type .eq. 2) then
      err_act=err_act+1
      write(*,*) "Failed umbrella equilibration trajectory..."
   else if (traj_type .eq. 3) then
      err_act=err_act+1
      write(*,*) "Failed umbrella sampling trajectory... Restart equilibration!"
   else if (traj_type .eq. 4) then
      write(*,*) "Failed recrossing parent trajectory... It will be restarted..."
      return
   end if
!
!     For umbrella equilibration or samplings: if more than act_err_max restarts failed,
!     choose new start structure for equilibration trajectory
!
   if (err_act .gt. err_act_max) then
      err_act=0
      write(*,*) "More than",err_act_max,"failures occured for this equilibration trajectory."
      write(*,*) "Another start structure will be taken from stored array."
      struc_reset=struc_reset+1
   end if
!
!     If in previous runs a new start structure index has been chosen:
!
   if (struc_reset .ge. 1) then
!
!     Up to four structure resets will be allowed! ideally, vary between increment and decrement
!     along the path, if boundaries are reached, continue on other sides
!
!     First go right  --x-- --> ---x-
      if (struc_reset .eq. 1) then   
         if (num_umbrella-umbr_run .lt. 1) then
            struc_index=struc_index-1
         else 
            struc_index=struc_index+1
         end if
!     Then go left   --x-- --> -x---

      else if (struc_reset .eq. 2) then
         if (umbr_run .lt. 1) then
            struc_index=struc_index+2
         else 
            if (num_umbrella-umbr_run .lt. 1) then
               struc_index=struc_index-2
            else 
               struc_index=struc_index-1
            end if
         end if
!     Then again go right --x-- --> ----x
      else if (struc_reset .eq. 3) then
         if (num_umbrella-umbr_run .lt. 2) then
            if (num_umbrella-umbr_run .lt. 1) then
               struc_index=struc_index-3
            else 
               struc_index=struc_index-2
            end if
         else 
            if (umbr_run .lt. 1) then 
               struc_index=struc_index+3
            else 
               struc_index=struc_index+2
            end if
         end if
!     Then again go left --x-- --> x----
      else if (struc_reset .eq. 4) then
          if (umbr_run .lt. 2) then
             if (umbr_run .lt. 1) then
                 struc_index=struc_index+4
             else 
                 struc_index=struc_index+3
             end if
          else 
             if (num_umbrella-umbr_run .lt. 2) then
                if (num_umbrella-umbr_run .lt. 1) then
                   struc_index=struc_index-4
                else 
                   struc_index=struc_index-3
                end if
             else 
                struc_index=struc_index-2
             end if
          end if 
      end if
!
!     Now store the structure which index is appointed
!      
      do j=1,nbeads
         q_i(1,:,j)=struc_equi(umbr_run+struc_index,1,:)
         q_i(2,:,j)=struc_equi(umbr_run+struc_index,2,:)
         q_i(3,:,j)=struc_equi(umbr_run+struc_index,3,:)
      end do
      write(*,*) "A new start structure will be taken: index",&
                  & umbr_run+struc_index,"of",num_umbrella,"."
   else 
!
!     Set back to actual start structure if number of errors wasnt reached
!
      do j=1,nbeads
         q_i(1,:,j)=struc_equi(umbr_run+struc_index,1,:)
         q_i(2,:,j)=struc_equi(umbr_run+struc_index,2,:)
         q_i(3,:,j)=struc_equi(umbr_run+struc_index,3,:)
      end do     
      write(*,*) "The actual start structure will be taken: index",umbr_run,"of",num_umbrella,"."
   end if
   return
end if
!
!     If also at the last step of the trajectory no error occured, the 
!     trajectory is fine; reset all local error parameters
!
if (istep .eq. allstep) then
   err_act=0 
   struc_reset=0
end if

return
end subroutine rpmd_check
