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
!     subroutine interpolate: For tri- and tetramolecular reactions:
!     first, take the TS structure as starting point, then determine 
!     the center of mass of the whole system, then, extrapolate the center
!     of masses of the single educts that they have a distance of R_inf to
!     each other. 
!
!     part of EVB
!
subroutine pre_sample(rank,n_all,gen_step,dt,coords,pre_equi,dt_info)
use general
use evb_mod
implicit none
integer::i,j,k   ! loop indices
integer::n_all  ! number of structures to be interpolated
integer::istep  ! the current verlet step number
integer::gen_step  ! number of MD steps per trajectory
integer::rank   ! the MPI process rank
real(kind=8)::dt   ! the timestep in atomic units
real(kind=8)::dt_info  ! time step in femtoseconds
real(kind=8)::en_act   ! the actual energy
real(kind=8)::xi_val,xi_real ! ideal and real xi values (only dummies here)
real(kind=8)::dxi_act(3,natoms)   ! derivative of the reaction coordinate (dummy)
real(kind=8)::coords(3,natoms),coord_save(3,natoms) ! the actual coordinates (TS)
real(kind=8)::derivs(3,natoms)   ! the derivative of the potential energy
real(kind=8)::pre_equi(n_all+1,3,natoms)
real(kind=8)::Red(3,sum_eds,sum_eds)  ! all possible educt-educt distances
real(kind=8)::r_eds(sum_eds,sum_eds) ! the actual distane between two educts
real(kind=8)::com(3,sum_eds)   ! array with center of mass coordinates
character(len=50)::subfname  ! name of the current xi value

!
!     store the actual TS coordinates and reset them at the end of 
!     this subroutine
!
coord_save=coords
!
!     First, calculate the rij values at the TS as a starting point
!
if (rank .eq. 0) then
   write(*,*)
   write(*,*) "-------------------PART   0-------------------"
   write(*,*) " Presampling the structures to geat educts with"
   write(*,*) " correct distances between their COMs."
   write(*,*) "----------------------------------------------"
   write(15,*)
   write(15,*) "-------------------PART   0-------------------"
   write(15,*) " Presampling the structures to geat educts with"
   write(15,*) " correct distances between their COMs."
   write(15,*) "----------------------------------------------"
end if

allocate(r_refs(sum_eds,sum_eds))

call calc_com(coords,com)
do i=1,sum_eds
   do j=i+1,sum_eds
      Red(1,i,j)=com(j,1)-com(i,1)
      Red(2,i,j)=com(j,2)-com(i,2)
      Red(3,i,j)=com(j,3)-com(i,3)
      r_eds(i,j)=sqrt(Red(1,i,j)*Red(1,i,j) + Red(2,i,j)*Red(2,i,j) &
         & + Red(3,i,j)*Red(3,i,j))
   end do
end do

q_i(:,:,1)=coords
!write(*,*) "start",n_all,r_eds(1,2),r_eds(1,3),r_eds(2,3),r_inf
!
!     perform the biased samplings for all windows between xi_min and xi_max!
!
do i=1,n_all+1
!
!     interpolate all r_eds values to the current point
!
   xi_val=1.d0-umbr_dist*(i-1)
!
!     initiate the dynamics
!
   call mdinit_bias(xi_val,dxi_act,derivs,1,1)

   do j=1,sum_eds
      do k=j+1,sum_eds
         r_refs(j,k)=xi_val*r_eds(j,k)+(1-xi_val)*r_inf 
      end do
   end do
!
!     write out info messages
!
   write(*,'(a,f8.4,a,f11.4,a)') " Begin structure interpolation at Xi= ",xi_val, &
            & " for",real(gen_step)*dt_info," ps."
   write(15,'(a,f8.4,a,f11.4,a)') " Begin structure interpolation at Xi= ",xi_val, &
            & " for",real(gen_step)*dt_info," ps."
   do istep=1,gen_step

      call verlet(istep,dt,derivs,en_act,0d0,0d0,xi_val,xi_real,dxi_act,i,3,.false.)
   end do
   pre_equi(i,:,:)=q_i(:,:,1)
end do
!
!    write the sampled structures to a file for later read in
!
open(unit=55,file="presampled_struc.xyz",status="unknown")
do i=1,n_all+1
   write(55,*) natoms
   write(55,*) "No.",i
   do j=1,natoms
      write(55,*) name(j),pre_equi(i,:,j)*bohr
   end do
end do
close(55)
!
!    print success message and asure that later starts of calculations
!    with the same settings do not restart it
!
call system("touch presampling_finished")
write(15,*) "----------------------------------------------"
write(15,*) "Presampling of COM distances finished!"
write(15,*)
write(*,*) "----------------------------------------------"
write(*,*) "Presampling of COM distances finished!"
write(*,*)
!
!    reset the actual structure to that of the TS
!
q_i(:,:,1)=coord_save
coords=coord_save

return
end subroutine pre_sample
