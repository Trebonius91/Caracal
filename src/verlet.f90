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
!     subroutine verlet: Main subroutine for performing a time step in an 
!      arbitrary MD trajectory calculation (unbiased or biased)
!
!     part of EVB
!
!     Disclaimer (parts (RPMD integration) were taken from):
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

subroutine verlet (istep,dt,derivs,epot,ekin,afm_force,xi_ideal,xi_real,&
          & dxi_act,round,constrain,analyze,rank)
use general 
use evb_mod
use qmdff
use pbc_mod
implicit none
!
!     include MPI library
!
include 'mpif.h'
!  ambigious reference of z-coord to qmdff z(94) array!!
integer::i,j,k,l,istep   ! loop indices and actual step number
integer::rank   ! the current MPI rank
real(kind=8)::dt,dt_2
real(kind=8)::etot,epot,ekin,epot1
real(kind=8)::eksum
logical::analyze  ! if temperature etc shall be calculated
real(kind=8)::temp,pres
real(kind=8),dimension(natoms)::charge  !dummy array for charges
!   for umbrella sampling
real(kind=8),dimension(3,natoms)::derivs_1d,q_1b
real(kind=8),dimension(3,natoms,nbeads)::derivs
real(kind=8),dimension(nat6)::act_int,int_ideal
real(kind=8)::poly(4,nbeads)  
real(kind=8)::beta_n    ! the inverse temperature per bead
real(kind=8)::twown   ! inverse inverse temperature
real(kind=8)::pi_n    ! half circle section for each bead
real(kind=8)::wk,wt,wm  ! circle section scaled parameters 
real(kind=8)::cos_wt,sin_wt   ! sinus and cosinus value of the circle section
!real(kind=8)::pi   ! the pi
real(kind=8)::infinity,dbl_prec_var  ! test if one of the coordinates is infinity
real(kind=8)::p_new ! the actual new momentum/coordinate component
real(kind=8),allocatable::q_old(:,:,:) ! for hard boxes and RPMD: store actual coord
real(kind=8)::centroid(3,natoms)  !the centroid with all the com's of all beads
integer::round  ! the actual window to sample (for calculation of avg/variance)
integer::constrain  ! if the trajectory shall be constrained to the actual xi value
                    ! (only for recrossing calculations)
                    ! -1 : no constrain or umbrella at all (dynamic.x)
                    !  0 : usual umbrella sampling
                    !  1 : constrain to dividing surface
                    !  2 : child trajectory: no restrain at all
                    !  3 : for pre sampling with tri/tetramolecular reactions
integer::const_good  ! 0:success, 1:failure (in constrain_q)
integer::num,modul
!   for AFM simulation
real(kind=8)::move_act(3),move_shall(3)
real(kind=8)::ekin1,ekin2
real(kind=8)::afm_bias(3),afm_force,dist
!     for MPI parallelization
integer::ierr,impi_error
integer::psize,source,status(MPI_STATUS_SIZE)
!   for mirror plane simulations
real(kind=8)::pos_first,pos_second
!   for umbrella samplings 
real(kind=8)::xi_ideal,xi_real
real(kind=8)::dxi_act(3,natoms)
!   for kinetic energy 
real(kind=8)::act_temp
!   for Nose-Hoover thermostat
real(kind=8)::massvec(3,natoms,nbeads)  ! vector with atom masses per dimension
real(kind=8)::ndoft ! number degrees of freedom
real(kind=8)::nose_zeta2
integer::tries  ! number of periodic corrections before throwing an error
integer::j_lower,j_upper
!   for VASP coordinate input
real(kind=8)::coord_mat(3,3),coord_inv(3,3)
real(kind=8)::q_act_frac(3)
real(kind=8)::q_full_frac(3,natoms,nbeads)
logical::act_fix
!   for Nose-Hoover barostat
real(kind=8)::e2,e4,e6,e8
real(kind=8)::eterm2,term,term2,resize,expterm
real(kind=8)::stress_ten(3,3),ekin_ten(3,3),factor
real(kind=8)::act_vol,act_dens
!   for Berendsen barostat 
real(kind=8)::beren_scale
real(kind=8)::press_avg  ! the average pressure for printout
!   for printout of coordinate values 
real(kind=8)::coord_act
real(kind=8)::ang,dihed
!parameter(pi=3.1415926535897932384626433832795029d0)


infinity = HUGE(dbl_prec_var)   ! set the larges possible real value
!
!     Calculate new box dimensions, if the NPT ensemble is used
!
if (periodic) then
   volbox=boxlen_x*boxlen_y*boxlen_z
end if
!
!     Reset the virial tensor if a pressure is applied
!
if (npt) then
   calc_vir=.true.
   vir_ten=0.d0
end if

!
!     Write out the current trajectory frame, if iwrite 
!
!if (istep .gt. 997) then
!   write(*,*) q_i(:,101,1)*bohr,xi_real
!end if
!write(*,*) "use",use_calc_rate,rank
if (constrain .lt. 0) then
   if (use_calc_rate .or. rank .eq. 0) then
      if (mod(istep,iwrite) .eq. 0) then
         if (.not. output_sparse .or. nbeads .eq. 1) then
            write(28,*) natoms*nbeads
            if (npt) then

               write(28,*) boxlen_x*bohr,boxlen_y*bohr,boxlen_z*bohr
            else
               write(28,*)
            end if
            do k=1,nbeads
               do i=1,natoms
                  write(28,*) name(i),q_i(:,i,k)*bohr
               end do
            end do
            flush(28)
         end if
!
!     If the VASP formate is used for input, write new CONTCAR and new frame 
!       in XDATCAR files!
!
         if (coord_vasp  .and. nbeads .eq. 1) then
            coord_mat(:,1)=vasp_a_vec(:)
            coord_mat(:,2)=vasp_b_vec(:)
            coord_mat(:,3)=vasp_c_vec(:)

            call matinv3(coord_mat,coord_inv)

            open(unit=50,file="CONTCAR",status="replace")
            write(50,'(a,i9,a)') "CONTCAR (step ",istep,"), written by Caracal (dynamic.x)"
            write(50,*) vasp_scale
            write(50,*) vasp_a_vec
            write(50,*) vasp_b_vec
            write(50,*) vasp_c_vec
            do i=1,nelems_vasp
               write(50,'(a,a)',advance="no") " ",trim(vasp_names(i))
            end do
            write(50,*)
            write(50,*) vasp_numbers(1:nelems_vasp)
            if (vasp_selective) then
               write(50,*) "Selective dynamics"
            end if
            write(50,*) "Direct"
            if (xdat_first) then
               write(51,'(a,i9,a)') "step ",istep,", written by Caracal (dynamic.x)"
               write(51,*) vasp_scale
               write(51,*) vasp_a_vec
               write(51,*) vasp_b_vec
               write(51,*) vasp_c_vec
               do i=1,nelems_vasp
                  write(51,'(a,a)',advance="no") " ",trim(vasp_names(i))
               end do
               write(51,*)
               write(51,*) vasp_numbers(1:nelems_vasp)
            xdat_first=.false.
            end if
            write(51,'(a,i9)') "Direct step ",istep
!  
!     As in usual CONTCAR files, give the positions in direct coordinates!
!     convert them back from cartesians, by using the inverse matrix
!  
            do i=1,natoms
               q_act_frac=matmul(coord_inv,q_i(:,i,1)*bohr)
               if (vasp_selective) then

                  if (at_move(i)) then
                     write(50,*) q_act_frac(:),"   T   T   T "
                     write(51,*) q_act_frac(:)!,"   F   F   F "
                  else
                     write(50,*) q_act_frac(:),"   F   F   F "
                     write(51,*) q_act_frac(:)!,"   T   T   T "
                  end if
               else
                  write(50,*) q_act_frac(:)
                  write(51,*) q_act_frac(:)
               end if
            end do
            close(50)
            flush(51)
         end if
      end if
   end if
end if
!
!     For the first timestep: set the pressure to the desired value
!
if (npt) then
   if (istep .eq. 1) press_act=pressure
end if
if (use_calc_rate .or. rank .eq. 0) then
   do k=1,nbeads
      do i=1,natoms
         do j=1,3
            massvec(j,i,k)=mass(i)
         end do
      end do
   end do
end if
!
!     For the Nose-Hoover chain thermostat: apply the full chain on the momenta
!     on the half-time step
!
if (use_calc_rate .or. rank .eq. 0) then
   if (constrain .ne. 2) then
      if (thermostat .eq. 2) then
         call get_centroid(centroid)
!
!     The NPT ensemble: correct momenta and pressure 
!
         if (npt) then
            if (barostat .eq. 2) then
               call nhc_npt(dt,centroid,volbox)
            end if
!      call berendsen(dt,centroid,volbox)
!
!     The NVT ensemble: correct only momenta
!
         else 
            call nhc(dt,centroid)
         end if
      end if
   end if
end if

!
!     update the momentum (half time step)
!     ---> (Nose-Hoover step D: second part of half time momentum update)
!
!if (thermostat .eq. 1 .or. thermostat .eq. 2) then
if (use_calc_rate .or. rank .eq. 0) then
   p_i=p_i-0.5d0*dt*derivs
end if
!end if
!

!
!     set momentum to zero for fixed atoms 
!
if (use_calc_rate .or. rank .eq. 0) then
   do i=1,natoms
      if (.not. at_move(i)) then
         p_i(:,i,:)=0.d0
      end if
   end do
end if
!
!     calculate averaged kinetic energy for subset of the system
!
if (use_calc_rate .or. rank .eq. 0) then
   if (calc_ekin) then
!
!     first, calculate the kinetic energy by applying the centroid-momenta approximation
!
      ekin1=0d0
      do i=1,ekin_num
         k=ekin_atoms(i)
         do j=1,nbeads
       !  write(677,*) i,j,p_i(:,i,j),mass(i)
            ekin1=ekin1+dot_product(p_i(:,k,j),p_i(:,k,j))/(2d0*mass(k))
         end do
      end do
      ekin_avg=ekin_avg+ekin1
    
!
!     second, calculate the kinetic energy by applying the virial theorem
!
      ekin2=0d0
      beta_n = beta / nbeads
      do i=1,ekin_num
         k=ekin_atoms(i)
         do j=1,nbeads
         
            ekin2=ekin2+dot_product(q_i(:,k,j)-centroid(:,k),derivs(:,k,j))!+&
                    !  & mass(i)/beta_n**2/nbeads*((2*q_i(:,i,j)-q_i(:,i,j_upper)-q_i(:,i,j_lower))))

         end do
      end do
      ekin2=-ekin2/(2.d0*nbeads)
      ekin2=ekin2+ekin_num*3d0/(2d0*beta)
      ekin2_avg=ekin2_avg+ekin2
      write(156,*) ekin1/(nbeads**2*ekin_num),ekin2/ekin_num
   end if
end if
!
!     for AFM simulatons: fix only the anchor atom!
!
if (use_calc_rate .or. rank .eq. 0) then
   if (afm_run) then
      p_i(:,afm_fix_at,:)=0.d0
   end if 
end if
!
!     If a barostat is used: prepare the volume update
!
if (use_calc_rate .or. rank .eq. 0) then
   if (npt) then
      if (barostat .eq. 2) then
         term = vbar*0.5d0*dt
         term2 = term*term
         expterm=exp(term)
         eterm2 = expterm * expterm
         e2 = 1.0d0 / 6.0d0
         e4 = e2 / 20.0d0
         e6 = e4 / 42.0d0
         e8 = e6 / 72.0d0
         resize = 1.0d0 + term2*(e2+term2*(e4+term2*(e6+term2*e8)))
         resize = expterm * resize * dt
      end if
   end if
end if
!
!     update the positions: for one bead, do the usual verlet procedure
!
!     ---> Nose-Hoover step E
if (use_calc_rate .or. rank .eq. 0) then
   if (nbeads .eq. 1) then
!
!     If a system with hard box walls is simulated: check, if one of the atoms 
!     would move outside the box this timestep; if this is the case, revert 
!     the momentum in the respective axis
!
      if (box_walls) then
         do i=1,natoms
            do j=1,3
               if (q_i(1,i,1)+p_i(1,i,1)*dt/massvec(1,i,1) .ge. walldim_x) then
                  p_i(1,i,1)=-p_i(1,i,1)
               else if (q_i(2,i,1)+p_i(2,i,1)*dt/massvec(2,i,1) .ge. walldim_y) then
                  p_i(2,i,1)=-p_i(2,i,1)
               else if (q_i(3,i,1)+p_i(3,i,1)*dt/massvec(3,i,1) .ge. walldim_z) then
                  p_i(3,i,1)=-p_i(3,i,1)
               else if (q_i(j,i,1)+p_i(j,i,1)*dt/massvec(j,i,1) .le. 0.d0) then 
                  p_i(j,i,1)=-p_i(j,i,1)
               end if
            end do
         end do
      end if
!
!     If mirror planes were introduced into the system, reflect the respective 
!     atoms on them if they would fly through the plane in the next time step
!
      if (mirrors) then
         do i=1,mirror_num
            pos_first=q_i(mirror_dims(i),mirror_ats(i),1)
            pos_second=pos_first+p_i(mirror_dims(i),mirror_ats(i),1)*dt/massvec(mirror_dims(i), &
                        & mirror_ats(i),1)
            if (pos_first .gt. mirror_pos(i) .and. pos_second .lt. mirror_pos(i)) then
               p_i(mirror_dims(i),mirror_ats(i),1)=-p_i(mirror_dims(i),mirror_ats(i),1)
            end if      
            if (pos_first .lt. mirror_pos(i) .and. pos_second .gt. mirror_pos(i)) then
               p_i(mirror_dims(i),mirror_ats(i),1)=-p_i(mirror_dims(i),mirror_ats(i),1)
            end if
  
         end do
      end if
!
!    For NPT (Nose-Hoover) or NVT/NVE ensembles
!
  
      if (npt) then
         if (barostat .eq. 2) then
            q_i=q_i*eterm2+p_i/massvec*resize
         else 
            q_i=q_i+p_i*dt/massvec
         end if
      else 
 !     q_i=q_i+p_i*dt/massvec
      do i=1,3
         do j=1,natoms
            q_i(i,j,1)=q_i(i,j,1)+p_i(i,j,1)*dt/mass(j)
         end do
      end do

      end if
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
!     For hard boxes and mirror planes: store the actual coordinates, 
!       if momenta shall be reversed
      if (box_walls .or. mirrors) then
         allocate(q_old(3,natoms,nbeads))
         q_old=q_i
      end if

      do i = 1, 3
         do j = 1, Natoms
!
!     For periodic systems: asure that all replicas of one atoms 
!     fulfil the minimum image convention! No single replicas moved 
!     on the other side of the box, generating a giant distance..
!     --> possible violations of the box size will be corrected 
!      automatically in the next section!
!

!         if (periodic) then
!            tries=0
!            do k=2,nbeads          
!               do while (abs(q_i(i,j,1)-q_i(i,j,k)) .gt. box_len2)
!                  q_i(i,j,k)=q_i(i,j,k)+sign(box_len,q_i(i,j,1)-q_i(i,j,k))
!                  tries=tries+1
!                  if (tries .gt. 20) then
!                     write(*,*) "Too many correction steps needed for periodic dynamics!"
!                     write(*,*) "The system seems to be exploded! Check your settings!"
!                     call fatal
!                  end if
!               end do
!            end do
!         end if
            call rfft(p_i(i,j,:), Nbeads)
            call rfft(q_i(i,j,:), Nbeads)
         end do
      end do
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
!
!    For NPT ensemble: resize the box dimensions!
!
               if (npt) then
                  if (barostat .eq. 2) then
                     q_i(i,j,k) = p_i(i,j,k) * poly(3,k)*resize + q_i(i,j,k) *  &
                              & poly(4,k) * eterm2
                  else 
                     q_i(i,j,k) = p_i(i,j,k) * poly(3,k) + q_i(i,j,k) * poly(4,k)
                  end if
               else 
                  q_i(i,j,k) = p_i(i,j,k) * poly(3,k) + q_i(i,j,k) * poly(4,k)
               end if
               p_i(i,j,k) = p_new
            end do
         end do 
      end do
!
!     Transform back to Cartesian space
!
      do i = 1, 3
         do j = 1, Natoms
            call irfft(p_i(i,j,:), Nbeads)
            call irfft(q_i(i,j,:), Nbeads)
         end do
      end do
!
!     If a system with hard box walls is simulated: check, if one of the atoms 
!     has moved outside the box this timestep; if this is the case, revert 
!     the momentum in the respective axis
!
      if (box_walls) then
         do k=1,nbeads
            do i=1,natoms
               do j=1,3
                  if (q_i(1,i,1) .ge. walldim_x) then
                     q_i(1,i,k)=2.d0*q_old(1,i,k)-q_i(1,i,k)
                     p_i(1,i,k)=-p_i(1,i,k)
                  else if (q_i(2,i,1) .ge. walldim_y) then
                     q_i(2,i,k)=2.d0*q_old(2,i,k)-q_i(2,i,k)
                     p_i(2,i,k)=-p_i(2,i,k)
                  else if (q_i(3,i,1) .ge. walldim_z) then
                     q_i(3,i,k)=2.d0*q_old(3,i,k)-q_i(3,i,k)
                     p_i(3,i,k)=-p_i(3,i,k)
                  else if (q_i(j,i,1)  .le. 0.d0) then
                     q_i(j,i,k)=q_old(j,i,k)-2.d0*q_i(j,i,k)
                     p_i(j,i,k)=-p_i(j,i,k)
                  end if
               end do
            end do
         end do
         deallocate(q_old)
      end if
!
!     If mirror planes exist, revert the momentum of the respective atoms at 
!     the respective axis
!
!
!     If mirror planes were introduced into the system, reflect the respective 
!     atoms on them if they would fly through the plane in the next time step
!
      if (mirrors) then
         do i=1,mirror_num
            do k=1,nbeads
               pos_first=q_old(mirror_dims(i),mirror_ats(i),k)
               pos_second=q_i(mirror_dims(i),mirror_ats(i),k)
               if (pos_first .gt. mirror_pos(i) .and. pos_second .lt. mirror_pos(i)) then
                  q_i(mirror_dims(i),mirror_ats(i),k)=2.d0*q_old(mirror_dims(i),&
                            & mirror_ats(i),k)-q_i(mirror_dims(i),mirror_ats(i),k)
                  p_i(mirror_dims(i),mirror_ats(i),k)=-p_i(mirror_dims(i),mirror_ats(i),k)
               end if
               if (pos_first .lt. mirror_pos(i) .and. pos_second .gt. mirror_pos(i)) then
                  q_i(mirror_dims(i),mirror_ats(i),k)=2.d0*q_old(mirror_dims(i), &
                            & mirror_ats(i),k)-q_i(mirror_dims(i),mirror_ats(i),k)
                  p_i(mirror_dims(i),mirror_ats(i),k)=-p_i(mirror_dims(i),mirror_ats(i),k)
               end if
            end do
         end do
      end if


   end if
!
!     For Nose-Hoover barostat in NPT ensemble, align the dimensions of the box
!
   if (npt) then
      if (barostat .eq. 2) then
         boxlen_x=boxlen_x*eterm2
         boxlen_y=boxlen_y*eterm2
         boxlen_z=boxlen_z*eterm2
         boxlen_x2=boxlen_x*0.5d0
         boxlen_y2=boxlen_y*0.5d0
         boxlen_z2=boxlen_z*0.5d0
         volbox=boxlen_x*boxlen_y*boxlen_z
      end if
   end if

!
!     For periodic systems: shift all atoms that were moved outside the 
!     box on the other side!  If several RPMD beads are involved, change the 
!     coordinates of all of them in the same way
! 
!     Special case for VASP coordinates: correct for arbitrary cell shape
!     in fractional coordinates! (might be extended to all calculations in the
!     future)
!
   if (periodic) then
      if (coord_vasp) then
         coord_mat(:,1)=vasp_a_vec(:)
         coord_mat(:,2)=vasp_b_vec(:)
         coord_mat(:,3)=vasp_c_vec(:)

         call matinv3(coord_mat,coord_inv)

         do i=1,nbeads
            do j=1,natoms
               q_full_frac(:,j,i)=matmul(coord_inv,q_i(:,j,i)*bohr)
            end do
         end do
         do i=1,nbeads
            do j=1,natoms
               do k=1,3
                  tries=0
                  do while (q_full_frac(k,j,i) .lt. 0.d0) 
                     q_full_frac(k,j,:)=q_full_frac(k,j,:) + 1.d0
                     tries = tries+1
                     if (tries .gt. 20) then
                        write(*,*) "Too many correction steps needed for periodic dynamics! (lower)"
                        write(*,*) "The system seems to be exploded! Check your settings!"
                        call fatal
                     end if                
                  end do
                  do while (q_full_frac(k,j,i) .gt. 1.d0)
                     q_full_frac(k,j,:)=q_full_frac(k,j,:) - 1.d0
                     tries = tries +1
                     if (tries .gt. 20) then
                        write(*,*) "Too many correction steps needed for periodic dynamics! (upper)"
                        write(*,*) "The system seems to be exploded! Check your settings!"
                        call fatal
                     end if
                  end do 
               end do    
            end do
         end do
         do i=1,nbeads
            do j=1,natoms
               q_i(:,j,i)=matmul(coord_mat,q_full_frac(:,j,i)/bohr)
            end do
         end do

      else 
         do i=1,nbeads
            do j=1,natoms
               do k=1,3
                  tries=0
                  do while (q_i(k,j,i) .lt. 0) 
                     if (k .eq. 1) then
                        q_i(k,j,:)=q_i(k,j,:)+boxlen_x
                     else if (k .eq. 2) then
                        q_i(k,j,:)=q_i(k,j,:)+boxlen_y
                     else if (k .eq. 3) then
                        q_i(k,j,:)=q_i(k,j,:)+boxlen_z
                     end if
                     tries=tries+1
                     if (tries .gt. 100) then
                        write(*,*) "Too many correction steps needed for periodic dynamics! (lower)"
                        write(*,*) "The system seems to be exploded! Check your settings!"
                        call fatal
                     end if
                  end do
                  if (k .eq. 1) then
                     do while (q_i(1,j,i) .gt. boxlen_x)
                        q_i(1,j,:)=q_i(1,j,:)-boxlen_x
                        tries=tries+1
                        if (tries .gt. 100) then
                           write(*,*) "Too many correction steps needed for periodic dynamics! (x, upper)"
                           write(*,*) "The system seems to be exploded! Check your settings!"
                           call fatal
                        end if
                     end do
                  else if (k .eq. 2) then
                     do while (q_i(2,j,i) .gt. boxlen_y)
                        q_i(2,j,:)=q_i(2,j,:)-boxlen_y
                        tries=tries+1
                        if (tries .gt. 100) then
                           write(*,*) "Too many correction steps needed for periodic dynamics! (y,upper)"
                           write(*,*) "The system seems to be exploded! Check your settings!"
                           call fatal
                        end if
                     end do
                  else if (k .eq. 3) then
                     do while (q_i(3,j,i) .gt. boxlen_z)
                        q_i(3,j,:)=q_i(3,j,:)-boxlen_z
                        tries=tries+1
                        if (tries .gt. 100) then
                           write(*,*) "Too many correction steps needed for periodic dynamics! (z,upper)"
                           write(*,*) "The system seems to be exploded! Check your settings!"
                           call fatal
                        end if
                     end do
                  end if
               end do
            end do
         end do
      end if
   end if
!
!     Calculate the centroid (center of masses) for the 
!     ring polymer and use its positions for application 
!     of the force constant
!
   call get_centroid(centroid)
end if
!
!     Write out the current centroid trajectory frame for more than one bead, if iwrite 
!
if (constrain .lt. 0) then
   if ((mod(istep,iwrite) .eq. 0) .and. (nbeads .gt. 1)) then
      write(128,*) natoms
      if (npt) then
         write(128,*) boxlen_x*bohr,boxlen_y*bohr,boxlen_z*bohr
      else
         write(128,*)
      end if
      do i=1,natoms
         write(128,*) name(i),centroid(:,i)*bohr
      end do
!
!     If the VASP formate is used for input, write new CONTCAR and new frame 
!       in XDATCAR files!
!     Only write centroid if more than one bead is present!
!
      if (coord_vasp) then
         open(unit=50,file="CONTCAR",status="replace")
         write(50,*) "CONTCAR (step ",istep,"), written by Caracal (dynamic.x)"
         if (xdat_first) then
            write(51,*) "XDATCAR written by Caracal (dynamic.x)"
            write(51,*) vasp_scale
            write(51,*) vasp_a_vec
            write(51,*) vasp_b_vec
            write(51,*) vasp_c_vec
            do i=1,nelems_vasp
               write(51,'(a,a)',advance="no") " ",trim(vasp_names(i))
            end do
            write(51,*)
            write(51,*) vasp_numbers(1:nelems_vasp)
            xdat_first=.false.
         end if
         write(50,*) vasp_scale
         write(50,*) vasp_a_vec
         write(50,*) vasp_b_vec
         write(50,*) vasp_c_vec
         do i=1,nelems_vasp
           write(50,'(a,a)',advance="no") " ",trim(vasp_names(i))
         end do
         write(50,*)
         write(50,*) vasp_numbers(1:nelems_vasp)
         write(50,*) "Direct"
         write(51,*) "Direct  step ",istep
         if (vasp_selective) then
            write(50,*) "Selective dynamics"
         end if
!  
!     As in usual CONTCAR files, give the positions in direct coordinates!
!     convert them back from cartesians, by using the inverse matrix
!  
         do i=1,natoms
            q_act_frac=matmul(coord_inv,centroid(:,i)*bohr)
            if (vasp_selective) then
               act_fix=.false.
               do j=1,fix_num
                  if (fix_list(j) .eq. i) then
                     act_fix=.true.
                  end if
               end do
               if (act_fix) then
                  write(50,*) q_act_frac(:),"   F   F   F "
                  write(51,*) q_act_frac(:)!,"   F   F   F "
               else 
                  write(50,*) q_act_frac(:),"   T   T   T "
                  write(51,*) q_act_frac(:)!,"   T   T   T "
               end if
            else
               write(50,*) q_act_frac(:)
               write(51,*) q_act_frac(:)
            end if  
         end do
         close(50)
      end if
   end if
end if
!
!     constrain the system to the xi value if desired
!
if (use_calc_rate .or. rank .eq. 0) then
   const_good=0
   if (constrain .eq. 1) then
      call constrain_q(centroid,xi_ideal,dxi_act,const_good,dt)
   end if
   if (const_good .eq. 0) then
      epot=0.d0
   else
      write(*,*) "The SHAKE algorithm for constraint exeeded the maximum number of"
      write(*,*) "iterations! Restarting the trajectory! (energy penalty given)"
      epot=100000.d0
   end if
end if
!
!     get the potential energy and atomic forces
!     for each bead at once: define its current structure and collect 
!     all results in a global derivs array thereafter
!
if (energysplit) then
   e_cov_split=0.d0
   e_noncov_split=0.d0
end if
xi_test=xi_real  ! TEST TEST
if (.not. use_calc_rate) then
   call mpi_barrier(mpi_comm_world,ierr)
   call mpi_bcast(q_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(p_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
end if
do i=1,nbeads 
   q_1b=q_i(:,:,i)
   call gradient (q_1b,epot1,derivs_1d,i,rank)
   derivs(:,:,i)=derivs_1d
   epot=epot+epot1
end do
!
!   write out the gradients to file (verbose)
!
if (use_calc_rate .or. rank .eq. 0) then
   if (verbose) then 
      if (mod(istep,iwrite) .eq. 0) then
         write(29,*) natoms*nbeads
         write(29,*)
         do k=1,nbeads
            do i=1,natoms
               write(29,*) i,k,derivs(:,i,k)
            end do
         end do
      end if

!
!    write out the velocities to file (verbose)
!

      if (mod(istep,iwrite) .eq. 0) then
         write(30,*) natoms*nbeads
         write(30,*)
         do k=1,nbeads
            do i=1,natoms
               write(30,*) i,k,p_i(:,i,k)/mass(i)
            end do
         end do
      end if
   end if
end if
!
!     For Mechanochemistry calculations: Add the additional forces here 
!     to the gradient of the reference method!
!
if (use_calc_rate .or. rank .eq. 0) then
   if (add_force) then
      do i=1,nbeads
         derivs(:,force1_at,i)=derivs(:,force1_at,i)+force1_v(:)*force1_k
         derivs(:,force2_at,i)=derivs(:,force2_at,i)+force2_v(:)*force2_k
      end do
   end if

!
!     For an AFM simulation run: determine the current position of the 
!     bias force and its vector to be applied to the system (and in order 
!     to calculate the effective AFM force)
!
   if (afm_run) then
      move_act=centroid(:,afm_move_at)
      move_shall=afm_move_first(:)+afm_move_v(:)*afm_move_dist*(real(istep)/real(afm_steps))

      afm_bias=(move_act-move_shall)*afm_k

      afm_force=sqrt(afm_bias(1)**2+afm_bias(2)**2+afm_bias(3)**2)

      do i=1,nbeads
         derivs(:,afm_move_at,i)=derivs(:,afm_move_at,i)+afm_bias
      end do
      afm_force=afm_force/newton2au
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
!   write(*,*) "dxii",dxi_act
!   stop "Bogougvo"
!   write(*,*) "derivs_after",derivs,xi_real
   else if (constrain .eq. 3) then
      call umbrella(centroid,xi_ideal,int_ideal,xi_real,dxi_act,derivs,0)
   end if
end if
!
!     update the momentum (full time step)
!
!     ---> (Nose-Hoover step F: first part of full time momentum update)
!
if (use_calc_rate .or. rank .eq. 0) then
   p_i=p_i-0.5d0*dt*derivs
end if

!
!     constrain the momentum (zero for dxi) if desired
!
if (use_calc_rate .or. rank .eq. 0) then
   if (constrain .eq. 1) then
      call constrain_p(dxi_act)
   end if
end if

!
!     For the Nose-Hoover chain thermostat: apply the full chain on the momenta
!     on the full-time step
!
!call get_centroid(centroid)
if (use_calc_rate .or. rank .eq. 0) then
   if (constrain .ne. 2) then
      if (thermostat .eq. 2) then
!
!     The NPT ensemble: correct momenta and pressure 
!
         if (npt) then
            if (barostat .eq. 2) then
               call nhc_npt(dt,centroid,volbox)
            end if
!
!     The NVT ensemble: correct only momenta
!
         else
            call nhc(dt,centroid)
         end if
      end if
   end if
!
!     Calculate the kinetic energy part of the stress tensor
!
   if (npt) then
      ekin_ten=0.d0
      do i=1,nbeads
         do j=1,natoms
            do k=1,3
               do l=1,3
                  ekin_ten(k,l)=ekin_ten(k,l)*p_i(k,j,i)*p_i(l,j,i)/mass(j)
               end do
            end do
         end do
      end do
!
!     Now calculate the actual pressure from the virial tensor / stress tensor
!
      factor = 1.d0 / volbox
      do i=1,3
         do j=1,3
            stress_ten(j,i)=factor*(2d0*ekin_ten(j,i)-vir_ten(j,i))
         end do
      end do
      press_act=(stress_ten(1,1)+stress_ten(2,2)+stress_ten(3,3)) / 3.d0
   end if
end if
!
!     If the Berendsen barostat is chosen, rescale the box and the coordinates
!
if (use_calc_rate .or. rank .eq. 0) then
   if (npt) then
      if (barostat .eq. 1) then
         beren_scale=(1.d0+(dt*compress/taupres)*(press_act-pressure))**(1.d0/3.d0)
  !   write(*,*) "scale",beren_scale,(dt*compress/taupres),(press_act-pressure)
         boxlen_x=boxlen_x*beren_scale
         boxlen_y=boxlen_y*beren_scale
         boxlen_z=boxlen_z*beren_scale
         boxlen_x2=boxlen_x*0.5d0
         boxlen_y2=boxlen_y*0.5d0
         boxlen_z2=boxlen_z*0.5d0
         volbox=boxlen_x*boxlen_y*boxlen_z
         q_i=q_i*beren_scale
  
      end if
   end if
!
!     Apply andersen thermostat to apply an random hit and 
!     change the momentum (only every andersen_step steps!)
!     --> Replace momenta with a fresh sampling from a Gaussian
!     distribution at the temperature of interest
!
   if (constrain .ne. 2) then
      if (thermostat .eq. 1) then
         if (mod(istep,andersen_step) .eq. 0) then
            call andersen
         end if
      end if
   end if
end if
!
!     calculate temperature and other dynamical parameters for each dump step
!
if (use_calc_rate .or. rank .eq. 0) then
   if (npt) then
      press_avg=press_avg+press_act
   end if
   if (analyze) then
      ekin1=0.d0
      ekin=0.d0
      do i=1,natoms
         do j=1,nbeads
            if (at_move(i)) then
               ekin1=ekin1+dot_product(p_i(:,i,j),p_i(:,i,j))/(2d0*mass(i))
            end if
         end do
      end do
      ekin=ekin1
      if (vasp_selective) then
         act_temp=2d0*ekin1/3d0/0.316679D-5/(natoms-fix_num)/nbeads/nbeads
      else 
         act_temp=2d0*ekin1/3d0/0.316679D-5/natoms/nbeads/nbeads
      end if
      if (npt) then
         act_vol=volbox*bohr**3
         act_dens=(mass_tot*1.6605402E-24*emass)/(act_vol*1E-24)
      

         write(*,'(i8,f12.5,a,f10.4,a,f12.4,f14.5,a,f12.6)') istep,epot/nbeads,"     ",act_temp,"  ",act_vol, & 
               & press_avg*prescon/iwrite,"    ",act_dens
      else if (nvt) then 
         write(*,'(a,i8,a,f12.5,a,f10.4,a)') " Step: ",istep,"  --  pot. energy: ", &
                   & epot/nbeads," Hartree  --  temperature: ",act_temp," K"
      else 
         write(*,'(a,i8,a,f12.5,a)') " Step: ",istep,"  --  pot. energy: ", &
                & epot/nbeads," Hartree "

      end if
      if (.not. nve) then
         write(236,*) istep,act_temp
         temp_test=temp_test+act_temp
      end if
      press_avg=0.d0
   end if
end if
!
!     Calculate values of evaluated coordinates in actual structure 
!
if (use_calc_rate .or. rank .eq. 0) then
   if (eval_coord) then
      if (mod(istep,eval_step) .eq. 0) then
         write(141,'(i10)',advance="no") istep
         do i=1,eval_number
            if (eval_inds(i,1) .eq. 1) then
               write(141,'(3f14.8)',advance="no") centroid(:,eval_inds(i,2))*bohr
            else if (eval_inds(i,1) .eq. 2) then 
               write(141,'(f14.8)',advance="no") dist(eval_inds(i,2),eval_inds(i,3), &
                    & centroid)*bohr
            else if (eval_inds(i,1) .eq. 3) then
               write(141,'(f14.8)',advance="no") ang(eval_inds(i,2),eval_inds(i,3), &
                    & eval_inds(i,4),centroid)*180.d0/pi
            else if (eval_inds(i,1) .eq. 4) then
                write(141,'(f14.8)',advance="no") dihed(eval_inds(i,2),eval_inds(i,3), &
                   & eval_inds(i,4),eval_inds(i,5),centroid)*180.d0/pi
            end if
         end do
         write(141,*)
      end if
   end if
end if
!
!     If an error occured in the trajectory and one of the entries is either NaN
!     or Infinity, the error variable is set to 1
!
!do i=1,3
!   do j=1,natoms
!      do k=1,nbeads 
!         if (q_i(i,j,k) .ne. q_i(i,j,k) .or. q_i(i,j,k) .gt. infinity) then
!            traj_error=1
!         end if 
!      end do
!   end do
!end do
!write(*,*) "p_i",p_i
!write(*,*) "derivs",derivs
!
!     Check if the current positions are infinity of NaN
!
do i=1,nbeads
   do j=1,natoms
      do k=1,3
         if (q_i(k,j,i) .ne. q_i(k,j,i)) then
            if (use_calc_rate .or. rank .eq. 0) then
               write(*,*) "Something went wrong during the dynamics."
               write(*,*) "One of the spatial coordinates is NaN!" 
            end if
            call fatal
         end if
         if (q_i(k,j,i) .gt. infinity) then
            if (use_calc_rate .or. rank .eq. 0) then
               write(*,*) "Something went wrong during the dynamics."
               write(*,*) "One of the spatial coordinates is infinity!"
            end if
            call fatal
         end if
      end do
   end do
end do
!
!     Calculate new box dimensions, if the NPT ensemble is used
!
if (periodic) then
   volbox=boxlen_x*boxlen_y*boxlen_z
   if (volbox .ne. volbox) then
      if (use_calc_rate .or. rank .eq. 0) then
         write(*,*) "Something went wrong during the dynamics."
         write(*,*) "The volume of the simulation box is NaN!"
      end if
      call fatal
   end if
   if (volbox .gt. infinity) then
      if (use_calc_rate .or. rank .eq. 0) then
         write(*,*) "Something went wrong during the dynamics."
         write(*,*) "The volume of the simulation box is infinity!"
      end if
      call fatal
   end if
end if
!
!     Remove total translation and rotation of the system during the dynamics 
!      (especially needed for Nose-Hoover thermostat!)
!
if (constrain .lt. 0) then
!   if (.not. coord_vasp) then
      if (use_calc_rate .or. rank .eq. 0) then
         call transrot(centroid)
      end if
!   end if
end if
return
end subroutine verlet


