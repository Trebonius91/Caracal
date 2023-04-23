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
!     subroutine optimize_dq: do multi start local search with 
!     Levenberg-Marquardt to optimize parameters for different dQ
!     EVB coupling terms 
!
!     part of EVB
!
subroutine optimize_dq(energies_qmdff,nat,fileenergy,filegeo,filets,&
           datapoints,energies_result)
use evb_mod
use general
use lm_module
implicit none
real(kind=8),dimension(datapoints,2)::energies_qmdff 
real(kind=8),dimension(datapoints)::energies_reference,qmdff_shift,dq_dq,energies_result
real(kind=8),dimension(3,nat)::coord,coord_ts
real(kind=8),dimension(3*nat)::dq
real(kind=8),dimension(:),allocatable::par_vec  ! vector with optimizable parameters
real(kind=8),dimension(:),allocatable::beta_grad !the beta_vector(gradient)
real(kind=8),dimension(:,:),allocatable::alpha_hess !the alpha_matrix(hessian)
real(kind=8),dimension(:),allocatable::deltaA,deriv
real(kind=8)::chi,sigma,E1,E2,deltaQ,deltaE,ep1,denominator
real(kind=8)::chi_old   !Levenberg-Marquardt-Parameter
real(kind=8)::min_ref,min_qmdff,func
real(kind=8),dimension(4*nat-9)::int_coord
real(kind=8)::ep2,ep3,numerator
integer::int_mode
integer::par_num
character(len=40)::fileenergy,filegeo,filets 
integer::nat,datapoints,parameters,i,fileenergy_unit,info,maxstep,j,k,worse,m
logical::has_next
logical::bad_run    ! decides, if a levenberg marquardt run should be broken off 
! random number generation:
INTEGER SEED,t
INTEGER :: i_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
INTEGER, DIMENSION(1:8) :: dt_seed
real(kind=8)::r
real(kind=8)::infinity,db1_prev_var  ! for check if any number is infinity
! for multistart local search  
real(kind=8),dimension(:),allocatable::par_opt
real(kind=8)::glob_chi
integer::msls_step,p
! for external Levenberg-Maquardt call
real(kind=8),allocatable::fvec(:),wa(:)
real(kind=8)::tol
integer::iflag,lwa
integer,allocatable::iwa(:)
external::lm_de_func,lm_dq_func

!
!     check, if any variable is infinity (define largest number)
!
infinity=HUGE(db1_prev_var)

fileenergy_unit=20
!
!  Write infos to logfile and initialize number of parameters
!
write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
write(*,*) "EVB-dQ-parameters are optimized with Multi Start Local Search."
write(*,*) "We will start",maxstart,"single optimization runs."

write(15,*) "--- dQ-EVB COUPLING PARAMETERIZATION ---"
write(15,*)
write(15,*) "Calculation initiated at: "
call timestamp ( )
write(15,*)
write(15,*) "You have requested a Multi Start local search optimization"
write(15,*) "for dQ energy coupling parameters."
write(15,*)
!
!    Define the set of internal coordinates: Read in from coord_def.inp or 
!    define it on the fly
!
int_mode=1
call init_int(filegeo,datapoints,1,int_mode)


!   Write out global settings for coupling and optimization
!
write(*,*) "Multi Start Local Search algorithm:"
write(*,*) "---------------------------------------------------"
write(15,*) "Multi Start Local Search algorithm:"
write(15,*)
write(15,*) "SETTINGS:"
write(15,'(a,i5)') " * Number of atoms in the system: ", natoms
write(15,'(a,i5)') " * Number of used internal coordinates: ", nat6
write(15,'(a,i5)') " * Number of energy points on the reactionpath: ", datapoints
write(15,'(a,i5)') " * Number of local searches: ", maxstart
write(15,'(a,f10.5,f10.5)') " * Bounds for spawning of new parameters: ", lower_bond, upper_bond

!
!   Read in the Structure of the transition-state
!   and convert it to internal coordinates for later use 
!
has_next=.true.
allocate(int_ts(nat6))
open(unit=33,file=filets,status='old')
call next_geo(coord_ts,natoms,33,has_next)
call xyz_2int(coord_ts,int_ts,nat)
close(33)
!
!   Calculate the geometrical shifts for all structures of the 
!   reaction path
!
open(unit=fileenergy_unit,file=fileenergy,status='unknown')
do i=1,datapoints
   read(fileenergy_unit,*) energies_reference(i)
end do
close(fileenergy_unit)

has_next=.true.
open(unit=34,file=filegeo,status='old')
do m=1,datapoints
   call next_geo(coord,natoms,34,has_next)
!
!     convert to internal coordinates
!
   call xyz_2int(coord,int_coord,nat)
   if (.not.has_next) exit
!
!     The geometrical coordinate dq=q-q_ts
!  
   do i=1,nat6
      dq(i)=int_coord(i)-int_ts(i)
      if (i .gt. 2*nat-3) then
         dq(i)=sqrt((int_coord(i)-int_ts(i))*(int_coord(i)-&
               &int_ts(i))+(int_coord(i+nat-3)-int_ts(i+nat-3))*&
               &(int_coord(i+nat-3)-int_ts(i+nat-3)))
      end if
   end do
   dq_dq(m)=0
   do i=1,nat6
!
!     take the squared value
!
      dq_dq(m)=dq_dq(m)+(dq(i)*dq(i))
      
   end do
   dq_dq(m)=dsqrt(dq_dq(m))/10
end do
!
!     define path for Levenberg-Marquardt
!
allocate(dq_dq_lm(datapoints))
dq_dq_lm=dq_dq
write(15,'(a)') " * Used off-diagonal-Term:"
write(*,*) "Used off-diagonal-Term:"
if (off_basis .eq. "1g") then
   write(15,*) "  a*exp(-b*deltaQ^2)  (1g)"
   write(*,*) "a*exp(-b*deltaQ^2)  (1g)"
   par_num=2
else if (off_basis .eq. "3g") then
   write(15,*) "  a*exp(-b*deltaQ^2)+c*exp(-d*(deltaQ+m)^2)+d*exp(-f*(deltaQ+n)^2)   (3g)"
   write(*,*) "a*exp(-b*deltaQ^2)+c*exp(-d*(deltaQ+m)^2)+d*exp(-f*(deltaQ+n)^2)   (3g)"
   par_num=8
else if (off_basis .eq. "sd2") then
   write(15,*) "  a*exp(-b*deltaQ^2)+c*deltaQ^2*exp(-d*deltaQ^2)+"
   write(15,*) "        e*deltaQ^2*exp(-f*deltaQ^2)   (sd2)"
   write(*,*) "a*exp(-b*deltaQ^2)+c*deltaQ^2*exp(-d*deltaQ^2)+"
   write(*,*) "              e*deltaQ^2*exp(-f*deltaQ^2)   (sd2)"
   par_num=6
else
   write(15,*) "This coupling term isn´t implemented!"
   write(15,*) "Look into the manual for further information."
   write(*,*) "This coupling term isn´t implemented!"
   write(*,*) "Look into the manual for further information."
   call fatal
end if

write(15,'(a,i5)') " * Number of parameters to optimize: ", par_num
allocate(par_vec(par_num))
write(15,'(a)') " * Used internal coordinates:"
write(15,*)
do i=1,nat6
   if (coord_def(i,1) .eq. 1) then
      write(15,'(a,i5,i5)') "   bond: ", coord_def(i,2:3)
   else if (coord_def(i,1) .eq. 2) then
      write(15,'(a,i5,i5,i5)') "   bend angle: ", coord_def(i,2:4)
   else if (coord_def(i,1) .eq. 3) then
      write(15,'(a,i5,i5,i5,i5)') "   out of plane angle: ", coord_def(i,2:5)
   else if (coord_def(i,1) .eq. 4) then
      write(15,'(a,i5,i5,i5,i5)') "   dihedral angle: ", coord_def(i,2:5)
   end if
end do
write(15,*)
write(15,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"


write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
write(15,*) "Parameters to be optimized:"
!
!    define the variables in lm_module to use them for Levenberg-Marquardt
!
allocate(energies_ref_lm(datapoints))
energies_ref_lm=energies_reference
allocate(ff_e1_lm(datapoints),ff_e2_lm(datapoints))
ff_e1_lm=energies_qmdff(:,1)
ff_e2_lm=energies_qmdff(:,2)
allocate(fvec(datapoints))
!
!    define other Levenberg-Marquardt parameters
!
lwa=datapoints*par_num+5*par_num+datapoints
allocate(wa(lwa),iwa(par_num))
!
!    The tolerance value for optimizations
!
tol=0.00001D+00


!
!   Generate random start numbers for all parameters
!   (taken from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html)
!
glob_chi=0d0
sigma=1d0
msls_step=1
do p=1,maxstart
   bad_run=.false.
   do i=1,par_num
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
! ----- Done setting up random seed -----

      CALL RANDOM_NUMBER(r)
!
!  Interval, in which the number should be
!
      r=r*(upper_bond-lower_bond)+lower_bond
      par_vec(i)=r
   end do
   write(15,'(a,i4)') " Multi start local search round ",p
   call r8vec_print (par_num, par_vec, '  Initial parameter values:' )

   write(15,*)
   write(15,*) "Do Levenberg-Marquardt optimization..."
   call lmdif1 ( lm_dq_func, datapoints, par_num, par_vec, fvec, tol, &
         & info, iwa, wa, lwa )
   if (info.eq.0) then
      write(15,*) "FATAL ERROR: The input parameters for Levenberg-Marquardt are improper!"
      write(15,*) "Contact Julien Steffen to solve this problem!"
      call fatal
   else if ((info.eq.1) .or. (info.eq.2) .or. (info .eq.3)) then
      write(15,*) "The Levenberg-Marquardt run has converged!"
   else if (info.eq.4) then
      write(15,*) "ERROR: The alpha-vector is orthogonal to the Jacobian columns to machine"
      write(15,*) "precision in this round!  No convergence could be reached!"
      bad_run=.true.
   else if (info.eq.5) then
      write(15,*) "ERROR: Too many Levenberg-Marquardt runs were needed in this round!"
      write(15,*) "No convergence could be reached!"
      bad_run=.true.
   else if ((info.eq.6) .or. (info.eq.7)) then
      write(15,*) "Levenberg-Marquardt run has converged, but the errorsum was higher than"
      write(15,*) "tolerance! Maybe try another coupling term?"
   end if
   call r8vec_print (par_num, par_vec, '  Optimized parameter values:' )
!
!    Evaluate the results of local optimization
!
   chi=sum(fvec)
!
!    catch any exception in optimization runs!
!  
   if (chi .ne. chi) then
      write(*,*) "ERROR: The value of the errorsum is NaN! Please check your input!"
      call fatal
   else if (chi .gt. infinity) then
      write(*,*) "ERROR: The value of the errorsum is infinity! Please check your input!"
      call fatal
   end if

   write(15,*) "Error sum for this round:",chi
   write(15,*)
!
!     Multi start local search: only if optimization has converged: 
!     if the actual fitness value is better than the best before, replace 
!     it with the actual value
!
   if (.not. bad_run) then
       if (glob_chi .eq. 0) then
          glob_chi=chi
          par_opt=par_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization finished, new best fitness:",glob_chi
          write(15,*) ">> Local optimization finished!"
          write(15,*) ">> First best fitness value:",glob_chi
          write(15,*)
       else if (chi .le. glob_chi) then
          glob_chi=chi
          par_opt=par_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization finished, new best fitness:",glob_chi
          write(15,*) ">> Local optimization finished!"
          write(15,*) ">> New best fitness value:",glob_chi
          write(15,*)
       else
!             write(*,*) "Local optimization finished, the best fitness is still:",glob_chi
          write(15,*) ">> Local optimization finished!"
          write(15,*) ">> No better value found..."
          write(15,*)
       end if
   else
       if (glob_chi .eq. 0) then
!
!      avoid endings of optimizations without optimized fitness value
!
          glob_chi=chi
          par_opt=par_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization not converged, new best fitness:",glob_chi
          write(15,*) ">> Local optimization not converged..."
          write(15,*) ">> However, this is the first measured fitness value."
          write(15,*) ">> First best fitness value:",glob_chi
          write(15,*)
       else if (chi .le. glob_chi) then
          glob_chi=chi
          par_opt=par_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization not converged, new best fitness:",glob_chi
          write(15,*) ">> Local optimization not converged..."
          write(15,*) ">> New best fitness value:",glob_chi
          write(15,*)
      end if
   end if
   write(15,*) "-----------------------------------------------------------------"
   msls_step=msls_step+1
end do
write(*,*)  ".. .----. -- -.. --- -. . "
write(*,*) "Multi start local optimization finished!"
write(*,*) "The best fitness-value is:", glob_chi
write(*,*) "Look into evbopt.log for further details."

write(15,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
write(15,*)
write(15,*)  ".. .----. -- -.. --- -. . "
write(15,*)
write(15,*) ">>> Multi start local optimization finished!"
write(15,*) ">>> The best fitness-value is:", glob_chi
write(15,*) ">>> The optimized parameter values are:"
do i=1,par_num
   write ( 15, '(2x,i8,2x,g16.8)' ) i, par_vec(i)
end do
write(15,*)
write(15,*) ">>> File evb_pars.dat was written"
write(15,*) ">>> This file is needed for specification of the coupling term!"
write(15,*)

!
!     Write the evb.info file with the optimized coupling parameters
!
open(unit=67,file="evb_pars.dat",status="unknown")
write(67,*) "** This file contains parameters for a EVB-dQ-calculation (",trim(off_basis),"-function) **"
do i=1,par_num
   write(67,*) par_vec(i)
end do
close(67)
write(*,*)
write(*,*) "Energies of the reactionpath are written to energies.qmdff."
write(*,*)
write(*,*) "evb_pars.dat file written. Use this file for further calculations."
write(*,*) "Exiting normally..."
write(*,*)

!
!   At last, calculate the path energies 
!   for the optimized alpha and b_vec values
!
if (maxstart .gt. 1) then
   par_vec=par_opt
end if
offa=par_vec(1)
offb=par_vec(2)
if (off_basis .eq. "3g") then
   offc=par_vec(3)
   offd=par_vec(4)
   offe=par_vec(5)
   offf=par_vec(6)
   offm=par_vec(7)
   offn=par_vec(8)
else if (off_basis .eq. "sd2") then
   offc=par_vec(3)
   offd=par_vec(4)
   offe=par_vec(5)
   offf=par_vec(6)
end if

do i=1,datapoints
   E1=energies_qmdff(i,1)
   E2=energies_qmdff(i,2)
   deltaQ=dq_dq(i)
   deltaE=abs(E1-E2)
   if (off_basis .eq. "1g") then
      ep1=exp(-offb*deltaQ*deltaQ)
      denominator=sqrt(deltaE*deltaE+4*offa*offa*ep1*ep1)
      energies_result(i)=0.5d0*((E1+E2-denominator))
   else if (off_basis .eq. "3g") then
      ep1=exp(-offb*deltaQ*deltaQ)
      ep2=exp(-offd*(deltaQ+offm)*(deltaQ+offm))
      ep3=exp(-offf*(deltaQ+offn)*(deltaQ+offn))
      energies_result(i)=0.5*(E1+E2-sqrt(deltaE*deltaE+4*(offa*&
            &ep1+offc*ep2+offe*ep3)*(offa*ep1+offc*ep2+offe*ep3)))
   else if (off_basis .eq. "sd2") then
      ep1=exp(-offb*deltaQ*deltaQ)
      ep2=exp(-offd*deltaQ*deltaQ)
      ep3=exp(-offf*deltaQ*deltaQ)
      energies_result(i)=0.5d0*(E1+E2-sqrt(deltaE*deltaE+4*(offa*&
            &ep1+offc*deltaQ*deltaQ*ep2+offe*deltaQ*deltaQ*ep3)*&
            &(offa*ep1+offc*deltaQ*deltaQ*ep2+offe*deltaQ*deltaQ*ep3)))
   end if
end do

return
end subroutine optimize_dq
