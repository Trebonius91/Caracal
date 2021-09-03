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
!     The subroutine opt_qmdff_ser optimizes the noncovalent parameters of a 
!     given QMDFF. It sets up the parameters 
!
!     The optimized parameters are written to [qmdffprefix].update
!     
!
!     part of EVB
!

subroutine opt_qmdff_ser(filegeo,fileenergy,points,fffile1)
use evb_mod
use qmdff
use general
use lm_module
implicit none
logical::has_next
integer::nat,nat3,i,j,k,l,p,m
integer::inc
integer::dg_evb_mode,mat_size
integer::points,mp
integer::coord_mode ! for init_int
! random number generation:
INTEGER SEED,t
INTEGER :: i_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
INTEGER, DIMENSION(1:8) :: dt_seed
real(kind=8)::r 
character(len=:),allocatable::char_num
character(len=30)::write_format
character(len=60)::filegeo,fileenergy
character(len=60)::fffile1
real(kind=8)::e1,e2,eref,V12,str_e1,str_e2
real(kind=8)::expo,d_p,root  ! short form for exp-part of the gaussians, dot-product
real(kind=8),allocatable::coord(:,:)
real(kind=8),dimension(3*natoms-6)::int_buf
real(kind=8),dimension(3*natoms)::g1,g2,gref
real(kind=8),dimension(3,natoms)::gqmdff1,gqmdff2
real(kind=8),dimension(points)::energies_ref,e_path,ff_e1,ff_e2,e_initial
real(kind=8),dimension(:,:),allocatable::path_xyz
real(kind=8),dimension(:),allocatable::int_path,int_coord  ! for internal coordinates
real(kind=8),dimension(:),allocatable::test_vec  ! TESTTEST
real(kind=8)::deldiscr,delsqrt,ediff,VdV,diffstep
integer::analytic  ! Test for EVB gradient calculation: Do it analytically
logical::bad_run    ! decides, if a levenberg marquardt run should be broken off 
! parameters for newton/LM optimizations
real(kind=8)::chi,chi_old,derivative,y1,y2,num_step,deltaE,func
real(kind=8)::s_rest,step
integer::maxstep,o,worse,INFO
real(kind=8),dimension(3*natoms-6)::g_ff1,g_ff2  ! for EVB-QMDFF gradients
real(kind=8),dimension(:,:,:),allocatable::xyz_path  ! xyz structures of the reactionpath
! for multi start local search
real(kind=8)::glob_chi,glob_chi_new
real(kind=8)::coord_test(3,7) 
real(kind=8)::dist,ang,oop,dihed 
real(kind=8),dimension(:),allocatable::rmsd_path_int 
real(kind=8)::minimum,temp  
real(kind=8),allocatable::temp3(:)  
integer,allocatable::temp2(:)  
integer::location,start,ihalf
integer,dimension(:,:),allocatable::coord_tmp,coord_reshape
integer::msls_step,line,readstat
integer,dimension(:),allocatable::coord_types
character(len=40)::coord_line
real(kind=8)::e_one
integer::par_num  ! number of QMDFF modification parameters 
real(kind=8),allocatable::par_vec(:)  ! vector of QMDFF modification parameters
real(kind=8),allocatable::par_opt(:)  ! vector with optimized QMDFF modification parameters
real(kind=8)::infinity,db1_prev_var  ! for check if any number is infinity
! for external Levenberg-Maquardt call
real(kind=8),allocatable::fvec(:),wa(:)
real(kind=8)::tol
integer::iflag,lwa
integer,allocatable::iwa(:)
character(len=60)::outfile_name
external::lm_qmdff ! the fitness function to be called in the Levenberg-Marquardt routine
integer::nat_mono  ! number of atoms of the monomer
integer::rank  ! MPI dummy parameter


rank=0
!
!     check, if any variable is infinity (define largest number)
!
infinity=HUGE(db1_prev_var)
!
!    Arrays for all needed input informations: energies (all_ens),
!    gradients (all_grads) and hessians (all_hess)
!    of the reference and the single QMDFFs (f1,f2)
!

!    short forms for variables
nat=natoms
nat3=3*nat
!
!    reduced dimension for internal variables
!
mp=dg_evb_points 
allocate(xyz_path(3,nat,points))
allocate(path_xyz(nat*3,points))
allocate(all_ens(dg_evb_points))
allocate(all_xyz(nat*3,dg_evb_points))
allocate(coord(3,natoms))
allocate(fvec(points))
!
!   allocate multi start local search parameters
!    First, read in the newly witten dimer QMDFF
!
call prepare("dimer.qmdff","dummy","dummy",1)
natoms=n_one
nat_mono=natoms/2  ! determine atom number of the monomer 
!
!     Open the trainset data and read them into arrays   
!
open (unit=126,file=fileenergy,status="old")
open (unit=127,file=filegeo,status='old')

do i=1,points
   read(126,*,iostat=readstat) energies_ref(i)
   if (readstat .ne. 0) then
      write(*,*) "The reference-energies file contains too few lines!"
      call fatal
   end if
   call next_geo(xyz_path(:,:,i),natoms,127,has_next)
   if (.not.has_next) exit
end do
allocate(xyz_ref_lm(3,natoms,points))
xyz_ref_lm=xyz_path/bohr
close(126)
close(127)

!
!    Calculate the initial QMDFF energies for comparison
!
chi=0.d0
par_num=1+7*nat_mono
allocate(mn_par(1+2*(par_num-1)))
mn_par=1.d0
do i=1,points

   call ff_eg(n_one,at,xyz_ref_lm(:,:,i),e_one,g_one)
   call ff_nonb(n_one,at,xyz_ref_lm(:,:,i),q,r0ab,zab,r094_mod,sr42,c6xy,e_one,g_one)
   call ff_hb(n_one,at,xyz_ref_lm(:,:,i),e_one,g_one)
   e_initial(i)=e_one
end do
e_initial=e_initial-(minval(e_initial)-minval(energies_ref))
do i=1,points
   chi=chi+(e_initial(i)-energies_ref(i))*(e_initial(i)-energies_ref(i))
end do

write(*,*) "The intial fitness is:",chi
chi=0.d0
!
!    determine number of parameters to optimize:
!    1 Coulomb (global scaling)
!    3*nat_mono Van der Waals (three parameters per atom)
!    nat_mono Hydrogen bond (one parameter per atom)
!
allocate(par_vec(par_num))
allocate(par_opt(par_num))
!
!   set important parameters
!
Lambda=lm_par
num_step=0.00001d0
glob_chi=0
msls_step=1
write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
write(*,*) "noncovalent QMDFF parameters are optimized with Multi Start Local Search."
write(*,*) "We will start",maxstart,"single optimization runs."

write(15,*) "--- QMDFF NONCOVALENT PARAMETERIZATION ---"
write(15,*)
write(15,*) "Calculation initiated at: " 
call timestamp ( )
write(15,*) 
write(15,*) "You have requested a Multi Start local search optimization"
write(15,*) "for noncovalent QMDFF parameters."
write(15,*)

write(15,*)

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
!   MULTI START LOCAL SEARCH: for all parameters individual values!
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
!   Write out global settings for coupling and optimization
!
write(*,*) "Multi Start Local Search algorithm:"
write(*,*) "---------------------------------------------------"

write(15,*) "Multi Start Local Search algorithm:"
write(15,*)
write(15,*) "SETTINGS:"
write(15,'(a,i5)') " * Number of atoms in the monomer: ", nat_mono
write(15,'(a,i5)') " * Number of optimizable parameters: ", par_num
write(15,'(a,i5)') " * Number of energy points in the trainset: ", points
write(15,'(a,i5)') " * Number of local searches ", maxstart
write(15,'(a,f10.5,f10.5)') " * Bounds for spawning of new params: ", lower_bond, upper_bond
write(15,*)
write(15,*) "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
!
!    define the variables in lm_module to use them for Levenberg-Marquardt
!
allocate(energies_ref_lm(points))
energies_ref_lm=energies_ref
!
!    define other Levenberg-Marquardt parameters
!
lwa=points*par_num+5*par_num+points
allocate(wa(lwa),iwa(par_num))
allocate(f_vec(points))
!
!    The tolerance value for optimizations
!
tol=0.00001D+00

do p=1,maxstart

!   Generate random start number
!
!   Generate random start numbers for all alphas
!   (taken from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html)
! ----- Set up random seed portably -----
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
!    Interval, in which the number should be
!
      r=r*(upper_bond-lower_bond)+lower_bond
      par_vec(i)=r
   end do
!   write(*,*) "hhh alpha manipulated!"
!   alph=40.d0

   write(15,'(a,i4)') " Multi start local search round ",p
!
!    Call the Levenberg-Marquardt routines from MINPACK
!
   call r8vec_print (par_num, par_vec, '  Initial parameter values:' )
   iflag = 1
!   tol = 0.000001D+00
   write(15,*)
   write(15,*) "Do Levenberg-Marquardt optimization..."
!
!     Fitness function is submitted via lm_function
!
   call lmdif1 ( lm_qmdff, points, par_num, par_vec, f_vec, tol, &
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
      write(15,*) "tolerance! Maybe try a new training set?"
   end if
   call r8vec_print ( par_num, par_vec, '  Optimized parameter values:' )
!
!    Evaluate the results of local optimization
!
!    TEST 25.10.2017: add penalty if wiggles are present..
!
!   if (extrema .gt. 3) then
   chi=sum(f_vec)
!   end if
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
write(15,*) ">>> The optimized alpha values are:" 
do i=1,par_num
   write ( 15, '(2x,i8,2x,g16.8)' ) i, par_opt(i)
end do
write(15,*)
write(15,*) ">>> File evb.info was written"
write(15,*) ">>> This file is needed for specification of the coupling term!"
write(15,*)
!
!     At least calculate the path energies for the optimized parameters 
!     and plot them with the inital energies for comparison
!
if (maxstart .gt. 1) then
   par_vec=par_opt
end if

par_opt=par_vec
chi=0.d0
open(unit=31,file="en_before_after.dat",status="unknown")
write(31,*) "# Reference energy    QMDFF energy (before)    QMDFF energy (after) "
mn_par(1)=par_vec(1)
mn_par(2:par_num)=par_vec(2:par_num)
mn_par(par_num+1:2*par_num)=par_vec(2:par_num)

do i=1,points

   call ff_eg(n_one,at,xyz_ref_lm(:,:,i),e_one,g_one)
   call ff_nonb(n_one,at,xyz_ref_lm(:,:,i),q,r0ab,zab,r094_mod,sr42,c6xy,e_one,g_one)
   call ff_hb(n_one,at,xyz_ref_lm(:,:,i),e_one,g_one)
   e_path(i)=e_one
end do
e_path=e_path-(minval(e_path)-minval(energies_ref_lm))
do i=1,points
   chi=chi+(e_path(i)-energies_ref_lm(i))*(e_path(i)-energies_ref_lm(i))
   write(31,*) energies_ref(i),e_initial(i),e_path(i)
end do
write(*,*)
write(*,*) "Energies of the reactionpath are written to energies.qmdff,"
write(*,*) "gradients are written to gradients.qmdff"
close(31)
!
!   write the par_opt reference file for further calculations
!
outfile_name=fffile1(1:len(trim(fffile1))-6) // "_mod.dat"
open(unit=67,file=outfile_name,status="unknown")
write(67,*) "** This file contains correction parameters for QMDFF nonbonded interactions **"
write(67,'(f14.8)') par_opt(1)
do i=1,natoms/2
   write(67,'(7f14.8)') par_opt(2+(i-1)*7:1+i*7)
end do
close(67)
write(*,*) 
write(*,*) trim(outfile_name)," file written. Use this file for further calculations."
write(*,*) "Exiting normally..."
write(*,*)

close(134)

!
!   deallocate arrays
!

return
end subroutine opt_qmdff_ser

