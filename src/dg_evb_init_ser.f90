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
!     The subroutine dg_evb_init_ser optimizes the unknown B_ijk-Parameters as well
!     as alpha_i-Parameters of the DG-EVB-calculation in a serial run (1 core).
!     It reads in all needed structures and energies, gradients and hessians,
!     sets up the DG-EVB coupling-function, solves the resulting set of linear 
!     equations and optimizes the variable alpha parameters with a multi start 
!     local search algorithm.
!     The optimized parameters are written into the evb_pars.dat-file, energies and gradients 
!     for the used structures are also plotted out 
!
!     All calculations are done using internal coordinates than can be chosen
!     by the user
!
!     part of EVB
!

subroutine dg_evb_init_ser(dg_evb_mode,filegeo,fileenergy,points)
use evb_mod
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
real(kind=8)::e1,e2,eref,V12,str_e1,str_e2
real(kind=8)::expo,d_p,root  ! short form for exp-part of the gaussians, dot-product
real(kind=8),dimension(3,natoms)::xyz2,coord
real(kind=8),dimension(3*natoms-6)::int_buf
real(kind=8),dimension(3*natoms)::g1,g2,gref
real(kind=8),dimension(3,natoms)::gqmdff1,gqmdff2
real(kind=8),dimension(3*natoms,3*natoms)::hess1,hess2
real(kind=8),dimension(dg_evb_points+add_alph)::alph
real(kind=8),dimension(points)::energies_ref,e_path,ff_e1,ff_e2
real(kind=8),dimension(1000)::e_diff
real(kind=8),dimension(:,:),allocatable::path_int
real(kind=8),dimension(:,:),allocatable::path_xyz
real(kind=8),dimension(:),allocatable::int_path,int_coord  ! for internal coordinates
real(kind=8),dimension(:),allocatable::test_vec  ! TESTTEST
real(kind=8)::deldiscr,delsqrt,ediff,VdV,diffstep
real(kind=8)::V_upper,V_lower,rootV,root2,off4
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
real(kind=8)::infinity,db1_prev_var  ! for check if any number is infinity
! for external Levenberg-Maquardt call
real(kind=8),allocatable::fvec(:),wa(:)
real(kind=8)::tol
integer::iflag,lwa
integer,allocatable::iwa(:)
external::lm_function
integer::rank  ! MPI dummy parameter

rank=0
!
!     check, if any variable is infinity (define largest number)
!
infinity=HUGE(db1_prev_var)
!
!     first, define the used set of internal coordinates!
!
coord_mode=2
if (read_coord) coord_mode=1
call init_int(filegeo,points,rank,coord_mode)
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
allocate(allf1_ens(dg_evb_points))
allocate(allf2_ens(dg_evb_points))
allocate(all_xyz(nat*3,dg_evb_points))
allocate(all_int(nat6,dg_evb_points))
allocate(path_int(nat6,points))
allocate(int_path(nat6))
allocate(point_int(nat6,dg_evb_points))
if (dg_evb_mode .ge. 2) then
   allocate(all_grad(nat6,points))
   allocate(allf1_grad(nat6,points))
   allocate(allf2_grad(nat6,points))
   if (dg_evb_mode .eq. 3) then
      allocate(all_hess(nat6,nat6,points))
      allocate(allf1_hess(nat6,nat6,points))
      allocate(allf2_hess(nat6,nat6,points))
   end if
end if 
!
!   allocate equation parameters
!
if (dg_evb_mode .eq. 1) then
   write(*,*) "Used DG-EVB mode: 1 - energies (E)"
   mat_size=dg_evb_points
else if (dg_evb_mode .eq. 2) then
   write(*,*) "Used DG-EVB mode: 2 - energies and gradients (E+G)"
   mat_size=dg_evb_points*(1+nat6)
else if (dg_evb_mode .eq. 3) then
   write(*,*) "Used DG-EVB mode: 3 - energies, gradients and hessians (E+G+H)"
   mat_size=dg_evb_points*(1+nat6+(nat6)*(nat6+1)/2)
end if
write(*,*) "The dimension of the linear equation F=DB is ",mat_size

allocate(f_vec(mat_size))
allocate(test_vec(mat_size)) ! TESTEST
allocate(b_vec(mat_size))

allocate(d_mat(mat_size,mat_size))
!
!   allocate multi start local search parameters
!
allocate(alph_opt(dg_evb_points+add_alph))
allocate(b_opt(mat_size))
!
!  open the grad_hess.dat file with informations about the reference points
!  read in the lines and convert geometries, gradients and hessians to 
!  internal coordinates
!
call read2int(dg_evb_mode)
open (unit=126,file=fileenergy,status="old")
open (unit=127,file=filegeo,status='old')
e_diff=0
!
!  calculate energies for the reaction path and convert all strucutres 
!  to internal coordinates
!

do i=1,points
   read(126,*,iostat=readstat) energies_ref(i)
   if (readstat .ne. 0) then
      write(*,*) "The reference-energies file contains too few lines!"
      call fatal
   end if
   call next_geo(coord,natoms,127,has_next)
   call eqmdff(coord,str_e1,str_e2)
   xyz_path(:,:,i)=coord
   ff_e1(i)=str_e1
   ff_e2(i)=str_e2
   call xyz_2int(coord,int_path,nat)
   if (.not.has_next) exit
   do j=1,nat6
      path_int(j,i)=int_path(j)
   end do
end do

close(126)
write(*,*) "##################  QMDFFs!"
do i=1,dg_evb_points
!
!   Read in the reference structures for every fixpoint
!
   call next_geo(coord,natoms,127,has_next)

   do j=1,natoms
      do k=1,3
         xyz2(k,j)=all_xyz((j-1)*3+k,i)*bohr
      end do
   end do
!
!   Calculate QMDFF energies of all reference-strucutures
!
   call eqmdff(xyz2,e1,e2)
   allf1_ens(i)=e1
   allf2_ens(i)=e2

!
!   Calculate QMDFF gradients and convert them to internal coordinates
!
   if (dg_evb_mode .ge. 2) then
      call egqmdff(xyz2,e1,e2,gqmdff1,gqmdff2)
      do j=1,nat
         do k=1,3
            g1((j-1)*3+k)=gqmdff1(k,j)
            g2((j-1)*3+k)=gqmdff2(k,j)
         end do
      end do
!
!   buffer array for internal geometry: int_buf
!
      call grad2int(xyz2,all_int(:,i),allf1_grad(:,i),g1)
      call grad2int(xyz2,all_int(:,i),allf2_grad(:,i),g2)
   end if
!
!   Calculate QMDFF hessians and convert them to internal coordinates
!
   if (dg_evb_mode .eq. 3) then
      call hessqmdff(xyz2,hess1,hess2)
      call hess2int(xyz2,all_int(:,i),allf1_hess(:,:,i),hess1,allf1_grad(:,i)) 
      call hess2int(xyz2,all_int(:,i),allf2_hess(:,:,i),hess2,allf2_grad(:,i))    
   end if
end do
point_int=all_int
!
!   Now set up the system of linear equations which shall be solved
!   Depending on which mode was chosen
!

!   LEFT SIDE (F-column-vector)
!   contained elements:
!   MODE1:
!   (V_12^2(q1), .., V_12^2(qn))
!   MODE2:
!   (V_12^2(q1), ..,V_12^2(qn),dV12^2(q1)/dq1, ..,dV12^2(qn)/dqn)
!   MODE3:
!   (V_12^2(q1), ..,V_12^2(qn),dV12^2(q1)/dq1, ..,dV12^2(qn)/dqn),
!   d^2V12(q1)/d^2q1, ..,d^2V12(qn)/d^2qn)
!
do i=1,dg_evb_points
   f_vec(i)=(allf1_ens(i)-all_ens(i))*(allf2_ens(i)-all_ens(i))
end do
if (dg_evb_mode .ge. 2) then
   do i=1,dg_evb_points
      inc=(i-1)*(nat6)+dg_evb_points
      do j=1,nat6
         f_vec(inc+j)=((allf1_grad(j,i)-all_grad(j,i))*&
                     &(allf2_ens(i)-all_ens(i)))+((allf1_ens(i)-all_ens(i))*&
                     &(allf2_grad(j,i)-all_grad(j,i)))
      end do
   end do
end if
if (dg_evb_mode .eq. 3) then
   do i=1,dg_evb_points
      line=0
      inc=(1+nat6)*mp+(i-1)*((nat6)*(nat6+1)/2)
      do j=1,nat6
         do k=j,nat6
            line=line+1
            f_vec(inc+line)=(allf2_grad(k,i)-all_grad(k,i))*(allf1_grad(j,i)-all_grad(j,i))+&
                        &(allf1_grad(k,i)-all_grad(k,i))*(allf2_grad(j,i)-all_grad(j,i))+&
                        &(allf2_ens(i)-all_ens(i))*(allf1_hess(j,k,i)-all_hess(j,k,i))+&
                        &(allf1_ens(i)-all_ens(i))*(allf2_hess(j,k,i)-all_hess(j,k,i))
         end do
      end do
   end do
end if

!
!   set important parameters
!
Lambda=lm_par
num_step=0.00001d0
glob_chi=0
msls_step=1
write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
write(*,*) "DG-EVB-parameters are optimized with Multi Start Local Search."
write(*,*) "We will start",maxstart,"single optimization runs."

write(15,*) "--- DG-EVB COUPLING PARAMETERIZATION ---"
write(15,*)
write(15,*) "Calculation initiated at: " 
call timestamp ( )
write(15,*) 
write(15,*) "You have requested a Multi Start local search optimization"
write(15,*) "for dg_evb alpha-parameters."
write(15,*)

if (read_coord) then
   write(15,*) "The internal coordinates used for DG-EVB coupling were read in from file."
else 
   write(15,*) "The internal coordinates used for DG-EVB coupling were found by"
   write(15,*) " analysis of the reaction path."
end if
if (double_alpha .and. dg_evb_mode.eq.3) then
   write(15,*)
   write(15,*) "You have requested DOUBLE_ALPHA, so now seperate alphas"
   write(15,*) "for the g(i,j)-gaussians are taken!"
   write(*,*) "You have requested DOUBLE_ALPHA, so now seperate alphas"
   write(*,*) "for the g(i,j)-gaussians are taken!"
end if
if (dg_mat_num) then
   write(*,*) "The coefficient matrix D for DG-EVB will be calculated numerically!"
else  
   write(*,*) "The coefficient matrix D for DG-EVB will be calculated analytically!"
end if
write(15,*)

! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
! ----------------------------------------------------------------------------------
!   MULTI START LOCAL SEARCH: for all alphas INDIVIDUAL values
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
write(15,'(a,i5)') " * Number of atoms in the system: ", nat
write(15,'(a,i5)') " * Number of used internal coordinates: ", nat6
if (dg_evb_mode .eq. 1) then
   write(15,'(a)') " * Used DG-EVB mode: 1 (energies)"
else if (dg_evb_mode .eq. 2) then
   write(15,'(a)') " * Used DG-EVB mode: 2 (energies + gradients)"
else 
   write(15,'(a)') " * Used DG-EVB mode: 3 (energies + gradients + hessians)"
end if
write(15,'(a,i5)') " * Number of DG-EVB points: ", dg_evb_points
write(15,'(a,i5)') " * Number of alpha coefficients: ", dg_evb_points+add_alph
write(15,'(a,i5)') " * Dimension of the D coefficient matrix: ", mat_size
write(15,'(a,i5)') " * Number of energy points on the reactionpath: ", points
write(15,'(a,i5)') " * Number of local searches ", maxstart
write(15,'(a,f10.5,f10.5)') " * Bounds for spawning of new alphas: ", lower_bond, upper_bond
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

!
!    define the variables in lm_module to use them for Levenberg-Marquardt
!
allocate(energies_ref_lm(points))
energies_ref_lm=energies_ref
dg_evb_mode_lm=dg_evb_mode
mat_size_lm=mat_size
allocate(ff_e1_lm(points),ff_e2_lm(points))
ff_e1_lm=ff_e1
ff_e2_lm=ff_e2
allocate(path_int_lm(nat6,points),fvec(points))
path_int_lm=path_int
!
!    define other Levenberg-Marquardt parameters
!
lwa=points*(dg_evb_points+add_alph)+5*(dg_evb_points+add_alph)+points
allocate(wa(lwa),iwa(dg_evb_points+add_alph))
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
   do i=1,dg_evb_points+add_alph
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
      alph(i)=r
   end do
!   write(*,*) "hhh alpha manipulated!"
!   alph=40.d0


   write(15,'(a,i4)') " Multi start local search round ",p
!
!    Call the Levenberg-Marquardt routines from MINPACK
!
   call r8vec_print ( dg_evb_points+add_alph, alph, '  Initial alpha values:' )
   iflag = 1
    
!   tol = 0.000001D+00
   write(15,*)
   write(15,*) "Do Levenberg-Marquardt optimization..."
!
!     Fitness function is submitted via lm_function
!
   call lmdif1 ( lm_function, points, dg_evb_points+add_alph, alph, fvec, tol, &
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
      write(15,*) "tolerance! Maybe try a new set of coordinates?"
   end if
   call r8vec_print ( dg_evb_points+add_alph, alph, '  Optimized alpha values:' )
!
!    Evaluate the results of local optimization
!
!    TEST 25.10.2017: add penalty if wiggles are present..
!
!   if (extrema .gt. 3) then
      chi=sum(fvec)!+(extrema-3)*0.1
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

   do i=1,dg_evb_points+add_alph
      if (alph(i) .le. 0.0d0) then
         write(15,*)
         write(15,*) "OOPS! At least one alpha value is lower than 0! This means that the corresponding"
         write(15,*) "gaussian tends to diverge at larger distances!"
         write(15,*) "Therefore we add a penalty of 1.0 to the error sum and count this run as"
         write(15,*) "not converged!"
         write(15,*)
         bad_run=.true.
         chi=chi+1d0
         exit
      end if
   end do
   write(15,*) "Error sum for this round:",chi
   write(15,*)
   if (.not. bad_run) then
       if (glob_chi .eq. 0) then
          glob_chi=chi
          alph_opt=alph
          b_opt=b_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization finished, new best fitness:",glob_chi
          write(15,*) ">> Local optimization finished!"
          write(15,*) ">> First best fitness value:",glob_chi
          write(15,*) 
       else if (chi .le. glob_chi) then
          glob_chi=chi
          alph_opt=alph
          b_opt=b_vec
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
          alph_opt=alph
          b_opt=b_vec
          write(*,*) "After Multi-Start-Local-Search step:",msls_step
          write(*,*) "Local optimization not converged, new best fitness:",glob_chi
          write(15,*) ">> Local optimization not converged..."
          write(15,*) ">> However, this is the first measured fitness value."
          write(15,*) ">> First best fitness value:",glob_chi
          write(15,*) 
       else if (chi .le. glob_chi) then
          glob_chi=chi
          alph_opt=alph
          b_opt=b_vec
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
do i=1,dg_evb_points+add_alph
   write ( 15, '(2x,i8,2x,g16.8)' ) i, alph_opt(i)
end do
write(15,*)
write(15,*) ">>> File evb_pars.dat was written"
write(15,*) ">>> This file is needed for specification of the coupling term!"
write(15,*)

!
!   At last, calculate the path energies AND GRADIENTS (!) 
!   for the optimized alpha and b_vec values
!
if (maxstart .gt. 1) then
   alph=alph_opt
   b_vec=b_opt
end if
123 open(unit=33,file="gradients.qmdff",status="unknown")

!
!    TEST TEST 12.07.2017 TEST TEST
!
alph_opt=alph
!
!  Write out energies and gradients for every structure 
!
call build_dmat(dg_evb_mode,mat_size,alph)
!
!  Do the matrix diagonalization to calculate b_vec parameters
!
call mat_diag(mat_size)

open(unit=31,file="energies.qmdff",status="unknown")

do i=1,points
   write(33,*) "structure No.",i
   write(33,*)

   V12=0
   call sum_v12(dg_evb_mode,mat_size,alph,path_int(:,i),V12)
   str_e1=ff_e1(i)
   str_e2=ff_e2(i)
   ediff=str_e1-str_e2
   root=(0.5d0*(str_e1-str_e2))*(0.5d0*(str_e1-str_e2))+V12
   if (root .le. 0) then
      e_path(i)=0.5*(str_e1+str_e2)
   else
      e_path(i)=0.5*(str_e1+str_e2)-sqrt(root)
   end if
   write(31,*)e_path(i) 
   chi=chi+(e_path(i)-energies_ref(i))*(e_path(i)-energies_ref(i))

end do
write(*,*)
write(*,*) "Energies of the reactionpath are written to energies.qmdff,"
write(*,*) "gradients are written to gradients.qmdff"
close(33)
close(31)
!
!   write the dg_evb output file for further calculations
!
open(unit=67,file="evb_pars.dat",status="unknown")
write(67,*) "** This file contains parameters for a DG-EVB-calculation **"
do i=1,mat_size
   write(67,*) b_vec(i)
end do
do i=1,dg_evb_points+add_alph
   write(67,*) alph(i)
end do
close(67)
write(*,*) 
write(*,*) "evb_pars.dat file written. Use this file for further calculations."
write(*,*) "Exiting normally..."
write(*,*)

close(134)

!
!   deallocate arrays
!
deallocate(all_ens,allf1_ens,allf2_ens,all_xyz,path_xyz,f_vec,b_vec,d_mat)
deallocate(alph_opt,b_opt)

return
end subroutine dg_evb_init_ser

