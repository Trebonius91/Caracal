!
!     The subroutine dg_evb_init_par optimizes the unknown B_ijk-Parameters as well
!     as alpha_i-Parameters of the DG-EVB-calculation in a MPI parallel run (x cores).
!     It reads in all needed structures and energies, gradients and hessians,
!     sets up the DG-EVB coupling-function, solves the resulting set of linear 
!     equations and optimizes the variable alpha parameters with a multi start 
!     local search algorithm.
!     The optimized parameters are written into the evb.info-file, energies and gradients 
!     for the used structures are also plotted out 
!
!     All calculations are done using internal coordinates than can be chosen
!     by the user
!
!     part of EVB
!

subroutine dg_evb_init_par(dg_evb_mode,filegeo,fileenergy,points,psize,rank)
use evb_mod
use general
use lm_module
implicit none
!
!     include MPI library
!
include 'mpif.h'
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
integer::clock
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
!
!     for MPI parallelization
!
integer::rank,psize  ! MPI : index of processor and their number
integer::ierr !MPI
real(kind=8)::run_stat  ! exit status (possible errors) for actual run
real(kind=8)::act_chi ! fitness for actual slave
real(kind=8)::act_alphas(dg_evb_points+add_alph) ! alpha values of slave
real(kind=8),allocatable::act_vec(:)  ! DG-EVB prefactor coefficients of slave
real(kind=8),allocatable::message(:) ! vector for results
integer::count  ! actual round
integer::dest,schedule,status,tag_mpi  ! destination number
integer::loop_large,loop_rest  ! master loops
integer::mess_len ! number of elements in result message vector
integer::source,numwork
real(kind=8)::numwork_slave
integer::address  ! source of finished output from slaves
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
!     Arrays for all needed input informations: energies (all_ens),
!     gradients (all_grads) and hessians (all_hess)
!     of the reference and the single QMDFFs (f1,f2)
!

!     short forms for variables
nat=natoms
nat3=3*nat
!
!     reduced dimension for internal variables
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
   if (rank .eq. 0) then
      write(*,*) "Used DG-EVB mode: 1 - energies (E)"
   end if
   mat_size=dg_evb_points
else if (dg_evb_mode .eq. 2) then
   if (rank .eq. 0) then
      write(*,*) "Used DG-EVB mode: 2 - energies and gradients (E+G)"
   end if
   mat_size=dg_evb_points*(1+nat6)
else if (dg_evb_mode .eq. 3) then
   if (rank .eq. 0) then
      write(*,*) "Used DG-EVB mode: 3 - energies, gradients and hessians (E+G+H)"
   end if
   mat_size=dg_evb_points*(1+nat6+(nat6)*(nat6+1)/2)
end if
if (rank .eq. 0) then
   write(*,*) "The (leading) dimension of the linear equation F=DB is ",mat_size
end if
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
!   allocate MPI message parameters
!
mess_len=4+dg_evb_points+add_alph+mat_size
allocate(act_vec(mat_size),message(mess_len))

!
!  open the ref.input file with informations about the reference points
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
      if (rank .eq. 0) then
         write(*,*) "The reference-energies file contains too few lines!"
         call fatal
      end if
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
do i=1,dg_evb_points
!
!   Read in the reference structures for every fixpoint
!
   call next_geo(coord,natoms,124,has_next)

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
!    set up random seed for all processors
!
call random_init_local(rank)
!
!   set important parameters
!
Lambda=lm_par
num_step=0.00001d0
glob_chi=0
msls_step=1
!     define default tag value
tag_mpi=0
!     define default count value
count=1
!
!    write out some useful informations regarding the calculation: only for process 1
!
if (rank .eq. 0) then
   write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
   write(*,*) "dg_evb-parameters are optimized with Multi Start Local Search."
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
   write(*,*) "The calculation will be done in parallel!",psize-1,"slaves and one"
   write(*,*) " master are initiated."
   write(15,*) "The calculation will be done in parallel!",psize-1,"slaves and one" 
   write(15,*) " master are initiated."

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
   write(15,'(a,i5)') " * Number of processors used for MPI: ", psize
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
end if
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

call mpi_barrier(mpi_comm_world,ierr)
!
!     The master process (rank=0) it sends the actual child sampling tasks
!     to the worker that will do the samplings for a bunch of child 
!     trajectories
!
loop_large=int(maxstart/(psize-1))
loop_rest=mod(maxstart,psize-1)
msls_step=1  ! total number of multi start local search steps
glob_chi=0.d0
if (rank .eq. 0) then
!
!     at first, give each slave a task until all full rounds are finished
!
   do i=1,loop_large+1
      do j=1,psize-1
         schedule=i
         if (i .eq. loop_large+1) then
            if (j .gt. loop_rest) then
               schedule=-1
            else
              ! cycle
            end if
         end if
         dest=j
!
!     actual number of local search step (include number of processors!)
!         
         msls_step=j+(i-1)*(psize-1)
!
!     if desired, print output for start of single runs
!
         if (schedule .ne. -1) then
            if (more_info) then
               write(*,*) "Local optimization No.",msls_step,"started.."
            end if
            write(15,'(a,i5,a)') "**** Multi start local search round: ",msls_step,"****"
            write(15,*)
            write(15,*) "Do Levenberg-Marquardt optimization..."
         end if
!
!     Due to useless MPI bugs on Elgar, the actual schedule number needs to be written
!     out, else it cannot appear at the slave
!
         write(99,*) "sended schedule",schedule
         numwork_slave=real(schedule)        
 
         call mpi_send(numwork_slave, 1, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
      end do
!
!     recieve results from all slaves for the i'th round
!     if number of procs isnt a divisor of total workload, ask only remaining ones!
!
      do j=1,psize-1
         if ((i .lt. loop_large+1) .or. (i.eq. loop_large+1 .and. j .le. loop_rest)) then
            message=0.d0
            call mpi_recv(message, mess_len, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,tag_mpi, &
               & MPI_COMM_WORLD,status,ierr)
            bad_run=.false.
!
!     guarantee that the same msls_step number of both starting and finishing operations for 
!     the same thread are obtained
!
            address=int(message(dg_evb_points+add_alph+mat_size+4))
            msls_step=address+(i-1)*(psize-1)
!
!     extract single results from large results array: 1. actual error sum,
!       2. actual values for alphas, 3. actual values for b prefactors,
!       4. status (if a problem arose, which kind was it..)
!
            act_chi=message(1)
            act_alphas=message(2:dg_evb_points+add_alph+1)
            act_vec=message(dg_evb_points+add_alph+2:dg_evb_points+add_alph+3+mat_size)
            run_stat=message(dg_evb_points+add_alph+mat_size+4)
!
!     check possible outcomings of the current local optimization (from Levenberg-Marquardt routine)
!
            if (((run_stat - 1.d0) .le. 1E-10) .or. ((run_stat - 2.d0) .le. 1E-10) .or.&
                   &  ((run_stat - 3.d0) .le. 1E-10)) then 
               write(15,'(a,i5,a)') " The Levenberg-Marquardt run has converged in round",msls_step,"!"
            else if ((run_stat - 4.d0) .le. 1E-10) then
               write(15,*) "ERROR: The alpha-vector is orthogonal to the Jacobian columns to machine"
               write(15,'(a,i5,a)') " precision in round No.",msls_step,"!  No convergence could be reached!"
               bad_run=.true.
            else if ((run_stat - 5.d0) .le. 1E-10) then
               write(15,'(a,i5,a)') " ERROR: Too many Levenberg-Marquardt runs were needed in round No.",msls_step,"!"
               write(15,*) "No convergence could be reached!"
               bad_run=.true.
            else if (((run_stat - 6.d0) .le. 1E-10) .or. ((run_stat - 7.d0) .le. 1E-10)) then
               write(15,*) "Levenberg-Marquardt run has converged, but the errorsum was higher than"
               write(15,'(a,i5,a)') " tolerance in round No.",msls_step,"! Maybe try a new set of coordinates?"
            end if
!
!     general checks: if fitness is infinity or NaN
!
            if (act_chi .ne. act_chi) then
               write(*,*) "ERROR: The value of the errorsum is NaN! Please check your input!"
               call fatal
            else if (act_chi .gt. infinity) then
               write(*,*) "ERROR: The value of the errorsum is infinity! Please check your input!"
               call fatal
            end if
!
!     check if one of the alphas is lower than 1 (i.e. it diverges)
!
            do k=1,dg_evb_points+add_alph
               if (act_alphas(k) .lt. 1.0d0) then
                  write(15,*)
                  write(15,*) "OOPS! At least one alpha value is lower than 1! This means that the corresponding"
                  write(15,*) "gaussian tends to diverge at larger distances!"
                  write(15,*) "Therefore we add a penalty of 1.0 to the error sum and count this run as"
                  write(15,'(a,i5,a)') " not converged! (in round No.",msls_step,")."
                  write(15,*)
                  bad_run=.true.
                  act_chi=act_chi+1d0
                  exit
               end if
            end do
            write(15,*) "Error sum for round No.:",msls_step,"is:",act_chi
            write(15,*)
            if (more_info) then
               write(*,*) "After Multi-Start-Local-Search step:",msls_step
               write(*,*) "Fitness value:",act_chi
            end if

!
!     after error handling, check if a new best run was found or not
!
            if (.not. bad_run) then
               if (glob_chi .eq. 0) then
                  glob_chi=act_chi
                  alph_opt=act_alphas
                  b_opt=act_vec
                  write(*,*) "After Multi-Start-Local-Search step:",msls_step
                  write(*,*) "Local optimization finished, new best fitness:",glob_chi
                  write(15,*) ">> Local optimization finished!"
                  write(15,*) ">> First best fitness value:",glob_chi
                  write(15,*)
               else if (act_chi .le. glob_chi) then
                  glob_chi=act_chi
                  alph_opt=act_alphas
                  b_opt=act_vec
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
               if (glob_chi .eq. 0.d0) then
!
!      avoid endings of optimizations without optimized fitness value
!
                  glob_chi=act_chi
                  alph_opt=act_alphas
                  b_opt=act_vec
                  write(*,*) "After Multi-Start-Local-Search step:",msls_step
                  write(*,*) "Local optimization not converged, new best fitness:",glob_chi
                  write(15,*) ">> Local optimization not converged..."
                  write(15,*) ">> However, this is the first measured fitness value."
                  write(15,*) ">> First best fitness value:",glob_chi
                  write(15,*)
               else if (act_chi .le. glob_chi) then
                  glob_chi=act_chi
                  alph_opt=act_alphas
                  b_opt=act_vec
                  write(*,*) "After Multi-Start-Local-Search step:",msls_step
                  write(*,*) "Local optimization not converged, new best fitness:",glob_chi
                  write(15,*) ">> Local optimization not converged..."
                  write(15,*) ">> New best fitness value:",glob_chi
                  write(15,*)
               end if
            end if
            write(15,*) "-----------------------------------------------------------------"
 
         end if
      end do
   end do
!
!     Switch off the remaining processors
!

   do j=1,loop_rest
      dest=j
      call mpi_send(-1.d0, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
   end do
else 
   do   
      source=rank 
!
!     Recieve message with actual workload instruction
!   
      call mpi_recv(numwork_slave,1, MPI_DOUBLE_PRECISION, 0,tag_mpi,MPI_COMM_WORLD,status,ierr)

!      write(*,*) "slave",rank,"has recieved",numwork_slave
!
!     If no more optimizations are to do, exit with the current worker
!
      if ((numwork_slave+1.d0)  .lt. 1D-6) then
         exit
      end if
!  
!     start actual calculation!
!

!   Generate random start number
!
!   Generate random start numbers for all alphas
!   (taken from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html)
! ----- Set up random seed portably -----
      bad_run=.false.
      do i=1,dg_evb_points+add_alph
!
!    call new random number
!
         CALL RANDOM_NUMBER(r)
!
!    Interval, in which the number should be
!
         r=r*(upper_bond-lower_bond)+lower_bond
         alph(i)=r
      end do

!
!    Call the Levenberg-Marquardt routines from MINPACK
!
      call r8vec_print ( dg_evb_points+add_alph, alph, '  Initial alpha values:' )
      iflag = 1
    
!
!     Fitness function is submitted via lm_function
!
      call lmdif1 ( lm_function, points, dg_evb_points+add_alph, alph, fvec, tol, &
            & info, iwa, wa, lwa )
      if (info.eq.0) then
         write(15,*) "FATAL ERROR: The input parameters for Levenberg-Marquardt are improper!"
         write(15,*) "Contact Julien Steffen to solve this problem!"
         call fatal
      end if
      call r8vec_print ( dg_evb_points+add_alph, alph, '  Optimized alpha values:' )
!   
!    Evaluate the results of local optimization
!
      chi=sum(fvec)
!
!     fill the message with results for the master 
!
      message(1)=chi
      message(2:dg_evb_points+add_alph+1)=alph
      message(dg_evb_points+add_alph+2:dg_evb_points+add_alph+2+mat_size)=b_vec
      message(dg_evb_points+add_alph+2+mat_size+1)=real(info)
      message(dg_evb_points+add_alph+2+mat_size+2)=real(rank) ! current rank of worker
      call mpi_send(message, mess_len, MPI_DOUBLE_PRECISION,0,tag_mpi,MPI_COMM_WORLD,ierr)
   end do
end if
call mpi_barrier(mpi_comm_world,ierr)
!
!     Write final informations and calculate path energies for the best parameters
!
if (rank .eq. 0) then
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
   write(15,*) ">>> File evb.info was written"
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
   open(unit=67,file="evb.info",status="unknown")
   write(67,*) "** This file contains parameters for a DG-EVB-calculation **"
   do i=1,mat_size
      write(67,*) b_vec(i)
   end do
   do i=1,dg_evb_points+add_alph
      write(67,*) alph(i)
   end do
   close(67)
   write(*,*) 
   write(*,*) "evb.info file written. Use this file for further calculations."
   write(*,*) "Exiting normally..."
   write(*,*)

   close(134)
end if 
call mpi_barrier(mpi_comm_world,ierr)

!
!   deallocate arrays
!
deallocate(all_ens,allf1_ens,allf2_ens,all_xyz,path_xyz,f_vec,b_vec,d_mat)
deallocate(alph_opt,b_opt)
!
!   remove useless fortran output file which is needed for MPI to work on Elgar
!
if (rank .eq. 0) then
   call system("rm fort.99")
end if
return
end subroutine dg_evb_init_par

