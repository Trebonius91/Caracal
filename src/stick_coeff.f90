!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CARACAL - Ring polymer molecular dynamics and rate constant calculations
!             on black-box generated potential energy surfaces
!
!   Copyright (c) 2025 by Julien Steffen (mail@j-steffen.org)
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
!     #################################################################
!     ##                                                             ##
!     ##  program stick_coeff -- sticking coefficients on surfaces   ##
!     ##                                                             ##
!     #################################################################
!
!     "stick_coeff" calculates the sticking coefficient of an arbitrary 
!     diatomic molecule on an arbitrary surface
!     
!
program stick_coeff
use general
use evb_mod
use omp_lib
use pbc_mod

implicit none
!
!     include MPI library
!
include 'mpif.h'
integer i,istep,nstep,asciinum,j,k,nat,l,m
integer mode,next,check,readstatus
integer::qmdffnumber,idum
logical::analyze  ! if temperature etc shall be calculated during verlet step
integer parameters!natoms,parameters  !mod
real*8 dt,dtdump,lengthUnits,energyUnits
real(kind=8),dimension(:),allocatable::int_coord  ! for internal coordinates
real(kind=8),dimension(:,:),allocatable::xyz2
real(kind=8),dimension(:,:),allocatable::ts_coordinates_a,ts_coordinates2_a
logical exist,query
real(kind=8)::etot,epot,ekin  ! the total energy for testing
real(kind=8),allocatable::centroid(:,:)  ! the RPMD centroid
real(kind=8)::dist,dist_first,dist_last
real(kind=8),allocatable::force_avg(:)
real(kind=8),allocatable::force_avgs(:)
real(kind=8),allocatable::len_avg(:),len_avgs(:)
integer::afm_avg
integer::afm_segment
integer::afm_segment_avg
logical::afm_second
character(len=120)::xyzfile

real(kind=8),dimension(:,:,:),allocatable::derivs  ! the derivatives of the potential 
                                             !   energy for all beads at once
real(kind=8)::k_B ! inverse temperature and boltzmann constant
real(kind=8)::afm_force   ! the returned force on the AFM atom
parameter(k_B=0.316679D-5)  ! the boltzmann constant in au/kelvin
real(kind=8)::dt_info ! time in ps for information
real(kind=8)::num_measure  ! how many times the temperature is measured..
integer::nbeads_store   ! the stored number of beads (RPMD)
character(len=20) keyword
character(len=80)::names,a80,coul_method
character(len=120) record
character(len=50)::fix_file ! file with numbers of fixed atoms
character(len=120) string
character(len=1)::qmdffnum
character(len=40)::commarg ! string for command line argument
integer::istat,readstat
character(len=10)::ensemble ! which thermodynamic ensemble is used
!     for MPI parallelization
integer::ierr,impi_error
integer::psize,source,status(MPI_STATUS_SIZE)
! For box simulations with periodic boundary conditions
real(kind=8)::e_pot_avg,e_kin_avg,e_tot_avg  ! averages for energies
real(kind=8)::e_cov_avg,e_noncov_avg 
real(kind=8),allocatable::q_result_act(:,:,:)
character(len=50)::baro_name,coul_name ! which barostat shall be used
!     periodic box size
real(kind=8)::xmin,ymin,zmin,xmax,ymax,zmax
!     for read in of local RPMD activation
character(len=1)::at_index(200)
character(len=20)::act_number
!   for umbrella samplings 
real(kind=8),allocatable::dxi_act(:,:)
real(kind=8)::xi_ideal,xi_real
integer::bias_mode,keylines_backup
integer::round,constrain
!   for sticking coefficient calculations
integer::num_traj  ! total number of trajectories
integer::count,dest,loop_large,loop_rest
integer::schedule,tag_mpi,traj_result
integer::numwork,traj_act
character(len=80)::sys_com
real(kind=8),allocatable::message(:)
real(kind=8)::coord_mat(3,3),coord_inv(3,3)
real(kind=8)::q_act_frac(3)
real(kind=8)::coll_ener
character(len=50)::ener_ref
logical::write_trajs
!   The OMP time measurement
real(kind=8)::time1_omp,time2_omp
!     the MPI rank 
integer::rank
!     the number of OMP threads 
integer::threads,id
real(kind=8)::duration_omp

!
!     Start MPI parallel computation:
!     Since the whole program is executed by all processes, all "serial"
!     parts will be executed only by processor #0
!
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,psize,ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
!
!     Activate the RPMD switch for correct application of gradient calculations
!
use_rpmd=.true.
!
!     The program calc_rate is not used..
!
use_calc_rate=.false.
!
!     The program stick_coeff is used! 
!
use_stick_coeff = .true.
!
!     The explore program is not used
!
use_explore = .false.
!
!     No frequency intensities as default
!
calc_freq_int = .false.

!
!
!     set up the structure and molecular mechanics calculation
!
call prog_initial(rank)
!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("stick_coeff")
      stop
   else
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'stick_coeff.x -help' or 'stick_coeff.x -h'."
      write(*,*)
   end if
end if

!
!     Read keywords from the key-file
!
call getkey(rank)
!
!      Option to change the atomic mass of one element in the system
!
change_mass=.false.
!do i = 1, nkey_lines
!   next = 1
!   record = keyline(i)
!   call gettext (record,keyword,next)
!   call upcase (keyword)
!   string = record(next:120)
!   if (keyword(1:20) .eq. 'CHANGE_MASS ') then
!       read(record,*) names,elem_mass,newmass
!       change_mass=.true.
!   end if
!end do

!if (change_mass) then
!   write(*,*) "The CHANGE_MASS option was activated!"
!   write(*,*) "All ", elem_mass," atoms will be assinged a mass of ",newmass, " a.u."
!end if



!
!     Read in the start structure for dynamics (keyword: XYZSTART)
!
fix_atoms=.false.
call getxyz(rank)
!
!     initialize the temperature, pressure and coupling baths
!
kelvin = 0.0d0
atmsph = 0.0d0
isothermal = .false.
isobaric = .false.
writestat=.true. ! print MD statistics!
allocate(centroid(3,natoms))

!-----------------------Dynamic-inputparameters--------------------------
!!
!!    All the following only for the rank 0 process if MPI us used   !!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Read the general dynamic input parameters   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   the default values:
!
nbeads=1  
ensemble="NVE"
nstep=0
dt=1.d0
dtdump=1.d0
bias_mode=0  ! no umbrella sampling as default
num_traj=0
coll_ener=0.d0
ener_ref="none"
write_trajs=.false.
allocate(dxi_act(3,natoms))
!
!    Read in the commands from the keyfile 
!
!    Both general MD commands and special sticking coefficient commands
!
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
!    Name of xyz input file (for information) 
   if (keyword(1:11) .eq. 'XYZSTART ') then
      read(record,*) names,xyzfile
!    Number of RPMD beads
   else if (keyword(1:20) .eq. 'RPMD_BEADS ') then
      read(record,*) names,nbeads
!    Which thermodynamic ensemble is used
   else if (keyword(1:13) .eq. 'ENSEMBLE ') then
      read(record,*) names,ensemble
!    How many MD steps are simulated 
   else if (keyword(1:11) .eq. 'STEPS ') then
      read(record,*) names,nstep
!    Time step of the Verlet algorithm
   else if (keyword(1:11) .eq. 'DELTAT ') then
      read(record,*) names,dt
!    Time between two xyz trajectory dumps
   else if (keyword(1:11) .eq. 'TDUMP ') then
      read(record,*) names,dtdump
!    Number of sticking coefficient trajectories
   else if (keyword(1:11) .eq. 'NUM_TRAJS ') then
      read(record,*) names,num_traj
!    The collision energy of the molecule in eV
   else if (keyword(1:11) .eq. 'COLL_ENER ') then
      read(record,*) names,coll_ener
!    If the collision energy shall be given with respect to the 
!     total velocity vector or with respect to the fraction of the 
!     vector normal to the surface
   else if (keyword(1:11) .eq. 'ENER_REF ') then
      read(record,*) names,ener_ref
!    If xyz files of all the trajectories shall be written into
!     the main calculation folder (else, only final structures are
!     written)
   else if (keyword(1:15) .eq. 'WRITE_TRAJS ') then
      write_trajs=.true.
   end if
end do
!
!     Check the the read in settings 
!
nvt=.false.
npt=.false.
if (nbeads .lt. 1) then
   if (rank .eq. 0) then
      write(*,*) "The number of RPMD beads must be positive!"
   end if
   call fatal
end if
call upcase(ensemble)
if (ensemble .eq. "NVT") then
   nvt=.true.
else if (ensemble .eq. "NPT") then
   npt=.true.
else if (ensemble .eq. "NVE") then
   nve=.true.
else 
   if (rank .eq. 0) then
      write(*,*) "No valid thermodynamic ensemble given! Choose either"
      write(*,*) " NVE, NVT or NPT."
   end if
   call fatal
end if
if (nstep .lt. 1) then
   if (rank .eq. 0) then
      write(*,*) "The number of MD timesteps must be greater than zero!"
   end if
   call fatal
end if
if (dt .lt. 1D-5) then
   if (rank .eq. 0) then
      write(*,*) "The MD timestep is too small to be useful!"
   end if
   call fatal
end if
if (num_traj .eq. 0) then
   if (rank .eq. 0) then
      write(*,*) "Please give the number of single sticking coefficient "
      write(*,*) " trajectories with the NUM_TRAJS command! (> 0)"
   end if
   call fatal 
end if
if (coll_ener .lt. 1D-9) then
   if (rank .eq. 0) then
      write(*,*) "Plase give the collision energy of the molecule with the "
      write(*,*) " COLL_ENER keyword in eV, larger than zero!"
   end if
   call fatal
end if
call upcase(ener_ref)
if (ener_ref .eq. "TOTAL") then
   if (rank .eq. 0) then
      write(*,*) "The collision energy is given for the total velocity."
   end if
else if (ener_ref .eq. "ORTHO") then
   if (rank .eq. 0) then
      write(*,*) "The collision energy is given orthogonal to the surface."
   end if
else
   if (rank .eq. 0) then
      write(*,*) "Please give the ENER_REF keyword!"
   end if
   call fatal
end if
if (write_trajs) then
   if (rank .eq. 0) then
      write(*,*) "xyz files of the single trajectories (trajxxx.xyz) will be written."
   end if
else
   if (rank .eq. 0) then
      write(*,*) "No trajectories will be written to xyz files."
   end if
end if

!
!     Read in the potential energy surface 
!
call read_pes(rank)

!
!    The write frequency
!
iwrite = nint(dtdump/(dt))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Read the parameter for the NVE ensemble     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Currently not needed
if (nve) then
   thermostat = 0
   kelvin=0
   thermo="none"
   ewald=.false.
   zahn=.true.
   coul_method="ZAHN"
   coul_cut=10.d0
   vdw_cut=10.d0
   andersen_step=70
   nose_q=100.d0
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'NVE {' .or. trim(adjustl(record(1:11))) &
                 &  .eq. 'NVE{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
!    For the periodic coulomb interactions
            if (keyword(1:13) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               if (rank .eq. 0) then
                  write(*,*) "The NVE section has no second delimiter! (})"
               end if
               call fatal
            end if
         end do
         exit
      end if
   end do
end if

!write(*,*) periodic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Read the parameter for the NVT ensemble     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if (nvt) then
!
!    The default values:
!
   kelvin=0
   thermo="none"
   andersen_step=70
   nose_q=100.d0
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'NVT {' .or. trim(adjustl(record(1:11))) &
                 &  .eq. 'NVT{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
!    The simulation temperature
            if (keyword(1:11) .eq. 'TEMP ') then
               read(record,*,iostat=readstat) names,kelvin
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "Correct format: TEMP [temperature (K)]"
                  end if
                  call fatal
               end if
!    The thermostat
            else if (keyword(1:20) .eq. 'THERMOSTAT ') then
               read(record,*,iostat=readstat) names,thermo
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "Correct format: THERMOSTAT [thermostat]"
                  end if
                  call fatal
               end if
               call upcase(thermo)
!    How often the Andersen thermostat velocity rescaling is applied
            else if (keyword(1:16) .eq. 'ANDERSEN_STEP ') then
               read(record,*,iostat=readstat) names,andersen_step
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "Correct format: ANDERSEN_STEP [No. of steps]"
                  end if
                  call fatal
               end if
!    The damping factor for the Nose-Hoover thermostat
            else if (keyword(1:13) .eq. 'NOSE_DAMP ') then
               read(record,*,iostat=readstat) names,nose_q
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "Correct format: NOSE_DAMP [integer (>1)]"
                  end if
                  call fatal
               end if
            end if
            if (keyword(1:13) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               if (rank .eq. 0) then
                  write(*,*) "The NVT section has no second delimiter! (})"
               end if
               call fatal
            end if
         end do
         exit
      end if
   end do
!
!     Check the the read in settings 
!
   if (kelvin .lt. 0.d0) then
      if (rank .eq. 0) then
         write(*,*) "Give a positive temperature!"
      end if
      call fatal
   end if
     
   if (trim(thermo) .eq. "ANDERSEN") then
      thermostat = 1
   else if (trim(thermo) .eq. "NOSE-HOOVER") then
      thermostat = 2
   else
      if (rank .eq. 0) then 
         write(*,*) "No availiable thermostat was chosen! Only the Andersen"
         write(*,*) "thermostat and the Nose-Hoov(e-Chain) thermostats are availiable at the moment." 
         write(*,*) "Add the keyword THERMOSTAT!"
      end if
      call fatal
   end if   
end if

if (mace_ase) then
   if (rank .gt. 0) then
      call system ("cp ../"//trim(xyzfile)//" .")
   end if
end if

!
!     define the beta parameter for thermodynamics
!
beta=1.d0/(kelvin*k_B)
!
!     convert timestep to atomic units!
!
dt_info=0.001*dt ! time in ps for information write out
dt = dt/2.41888428E-2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     Read more detailed/special MD parameters    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!      Evaluate time-dependent values of chosen coordinates
!
eval_coord=.false.
eval_step = 1
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:20) .eq. 'EVAL_COORD ') then
      read(record,*,iostat=readstat) names,eval_step
      if (readstat .ne. 0) then
         if (rank .eq. 0) then
            write(*,*) "Correct format: EVAL_COORD [write frequency (number)]"
         end if
         call fatal
      end if
      if (rank .eq. 0) then
         write(*,*) "The keyword EVAL_COORD was activated! All coordinates listed "
         write(*,*) " in the file 'coord_eval.inp' will be read in and their values "
         write(*,*) " will be written out every ",eval_step," steps to the file "
         write(*,*) " coord_eval.dat"
      end if
      call ev_coord_init  
      eval_coord=.true.
      exit
   end if
end do

!
!      Set random seed of Andersen thermostat to specific value
!
eval_coord=.false.
eval_step = 1
seed_manual = .false.
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:20) .eq. 'RANDOM_SEED ') then
      seed_manual = .true.
      read(record,*,iostat=readstat) names,seed_value
      if (readstat .ne. 0) then
         if (rank .eq. 0) then
            write(*,*) "Correct format: RANDOM_SEED [value (integer)]"
         end if
         call fatal
      end if
      if (rank .eq. 0) then
         write(*,*) "The random seed for the thermostat is set manually to ",seed_value
      end if
      exit
   end if
end do


!
!      Read in velocities from file as starting for momentum
!
read_vel=.false.
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:20) .eq. 'READ_VEL ') then
      read_vel=.true.
      if (rank .eq. 0) then
         write(*,*) "Starting velocities/momenta will be read in from file"
         write(*,*) " 'velocity_start.dat'!"
      end if
   end if
end do

if (read_vel) then
   allocate(vel_start(natoms,3))
   if (rank .eq. 0) then
      open(unit=291,file="velocity_start.dat",status="old")
      read(291,*) 
      read(291,*)
      do i=1,natoms
         read(291,*) idum,idum,vel_start(i,:)
      end do

      close(291)
   end if
end if

!
!     store the number of beads during structure gen.
!
nbeads_store=nbeads


!-----------------------Initialize-MD-Run--------------------------
!
!
!     If the system shall be described periodic of with hard box walls : 
!     Shift the system to the first octand
!     and check if the given box lengths are large enough to contain it     
!
if (periodic .or. box_walls) then
   if (.not. coord_vasp) then
      xmin=minval(x(:))
      ymin=minval(y(:))
      zmin=minval(z(:))

      x(:)=x(:)-xmin
      y(:)=y(:)-ymin
      z(:)=z(:)-zmin
   
      xmax=maxval(x(:))
      ymax=maxval(y(:))
      zmax=maxval(z(:))
      if (xmax .gt. boxlen_x) then
         if (rank .eq. 0) then
            write(*,*) "ERROR! You have chosen a periodic/hard box, but the elongation"
            write(*,*) " of the given structure in x-direction is too large!"
         end if
         call fatal
      else if (ymax .gt. boxlen_y) then
         if (rank .eq. 0) then
            write(*,*) "ERROR! You have chosen a periodic/hard box, but the elongation"
            write(*,*) " of the given structure in y-direction is too large!"
         end if
         call fatal
      else if (zmax .gt. boxlen_z) then
         if (rank .eq. 0) then
            write(*,*) "ERROR! You have chosen a periodic/hard box, but the elongation"
            write(*,*) " of the given structure in z-direction is too large!"
         end if
         call fatal
      end if
!
!     Place the initial structure(s) in the center of the box!
!     Deactivated since it tends to break calculations 
!  
!      x(:)=x(:)+0.5d0*(boxlen_x-xmax)
!      y(:)=y(:)+0.5d0*(boxlen_y-ymax)
!      z(:)=z(:)+0.5d0*(boxlen_z-zmax)
!
!     If a barostat shall be applied, set the initial values 
!
      vel_baro=0.d0
    end if
end if
if (npt .and. coord_vasp) then
   if (rank .eq. 0) then
      write(*,*) "Currently, the NpT ensemble cannot be used with VASP coordinates!"
   end if
   call fatal
end if
if (npt .and. .not. periodic) then
   if (rank .eq. 0) then
      write(*,*) "The NPT ensemble can only be used for periodic systems!"
   end if
   call fatal
end if
!
!     initialize and setup dynamics
!     --> allocate arrays for positions (q_i) and momenta (p_i) and 
!     convert to bohr 
!
allocate(q_i(3,natoms,nbeads_store),p_i(3,natoms,nbeads_store))
do k=1,nbeads
   do i=1,natoms
      q_i(1,i,k)=x(i)/bohr
      q_i(2,i,k)=y(i)/bohr
      q_i(3,i,k)=z(i)/bohr    
   end do
end do
!
!     Determine all fixed atoms for the dynamics
!
if (allocated(at_move)) deallocate(at_move)
allocate(at_move(natoms))
at_move=.true.
if (fix_atoms) then
   do i=1,natoms
      do j=1,fix_num
         if (fix_list(j) .eq. i) then
            at_move(i) = .false.
         end if
      end do
   end do
end if
!
!     For NPT ensemble: convert the pressure unit
!
if (npt) then
   pressure=pressure/prescon
end if
! 
allocate(derivs(3,natoms,nbeads))
!
!     For the AFM experiment: store the initial position of anchor atom
!
if (afm_run) then
   afm_move_first(:)=q_i(:,afm_move_at,1)
end if

!
!     If the internal coordinates shall be written to file, open it
!
if (int_coord_plot) open(unit=191,file="int_coord.out",status="unknown")
!
!-----------------------Print-Information-for-run------------------------
!
!
if (rank .eq. 0) then
   write(*,*) "----------------------------------------------------------"
   write(*,*) "The following molecular dynamics job will be started:"
   write(*,*) " - Input geometry: ",trim(xyzfile)
   if (nbeads .eq. 1) then
      write(*,*) " - Number of RPMD beads: 1 (classical)"
   else 
      write(*,'(a,i4)') "  - Number of RPMD beads: ",nbeads
   end if
   if (ensemble .eq. "NVT") then
      write(*,*) " - Thermodynamic ensemble: NVT "
   else if (ensemble .eq. "NPT") then
      write(*,*) " - Thermodynamic ensemble: NpT "
   else if (ensemble .eq. "NVE") then
      write(*,*) " - Thermodynamic ensemble: NVE "
   end if
   write(*,'(a,i10)') "  - Number of simulated MD steps:",nstep
   write(*,'(a,f12.5)') "  - MD time step (fs):",dt*2.41888428E-2
   write(*,'(a,f15.5)') "  - Total simulated timescale (ps):",dt*2.41888428E-2*nstep/1000.d0
   write(*,'(a,f13.5,a)') "  - Trajectory frames will be written every",dtdump," fs"
   if (nvt) then
      write(*,'(a,f12.5)') "  - Desired temperature (K): ",kelvin
      write(*,*) " - The temperature will be applied by the ",trim(thermo)," thermostat"
      if (periodic) then
      end if
      if (box_walls) then
         write(*,*) " - The hard walled box has the dimensions (x,y,z) (Angstron):"
         write(*,'(a,f14.6,f14.6,f14.6)') "    ",boxlen_x*bohr,boxlen_y*bohr,boxlen_z*bohr
      end if
      if (thermo .eq. "ANDERSEN") then
         write(*,'(a,i6,a)') "  - The Andersen velocity rescaling is done every ",andersen_step," MD steps"
      end if 
      if (thermo .eq. "NOSE-HOOVER") then
         write(*,'(a,f12.5,a)') "  - The Nose-Hoover thermostat damping Q is: ",nose_q," time steps"    
      end if 
   end if

   if (npt) then
      write(*,'(a,f12.5)') "  - Desired temperature (K): ",kelvin
      write(*,'(a,f12.5)') "  - Desired pressure (atm): ",pressure*prescon
      write(*,*) " - The temperature will be applied by the ",trim(thermo)," thermostat"
      write(*,*) " - The pressure will be applied by the ",trim(baro_name), " barostat"
      write(*,*) " - The periodic box has the initial dimensions (x,y,z) (Angstron):"
      write(*,'(a,f14.6,f14.6,f14.6)') "    ",boxlen_x*bohr,boxlen_y*bohr,boxlen_z*bohr
      if (box_walls) then
         write(*,*) " - The hard walled box has the dimensions (x,y,z) (Angstron):"
         write(*,'(a,f14.6,f14.6,f14.6)') "    ",boxlen_x*bohr,boxlen_y*bohr,boxlen_z*bohr
      end if
      if (thermo .eq. "ANDERSEN") then
         write(*,'(a,i6,a)') "  - The Andersen velocity rescaling is done every ",andersen_step," MD steps"
      end if
      if (thermo .eq. "NOSE-HOOVER") then
         write(*,'(a,f12.5,a)') "  - The Nose-Hoover thermostat damping Q is: ",nose_q," time steps"
      end if
   end if
   if (mirrors) then
      write(*,*) " - The following mirror planes will be present in the system:"
      do i=1,mirror_num
         if (mirror_dims(i) .eq. 1) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : x at",mirror_pos(i)*bohr," Ang."
         else if (mirror_dims(i) .eq. 2) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : y at",mirror_pos(i)*bohr," Ang."
         else if (mirror_dims(i) .eq. 3) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : z at",mirror_pos(i)*bohr," Ang."
         end if
      end do
   end if
end if

!
!     Initialize random seed for initial momenta
!
call random_init_local(0)
!
!     call the initialization routine
!
!     Broadcast positions and momenta to all threads
!
call mpi_bcast(q_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(p_i,natoms*3*nbeads,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_barrier(mpi_comm_world,ierr)
call mpi_barrier(mpi_comm_world,ierr)
!
!    Morse for "here is Caracal"
if (rank .eq. 0) then
   write(*,*) 
   write(*,*)  " .... . .-. .  .. ...  -.-. .- .-. .- -.-. .- .-.."
end if
!
!    Measure the needed time for dynamic steps
!
if (rank .eq. 0) then
   call cpu_time(time1)
   time1_omp = omp_get_wtime()
end if
call mpi_barrier(mpi_comm_world,ierr)
!
!     integrate equations of motion to take a time step
!
if (afm_run) afm_steps=nstep
l=1
j=1
if (afm_run) force_avg(:)=0
if (afm_run) len_avg(:)=0

allocate(message(1+nbeads*natoms*3))


!
!    Now enter the main loop for the trajectories
!
!    Distribute the jobs on all MPI slave processes
!
loop_large=int((num_traj)/(psize-1))
loop_rest=mod(num_traj,psize-1)
!     define default count value
count=1
tag_mpi=0
allocate(q_result_act(3,natoms,psize-1))
open(unit=57,file="XDATCAR_final",status="replace")
xdat_first=.true.
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
               cycle
            end if
         end if
         dest=j
!
!     Send number of trajectory to slave process
!
         call mpi_send(schedule, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
      end do
!
!     recieve results from all slaves for the i'th round
!     Only consider first bead now (layer maybe centroid?)
!
      do j=1,psize-1
         if (i .lt. loop_large+1) then
            call mpi_recv(message, 3*natoms*nbeads+1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,tag_mpi, &
               & MPI_COMM_WORLD,status,ierr)
            traj_result=int(message(1))
            do k=1,1
               do l=1,natoms
                  do m=1,3
                     q_result_act(m,l,j)=message((l-1)*3+m+1)
                  end do
               end do
           end do
         end if
      end do
!
!     Append final structures of current run to XDATCAR_end file
!
      do j=1,psize-1
         if (i .lt. loop_large+1) then
            traj_act=(schedule-1)*(psize-1)+j
            if (xdat_first) then
               write(57,'(a,i9,a)') "Final sticking coefficient structures, step=trajectory number"
               write(57,*) vasp_scale
               write(57,*) vasp_a_vec
               write(57,*) vasp_b_vec
               write(57,*) vasp_c_vec
               do k=1,nelems_vasp
                  write(57,'(a,a)',advance="no") " ",trim(vasp_names(k))
               end do
               write(57,*)
               write(57,*) vasp_numbers(1:nelems_vasp)
               xdat_first=.false.
            end if
            write(57,'(a,i9)') "Direct final structure of trajectory ",traj_act
!
!     As in usual CONTCAR files, give the positions in direct coordinates!
!     convert them back from cartesians, by using the inverse matrix
!
            coord_mat(:,1)=vasp_a_vec(:)
            coord_mat(:,2)=vasp_b_vec(:)
            coord_mat(:,3)=vasp_c_vec(:)
         
            call matinv3(coord_mat,coord_inv)

            do k=1,natoms
               q_act_frac=matmul(coord_inv,q_result_act(:,k,j)*bohr)
               write(57,*) q_act_frac(:)
            end do
         end if
      end do
   end do
!
!     If the total number of jobs is no multiple of the number of processors,
!      do the remaining jobs separately
!
   do j=1,loop_rest
      dest=j
!
!     Send number of trajectory to slave process
!
      call mpi_send(i-1, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
   end do

!
   do j=1,loop_rest
      call mpi_recv(message, 3*natoms*nbeads+1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE,tag_mpi, &
         & MPI_COMM_WORLD,status,ierr)
      traj_result=int(message(1))
      do k=1,1
         do l=1,natoms
            do m=1,3
               q_result_act(m,l,j)=message((l-1)*3+m+1)
            end do
         end do
      end do
   end do
!
!     Append final structures to XDATCAR_end file
!
   do j=1,loop_rest
      traj_act=(i-2)*(psize-1)+j
      if (xdat_first) then
         write(57,'(a,i9,a)') "Final sticking coefficient structures, from stick_coeff.x"
         write(57,*) vasp_scale
         write(57,*) vasp_a_vec
         write(57,*) vasp_b_vec
         write(57,*) vasp_c_vec
         do k=1,nelems_vasp
            write(57,'(a,a)',advance="no") " ",trim(vasp_names(k))
         end do
         write(57,*)
         write(57,*) vasp_numbers(1:nelems_vasp)
         xdat_first=.false.
      end if
      write(57,'(a,i9)') "Direct final structure of trajectory ",traj_act
!
!     As in usual CONTCAR files, give the positions in direct coordinates!
!     convert them back from cartesians, by using the inverse matrix
!
      coord_mat(:,1)=vasp_a_vec(:)
      coord_mat(:,2)=vasp_b_vec(:)
      coord_mat(:,3)=vasp_c_vec(:)

      call matinv3(coord_mat,coord_inv)

      do k=1,natoms
         q_act_frac=matmul(coord_inv,q_result_act(:,k,j)*bohr)
         write(57,*) q_act_frac(:)
      end do
   end do

   do j=1,loop_rest
      dest=j
      call mpi_send(-1, count, MPI_DOUBLE_PRECISION, dest,tag_mpi,MPI_COMM_WORLD,ierr)
   end do

else 
   source=0
   message=rank
   do
      call mpi_recv(numwork, count, MPI_DOUBLE_PRECISION, 0,tag_mpi,MPI_COMM_WORLD,status,ierr)
!
!     If no more samplings are to do, exit with the current worker
!       
      if (numwork .eq. -1) then
         exit
      end if
      traj_act=(numwork-1)*(psize-1)+rank
      write(*,'(a,i7,a,i7,a)') " Start trajectory",(numwork-1)*(psize-1)+rank," of ",num_traj," ..."
!
!     File name for trajectory file of each trajectory (only is write_traj is activated!)
!
      if (traj_act .lt. 10) then
         write(sys_com,'(a,i1,a)') "../traj",traj_act,".xyz"
      else if (traj_act .lt. 100) then
         write(sys_com,'(a,i2,a)') "../traj",traj_act,".xyz"
      else if (traj_act .lt. 1000) then
         write(sys_com,'(a,i3,a)') "../traj",traj_act,".xyz"
      else if (traj_act .lt. 10000) then
         write(sys_com,'(a,i4,a)') "../traj",traj_act,".xyz"
      else
         write(sys_com,'(a,i4,a)') "../traj",traj_act,".xyz"
      end if
      if (write_trajs) then
         open(unit=28,file=sys_com,status="replace")
      end if
!
!     Initialize momenta of atoms
!
      call mdinit(derivs,xi_ideal,dxi_act,bias_mode,rank)
!
!     Ensure that all initial momenta of the two last atoms of the structure 
!     (the gas phase molecule) are directed to the surface:
!     If the z-momentum is positive, add a minus sign
!
      do j=natoms-1,natoms
         if (p_i(3,j,1) .gt. 0) then
            p_i(3,j,1)=-p_i(3,j,1)
         end if
      end do

      do istep = 1, nstep
         analyze=.false.
!
!     Write the trajectory frame every iwrite timesteps
!
        
         if (mod(istep,iwrite) .eq. 0 .and. write_trajs) then
            analyze=.true.
         end if
         round=0
         constrain=-1
         epot=0.d0
         call verlet (istep,dt,derivs,epot,ekin,afm_force,xi_ideal,xi_real,dxi_act, &
                  & round,constrain,analyze,rank)
      end do
      if (write_trajs) then
         close(28)
      end if

!
!     Send back the message with the result: does the molecule stick to the surface?
!
      message(1)=1    
      do i=1,nbeads
         do j=1,natoms
            do k=1,3
               message((i-1)*natoms*3+(j-1)*3+k+1)=q_i(k,j,i)
            end do
         end do
      end do
      call mpi_send(message, natoms*3*nbeads+1, MPI_DOUBLE_PRECISION, source,tag_mpi,MPI_COMM_WORLD,ierr)
   end do
end if
call mpi_barrier(mpi_comm_world,ierr)
!
!    If MACE is called externally, terminate the python ASE script
!
if (mace_ase) then
   call system("touch end_cycle")
end if


!
!    If MACE is called externally, terminate the python ASE script
!
if (mace_ase) then
   if (rank .gt. 0) then
      if (rank .lt. 10) then
         write(sys_com,'(a,i1)') "rank",rank
      else if (rank .lt. 100) then
         write(sys_com,'(a,i2)') "rank",rank
      else
         write(sys_com,'(a,i3)') "rank",rank
      end if
      call system("touch end_cycle")
      call sleep(1)
      call system("cd ..")
      call chdir ("..")
      call system("rm -r "//trim(sys_com))
   end if

   call system("touch end_cycle")
end if
!
!    calculate the needed time for dynamic steps
!
if (rank .eq. 0) then
   call cpu_time(time2)
   time2_omp = omp_get_wtime()

   duration=time2-time1
   duration_omp=time2_omp-time1_omp
   write(*,*) 
   write(*,*) "Final structures of all trajectories are written to "
   write(*,*) " the file XDATCAR_final."
   write(*,*)
   write(*,'(A, F12.3, A)') " The calculation needed a time of",duration_omp," seconds."
!write(*,'(A, i10, A)') " The calculation needed a time of ",duration_omp," seconds."
   write(*,*)
   write(*,*) ".. .----. -- -.. --- -. . "
   write(*,*) "Sticking coefficient calculation successfully finished!"
end if
call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)


end program stick_coeff

