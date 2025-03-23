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
!     #################################################################
!     ##                                                             ##
!     ##  program explore  --  pseudo QM package                     ##
!     ##                                                             ##
!     #################################################################
!
!
!     "explore" manages different types of calculations for one 
!     input structure. Availiable are:
!     - energy (always)
!     - gradient 
!     - hessian and normal mode analysis 
!     - optimization of minima 
!     - optimization of transition states
!     - optimization of reaction paths (IRC)
!     

program explore
use general
use evb_mod
use pbc_mod
use debug
use omp_lib

implicit none
integer::num_arg,input_unit,i,qmdff_energies_unit,asciinum
integer::qmdffnumber,nat,nat3,l,m
integer::readstat  ! for error handling
integer::state_open  ! if a file was opened successfully
integer::reference_counter,parameters  !zahl der Strukturen im xyz-file
real(kind=8)::energy,e_qmdff1,e_qmdff2,e_evb
real(kind=8),dimension(:,:),allocatable::coord,energies_qmdff 
real(kind=8),dimension(:,:),allocatable::energies_tmp  ! für temporäres zwischenspeichern
real(kind=8),dimension(:,:),allocatable::g_evb
real(kind=8),dimension(:),allocatable::int_coord,geo_int,geo_xyz1  ! for internal coordinates
real(kind=8),dimension(:,:),allocatable::xyz2,geo_xyz
real(kind=8),dimension(:,:),allocatable::eigvecs   ! matrix with eigenvector amplitudes
integer,dimension(:),allocatable::atind   ! array with element numbers
character(len=70)::fffile1,fffile2,fffile3,fileinfo,filegeo,fileenergy
character(len=70)::filets,filets2,xyzfile
character(len=60)::test,names,soschl_pre1,soschl_hess
character(len=1)::qmdffnum
character(len=20) keyword
character(len=120) record
character(len=120) string
character(len=120) coupling
character(len=80)::jobtype
character(len=80) date
! for former egrad.x capabilities (energies+gradient+QMDFF debug)
!  for Wilson matrix printout
real(kind=8),allocatable::internal(:),B_mat(:,:),dB_mat(:,:,:)
logical::path_struc,print_wilson
integer::wilson_mode,ref_count
integer::time_int1,time_int2
! for manual internal coordinates
integer,dimension(:,:),allocatable::coord_tmp
integer,dimension(:),allocatable::coord_types
character(len=40)::coord_line
! Time measurement
real(kind=8)::time1_omp,time2_omp
!
integer mode,next,j,k,readstatus,dg_evb_mode,mat_size
integer::int_mode  ! method for defining internal coordinates (if used)
logical::exist,exists,has_next,coupl1,par_soschl

logical::grad,frequency,opt_min,orca_fake,ts_opt,calc_irc,calc_egrad
! for EVB-QMDFF hessian calculation
real(kind=8),dimension(:,:),allocatable::hess,mass_mat
real(kind=8),dimension(:),allocatable::h_out,freqs
integer::maxline
character(len=40)::commarg ! string for command line argument
! for fragment calculations
integer::indi_tmp(100)
character(len=80) a80,adum
integer::maxcycle,rest
!     the MPI rank (here always 0)
integer::rank
integer::nkey_test

real(kind=8),allocatable::xyz_init(:,:)  ! geometry for pGFN-FF init
logical::read_init_struc
character(len=2),allocatable::names_init(:)   ! element symbols for pGFN-FF init
real(kind=8),allocatable::grad_act(:,:)
real(kind=8)::energy_act
real(kind=8),allocatable::rab_ref(:,:)  ! for bond determinations in PES TOPOL mode
real(kind=8)::r_act,dist,r_ref
integer::ind_act,bond_sum,ang_sum
integer::nel,nbf ! for EHT basis set
real(kind=8),allocatable::at_charges(:) ! partial charges of all atoms
real(kind=8),allocatable::wbo(:,:)  ! Wiberg-Mayer bond orders
logical::okbas
real(kind=8)::eel
real(kind=8)::bond_order(100000)
real(kind=8)::wbo_cutoff
integer,allocatable::el_inds(:)  ! element indices 

!
!     Set MPI rank to zero for this program
!
rank=0
!
!     no MPI is used
!
use_mpi=.false.
!
!     The program calc_rate is not used..
!
use_calc_rate=.false.
!
!     The explore program is used
!
use_explore = .true.
!
!     No frequency intensities as default
!
calc_freq_int = .false.
!
!     set up the structure and mechanics calculation
!
call prog_initial(rank)

!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("explore")
      stop
   else
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'explore.x -help' or 'explore.x -h'."
   end if
end if
!
!     Open logfile for EXPLORE calculation on single structure
!     print header and time
!

open (unit=15,file="explore.log",status="unknown")
write(15,*) "--- EXPLORE CALCULATION  ---"
write(15,*)
write(15,*) "Calculation initiated at: "
call timestamp ( )

!
!     Read keywords from the key-file
!
call getkey(rank)
!
!     Read in the (start) structure for calculations
!
nkey_test=nkey_lines
call getxyz(rank)
!
!     Read in the potential energy surface parameters
!
call read_pes(rank)

!
!     Read in the relevant keywords for explore.x
!
grad=.true. ! calculate always gradients...
geomax=100
calc_frag=.false.
calc_egrad=.false.
frequency=.false.  
opt_min=.false.
ts_opt=.false.
calc_irc=.false.
! 
!    the default values
!
maxiter=500
stepmax=0.05d0
ethr=1D-6
gthr=5D-5
dthr=5D-5
gmaxthr=2D-4
dmaxthr=2D-4
newton_raphson=.false.
geoopt_algo = "cg"
irc_maxstep=200
irc_steplen=0.1d0
irc_eulerlen=0.005d0
irc_ethr=1D-6
irc_gthr=1D-6

!
!     Which job shall be done?
!
nkey_lines=nkey_test
jobtype="NONE"
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'JOB ') then
      read(record,*,iostat=readstat) names,jobtype
      if (readstat .ne. 0) then
         write(*,*) "Correct format: JOB [job-type]"
         call fatal
      end if
      call upcase(jobtype)
   else if (keyword(1:11) .eq. 'XYZSTART ') then
      call getword (record,xyzfile,next)
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      inquire(file=xyzfile,exist=exist)
   end if
end do

if (trim(jobtype) .eq. "OPT_MIN") then
   write(*,*) "OPT_MIN: A local geometry optimization will be done."
else if (trim(jobtype) .eq. "OPT_TS") then
   write(*,*) "OPT_TS: A transition state optimization will be done."
   ts_opt=.true.
else if (trim(jobtype) .eq. "EGRAD") then
   write(*,*) "EGRAD: Energies and gradients for a number of structures"
   write(*,*) " in the XYZSTART file will be calculated."
   calc_egrad=.true.
else if (trim(jobtype) .eq. "IRC") then
   write(*,*) "IRC: A reaction path will be optimized."
   calc_irc=.true.
else if (trim(jobtype) .eq. "FREQ") then
   write(*,*) "FREQ:  The Frequencies for the current structure will & 
                 & be calculated."
else if (trim(jobtype) .eq. "OPTFREQ") then
   write(*,*) "OPTFREQ: A local minimum optimization followed by a "
   write(*,*) "    frequency calculation will be done."
else 
   write(*,*) "No valid job was chosen! Choose one of the following: "
   write(*,*) " OPT_MIN, OPT_TS, IRC, FREQ, OPTFREQ"
   call fatal
end if

if (pes_topol) then
   if (.not. calc_egrad) then
      write(*,*) "With PES TOPOL, only the EGRAD option can be chosen!"
      call fatal
   end if
end if
!
!    Read in sections for different job types 
!
!    A: The OPT_MIN and OPT_TS jobs
!

if (trim(jobtype) .eq. "OPT_MIN" .or. trim(jobtype) .eq. "OPT_TS" &
            &  .or. trim(jobtype) .eq. "OPTFREQ") then
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'OPT {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'OPT{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)

            if (keyword(1:20) .eq. 'MAXITER ') then
               read(record,*,iostat=readstat) names,maxiter
               if (readstat .ne. 0) then
                  write(*,*) "Correct format MAXITER [No. of iterations]"
                  call fatal
               end if
!     which optimization algorithm shall be used 
            else if (keyword(1:15) .eq. 'ALGO ') then
               read(record,*,iostat=readstat) names,adum
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ALGO [specifier]"
                  call fatal
               end if
               call upcase(adum)
               if (trim(adum) .eq. "CONJ_GRAD") then
                  geoopt_algo = "cg"
               else if (trim(adum) .eq. "BFGS") then
                  geoopt_algo = "bfgs"
               end if

!     maximum length of a single optimization step
            else if (keyword(1:20) .eq. 'STEPSIZE ') then
               read(record,*,iostat=readstat) names,stepmax
               if (readstat .ne. 0) then
                  write(*,*) "Correct format STEPSIZE [Opt. stepsize (A)]"
                  call fatal
               end if
!     energy change convergence criterion
            else if (keyword(1:20) .eq. 'CRIT_DE ') then
               read(record,*) names,ethr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ETHR [energy (Hartree)]"
                  call fatal
               end if
!     gradient norm convergence criterion
            else if (keyword(1:20) .eq. 'CRIT_GNORM ') then
               read(record,*) names,gthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ETHR [gradient comp.]"
                  call fatal
               end if
!     geometry step norm convergence criterion
            else if (keyword(1:20) .eq. 'CRIT_QNORM ') then
               read(record,*) names,dthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format DTHR [geom. step]"
                  call fatal
               end if
!     largest component in gradient vector
            else if (keyword(1:20) .eq. 'CRIT_GMAX ') then
               read(record,*) names,gmaxthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format GMAXTHR [gradient comp.]"
                  call fatal
               end if
!     largest component in geometry change vector
            else if (keyword(1:20) .eq. 'CRIT_QMAX ') then
               read(record,*) names,dmaxthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format GMAXTHR [geom. change]"
                  call fatal
               end if
!     If newton raphson shall be used (else: the better P-RFO method)
!      ---> Currently: always Newton-Raphson!
!            else if (keyword(1:20) .eq. 'NEWTON_RAPHSON ') then
!               if (.not. newton_raphson) then
!                  write(*,*) "The simple Newton-Raphson method for PES searching &
!                               &will be used."
!                  write(15,*) "The simple Newton-Raphson method for PES searching &
!                               &will be used."
!               end if
!               newton_raphson=.true.
!     If internal coordinates shall be defined manually
            else if (keyword(1:16) .eq. 'READ_COORD ') then
               read_coord=.true.
            end if
         end do
      end if
   end do
end if
newton_raphson=.true.
!
!    B: The IRC jobs
!

if (trim(jobtype) .eq. "IRC") then
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'IRC {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'IRC{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
!     Desired number of maximum IRC steps in each direction
            if (keyword(1:20) .eq. 'MAXSTEP ') then
               read(record,*) names,irc_maxstep
               if (readstat .ne. 0) then
                  write(*,*) "Correct format MAXSTEP [step number]"
                  call fatal
               end if
!     Desired length of an IRC step ((amu)^1/2*bohr)
            else if (keyword(1:20) .eq. 'STEPLEN ') then
               read(record,*) names,irc_steplen
               if (readstat .ne. 0) then
                  write(*,*) "Correct format STEPLEN [Opt. stepsize (A)]"
                  call fatal
               end if
!     Desired length of an elementary Euler step ((amu)^1/2*bohr)
!     --> should be much smaller than the IRC step itself!
            else if (keyword(1:20) .eq. 'EULERLEN ') then
               read(record,*) names,irc_eulerlen
               if (readstat .ne. 0) then
                  write(*,*) "Correct format EULERLEN [size (A)]"
                  call fatal
               end if
!     IRC energy change convergence criterion
            else if (keyword(1:20) .eq. 'ETHR ') then
               read(record,*) names,irc_ethr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ETHR [energy (Hartree)]"
                  call fatal
               end if
!     IRC gradient norm convergence criterion
            else if (keyword(1:20) .eq. 'GTHR ') then
               read(record,*) names,irc_gthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format GTHR [gradient (Ht/A]"
                  call fatal
               end if
!     If internal coordinates shall be defined manually
            else if (keyword(1:16) .eq. 'READ_COORD ') then
               read_coord=.true.
            end if

         end do
      end if
   end do
end if

!
!    C: The frequency calculation jobs 
!

calc_frag=.false.
orca_fake=.false.
if (trim(jobtype) .eq. "FREQ" .or. trim(jobtype) .eq. "OPTFREQ") then
   frequency=.true.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'FREQ {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'FREQ{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
!    If only parts of the system (a fragment) shall be calculated
            if (keyword(1:20) .eq. 'CALC_FRAG ') then
               calc_frag=.true.
            else if (keyword(1:20) .eq. 'ORCA_FAKE ') then
               orca_fake=.true.
            end if
         end do
      end if
   end do
end if

!
!  If a fragment is used, read in the fragment atoms (indi-array)
!


nats=natoms
nat=natoms

if (calc_frag) then
   inquire(file="fragment.inp",exist=exist)
   if (.not. exist) then
      write(*,*) "The file fragment.inp could not be found!"
      write(*,*) "This file is needed for EVB-QMDFF fragment calculations!"
      call fatal
   end if
   i=0
   open(unit=69,file="fragment.inp",status="old")
      do 
         i=i+1
         read(69,*,iostat=readstat) indi_tmp(i)
         if (readstat .ne. 0) then   
            i=i-1         
            exit
         end if
      end do
   close(69)
   nats=i
   allocate(indi(nats)) 
   allocate(indi3(nats*3))
!
!  fill the final array and control if the indices are legal
!
   do j=1,nats
      indi(j)=indi_tmp(j)
      if ((indi(j) .gt. nat) .or. (indi(j) .le. 0)) then
         write(*,*) "ERROR: one index stated in fragment.inp isn´t part of"
         write(*,*) "the molecule!"
         write(*,*) "line i::",j,", indi(i):",indi(j)
         call fatal
      end if
   end do
!
!  define array of a shape (3*N)
!
   do i=1,nats
      do j=1,3
         indi3((i-1)*3+j)=(indi(i)-1)*3+j
      enddo
   enddo
   write(15,*)
   write(15,*) "A CALC_FRAG calculation was requested!"
   write(15,*) "Only a part of the full structure was calculated"
   write(15,*) "The picked atoms are:"
   write(15,*)
   do i=1,natoms
      do j=1,nats
         if (indi(j) .eq. i) then
            write(15,*) i,name(i),"--->",j
         end if
      end do
   end do
   write(15,*)
else 
   allocate(indi(nats))
   do i=1,nats
     indi(i)=i
   end do
end if 


call system_clock(time_int1)
time1_omp = omp_get_wtime()

!
!    For energy+gradient calculations: start here!
!
if (calc_egrad) then
!
!     If the PES TOPOL mode is activated, initialize the reference bond 
!      length array (r0ab)
!
   if (pes_topol) then
      if (topol_bonds .eq. "VDW") then
         allocate(rab_ref(94,94))
         call setr0 (94,rab_ref)
      end if
!
!    For EHT calculations: set zeta basis functions and valence electrons
!
      if (topol_bonds .eq. "EHT") then
         call setvalel
         call setZETAandIP
         wbo_cutoff=0.5d0*topol_eht_cutoff
      end if
      allocate(coord_tmp(1000000,5))
   end if
!
!     If the debugging mode shall be started to print out more details!
!
   do_debug=.false.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:16) .eq. 'DEBUG ') then
         do_debug=.true.
      end if
   end do
!
!     If, for coordinate analysis etc., the Wilson matrix for each structure shall
!     be written to file 
!
   print_wilson=.false.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:16) .eq. 'PRINT_WILSON ') then
         read(record,*,iostat=readstat) names,wilson_mode
         if (readstat .ne. 0) then
            write(*,*) "The keyword PRINT_WILSON is incomplete!"
            call fatal
         end if
         print_wilson=.true.
         if ((wilson_mode .ne. 1) .and. (wilson_mode .ne. 2)) then
            write(*,*) "ERROR! No valid number for PRINT_WILSON given!"
            call fatal
         end if
      end if
   end do


!
!     Print debug information message
!     Also allocate debug arrays
!
   if (do_debug) then
      write(*,*) "THE DEBUGGING MODE WAS REQUESTED; ADDITIONAL INFOS WILL & 
          & BE PRINTED (if QMDFFs are used)"

      open(unit=288,file="bond_debug.dat",status="replace")
      open(unit=289,file="ang_debug.dat",status="replace")
      open(unit=290,file="dihed_debug.dat",status="replace")
      open(unit=291,file="dispersion_debug.dat",status="replace")
      open(unit=292,file="coulomb_debug.dat",status="replace")
!
!     QMDFF energy components: (no QMDFFs, no structures, no energyterms)
!
      allocate(qmdff_parts(2,1000,10000))
!
!     the labels for the different QMDFF energy components 
!
      allocate(parts_labels(2,10000))
      parts_labels=" "
   end if
!
!     If the PES TOPOL option is used, open a file for internal coordinate values 
!
   if (pes_topol) then
      open(unit=214,file="pes_topol.dat",status="replace")
      write(214,'(a)') "#This file contains the definition and actual values of "
      write(214,'(a)') "# all internal coordinates of each given trajectory frame."
      write(214,'(a)') "#Bond lenghts are given in bohr, angles in radians"
      if (topol_bonds .eq. "VDW") then
         write(214,'(a)') "# No.     type         atom1    atom2    atom3    atom4    value"
      else if (topol_bonds .eq. "EHT") then
         write(214,'(a)') "# No.     type         atom1    atom2    atom3    atom4  &
                  &   value      (bond-order)"
      end if
      write(*,*) "The definitions and values of all internal coordinates are written &
              &to 'pes_topol.dat'"
   end if  

!
!     Read in the internal coordinates from file if the Wilson matrix 
!     shall be calculated for each structure
!     Perform no allocation for the PES TOPOL mode
!
   if (print_wilson) then
      open(unit=215,file="wilson_mat.dat",status="replace")
      write(215,'(a)') "# This file contains the Wilson matrix values of each structure"
      write(215,'(a)') "# Values given in bohrs or radians"
      write(215,'(a)') "# No. int. coord.(i)   No. cart. coord(j)        B(i,j)"
      if (wilson_mode .eq. 2) then
         open(unit=216,file="wilson_deriv.dat",status="replace")
         write(216,'(a)') "# This file contains the Wilson matrix derivatives of each structure"
         write(216,'(a)') "# Values given in bohrs or radians"
         write(216,'(a)') "# No. int. coord.(i)   No. cart. coord(j)    No.cart.coord(k)     B'(i,j,k)"
      end if
      write(*,*) "The Wilson matrices of the structures are written to 'wilson_mat.dat'"
      if (wilson_mode .eq. 2) then
         write(*,*) "The Wilson matrix derivatives are written to 'wilson_derivs.dat'"
      end if
   end if

!
!     produce the EVB-QMDFF-energies for all points
!
   ref_count=0
   open(unit=166,file=xyzfile,status='old',iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",xyzfile," with the geometries to be calculated is missing!"
      write(*,*) "Please edit the XYZSTART keyword!"
      call fatal
   end if
   if (.not. pes_topol) then
      open(unit=44,file="energies.dat",status='unknown')
      open(unit=45,file="gradients.dat",status='unknown')
      write(*,*)  " -.-. .- .-. .- -.-. .- .-.. "
      write(*,*) "Calculating energy and gradient ..."
   else
      write(*,*)  " -.-. .- .-. .- -.-. .- .-.. "
      write(*,*) "Perform coordinate analyses  ..."
   end if
   allocate(g_evb(3,natoms))
   allocate(coord(3,natoms))
!
!     open file with single QMDFF energies
!
   if (.not. pes_topol) then
      open (unit=99,file="single_qmdff.dat",status='unknown')
   end if
   if (treq) open (unit=48,file="treq.out",status="unknown")
   if (int_grad_plot) open (unit=192,file="int_grad.out",status="unknown")
   do
      ref_count=ref_count+1
      call next_geo(coord,natoms,166,has_next)
      if (.not.has_next) exit
!
!     If the PES TOPOL mode is activated, determine the internal coordinates 
!
      if (pes_topol) then
!
!     First, determine all element indices
!
         if (allocated(el_inds)) deallocate(el_inds)
         if (allocated(coord_def)) deallocate(coord_def)
         
         allocate(el_inds(natoms))
         do i=1,natoms
            call elem(el_names(i),el_inds(i))
         end do
!
!     Now, determine the list of bonds in the system
!         
!     A: comparison of atom distances with VDW pair radii
!
         ind_act=0
         if (topol_bonds .eq. "VDW") then
            do i=1,natoms
               do j=i+1,natoms
                  r_act=dist(i,j,coord)
!
!     Optional: scale the VDW reference radii with the topol_vdw_scale parameter
!      the larger the parameter, the larger the reference bonds get, and the 
!      more bonds in the system will be detected
!
                  r_ref=rab_ref(el_inds(i),el_inds(j))*bohr/2*1.2*topol_vdw_scale
                  if (r_act .lt. r_ref) then
                     ind_act=ind_act+1
                     coord_tmp(ind_act,1) = 1
                     coord_tmp(ind_act,2) = i
                     coord_tmp(ind_act,3) = j
                  end if
               end do
            end do
         end if
!
!     B: Perform an Extended Hückel Theory (EHT) calculation and get
!       Wiberg-Mayer bond orders. Define all bonds which bond orders are 
!       above value topol_eht_cutoff
!
         if (topol_bonds .eq. "EHT") then
            call basis0_eht(natoms,el_inds,nel,nbf)
            call basis_eht (natoms,el_inds,nbf,okbas) 

            if (.not. okbas) then
               write(*,*) "Problem with EHT basis set! Maybe try BONDS VDW?"
               call fatal
            end if
            if (allocated(at_charges)) deallocate(at_charges)
            if (allocated(wbo)) deallocate(wbo)
            allocate(at_charges(natoms))
            allocate(wbo(natoms,natoms))
            at_charges=0.d0 ! set partial charges of atoms to zero (seems to work well)
            call ehtfull(natoms,el_inds,coord/bohr,at_charges,nel,nbf,eel,wbo,1.0d0)
!
!     Check all bond orders: if one is above the cutoff, set the bond
!
            do i=1,natoms
               do j=i+1,natoms
                  if (wbo(i,j) .gt. wbo_cutoff) then
                     ind_act=ind_act+1
                     coord_tmp(ind_act,1) = 1
                     coord_tmp(ind_act,2) = i
                     coord_tmp(ind_act,3) = j
                     bond_order(ind_act) = wbo(i,j)
                  end if
               end do
            end do
         end if
         bond_sum=ind_act

!
!     Determine the bond angles: all pairs of neighbored bonds!
!
         do i=1,bond_sum
            do j=i+1,bond_sum
!
!     i2=j2: i3-i2-j3
!
               if (coord_tmp(i,2) .eq. coord_tmp(j,2)) then
                  ind_act = ind_act +1   
                  coord_tmp(ind_act,1) = 2
                  coord_tmp(ind_act,2) = coord_tmp(i,3)
                  coord_tmp(ind_act,3) = coord_tmp(i,2)
                  coord_tmp(ind_act,4) = coord_tmp(j,3)
!
!     i2=j3: i3-i2-j2
!
               else if (coord_tmp(i,2) .eq. coord_tmp(j,3)) then
                  ind_act = ind_act +1
                  coord_tmp(ind_act,1) = 2
                  coord_tmp(ind_act,2) = coord_tmp(i,3)
                  coord_tmp(ind_act,3) = coord_tmp(i,2)
                  coord_tmp(ind_act,4) = coord_tmp(j,2)
!
!     j2=i3: i2-i3-j3
!
               else if (coord_tmp(i,3) .eq. coord_tmp(j,2)) then
                  ind_act = ind_act +1
                  coord_tmp(ind_act,1) = 2
                  coord_tmp(ind_act,2) = coord_tmp(i,2)
                  coord_tmp(ind_act,3) = coord_tmp(i,3)
                  coord_tmp(ind_act,4) = coord_tmp(j,3)
!
!     j3=i3: i2-i3-j2
!
               else if (coord_tmp(i,3) .eq. coord_tmp(j,3)) then
                  ind_act = ind_act +1
                  coord_tmp(ind_act,1) = 2
                  coord_tmp(ind_act,2) = coord_tmp(i,2)
                  coord_tmp(ind_act,3) = coord_tmp(i,3)
                  coord_tmp(ind_act,4) = coord_tmp(j,2)
               end if
            end do
         end do

         ang_sum=ind_act-bond_sum
!
!     Determine the dihedral or out of plane angles: all pairs of angles
!      with two identical atoms!
!     dihedral: third atoms must be bound at two different center atoms
!     out-of-plane: third atoms must be bound at the same center atom
!
         do i=1,ang_sum
            do j=i+1,ang_sum
!       
!     i2=j2 and i3=j3 (rest: i4, j4): out-of-plane
!
               if ((coord_tmp(i+bond_sum,2) .eq. coord_tmp(j+bond_sum,2)) .and.  &
                &   (coord_tmp(i+bond_sum,3) .eq. coord_tmp(j+bond_sum,3))) then 
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 4
                  coord_tmp(ind_act,2) = coord_tmp(i+bond_sum,3)
                  coord_tmp(ind_act,3) = coord_tmp(i+bond_sum,2)
                  coord_tmp(ind_act,4) = coord_tmp(i+bond_sum,4)
                  coord_tmp(ind_act,5) = coord_tmp(j+bond_sum,4)
!
!     i3=j3 and i4=j4 (rest: i2, j2): out-of-plane
!
               else if ((coord_tmp(i+bond_sum,3) .eq. coord_tmp(j+bond_sum,3)) .and.  &
                &   (coord_tmp(i+bond_sum,4) .eq. coord_tmp(j+bond_sum,4))) then
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 4
                  coord_tmp(ind_act,2) = coord_tmp(i+bond_sum,3)
                  coord_tmp(ind_act,3) = coord_tmp(i+bond_sum,4)
                  coord_tmp(ind_act,4) = coord_tmp(i+bond_sum,2)
                  coord_tmp(ind_act,5) = coord_tmp(j+bond_sum,2)
!
!     i3=j2 and i4=i3 (rest: i2, j4): dihedral
!
               else if ((coord_tmp(i+bond_sum,3).eq. coord_tmp(j+bond_sum,2)) .and.  &
                &   (coord_tmp(i+bond_sum,4) .eq. coord_tmp(j+bond_sum,3))) then
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 3
                  coord_tmp(ind_act,2) = coord_tmp(i+bond_sum,2)
                  coord_tmp(ind_act,3) = coord_tmp(i+bond_sum,3)
                  coord_tmp(ind_act,4) = coord_tmp(i+bond_sum,4)
                  coord_tmp(ind_act,5) = coord_tmp(j+bond_sum,4)
!
!     i2=j3 and i3=j4 (rest: i4, i2): dihedral
!
               else if ((coord_tmp(i+bond_sum,2) .eq. coord_tmp(j+bond_sum,3)) .and.  &
                &   (coord_tmp(i+bond_sum,3) .eq. coord_tmp(j+bond_sum,4))) then
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 3
                  coord_tmp(ind_act,2) = coord_tmp(j+bond_sum,2)
                  coord_tmp(ind_act,3) = coord_tmp(j+bond_sum,3)
                  coord_tmp(ind_act,4) = coord_tmp(j+bond_sum,4)
                  coord_tmp(ind_act,5) = coord_tmp(i+bond_sum,4)
!
!     i3=j4 and i4=j3 (rest: i2 and j2): dihedral
!
               else if ((coord_tmp(i+bond_sum,3) .eq. coord_tmp(j+bond_sum,4)) .and.  &
                &   (coord_tmp(i+bond_sum,4) .eq. coord_tmp(j+bond_sum,3))) then
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 3
                  coord_tmp(ind_act,2) = coord_tmp(i+bond_sum,2)
                  coord_tmp(ind_act,3) = coord_tmp(i+bond_sum,3)
                  coord_tmp(ind_act,4) = coord_tmp(i+bond_sum,4)
                  coord_tmp(ind_act,5) = coord_tmp(j+bond_sum,2)
!
!     i3=j2 and i2=j3 (rest: i4 and j4): dihedral
!
               else if ((coord_tmp(i+bond_sum,3) .eq. coord_tmp(j+bond_sum,2)) .and.  &
                &   (coord_tmp(i+bond_sum,2) .eq. coord_tmp(j+bond_sum,3))) then
                  ind_act = ind_act + 1
                  coord_tmp(ind_act,1) = 3
                  coord_tmp(ind_act,2) = coord_tmp(i+bond_sum,4)
                  coord_tmp(ind_act,3) = coord_tmp(i+bond_sum,3)
                  coord_tmp(ind_act,4) = coord_tmp(i+bond_sum,2)
                  coord_tmp(ind_act,5) = coord_tmp(j+bond_sum,4)
               end if
            end do
         end do      
!
!     The number of internal coordinates
!   
         nat6=ind_act
         allocate(coord_def(nat6,5))
         coord_def(1:nat6,:)=coord_tmp(1:nat6,:)
      end if
! 
!     Distinguish between debug and normal calculation
!
      coord=coord/bohr
      call gradient(coord,e_evb,g_evb,1,1)  ! else calculate the usual gradient
!
!     If desired, calculate the Wilson matrix the structure and print it to file
!
!     Determine internal coordinates, if no pes_topol keyword is used

      if (print_wilson) then
         if (.not. pes_topol) then
            call init_int("wilson",0,0,1)
         end if
         if (allocated(internal)) deallocate(internal)
         if (allocated(B_mat)) deallocate(B_mat)
         allocate(internal(nat6),B_mat(nat6,3*natoms))
         if (wilson_mode .eq. 2) then
            if (allocated(dB_mat)) deallocate(dB_mat)
            allocate(dB_mat(nat6,3*natoms,3*natoms))
         end if

         call xyz_2int(coord,internal,natoms)
!
!     If PES TOPOL option is activated, write internal coordinate definition and values 
!         to file pes_topol.dat
!
         write(214,'(a,i6,a)') "#-------- structure",ref_count,": -------------------------"
         do i=1,nat6
            if (coord_def(i,1) .eq. 1) then
               if (topol_bonds .eq. "VDW") then               
                  write(214,'(i9,a,i9,i9,a,f14.7)') i,"  bond  ",coord_def(i,2),coord_def(i,3),&
                     & "                  ",internal(i)
               else if (topol_bonds .eq. "EHT") then
                  write(214,'(i9,a,i9,i9,a,2f14.7)') i,"  bond  ",coord_def(i,2),coord_def(i,3),&
                     & "                  ",internal(i),bond_order(i)
               end if
            else if (coord_def(i,1) .eq. 2) then
               write(214,'(i9,a,i9,i9,i9,a,f14.7)') i,"  angle ",coord_def(i,2),coord_def(i,3), &
                  & coord_def(i,4),"         ",internal(i)
            else if (coord_def(i,1) .eq. 3) then
               write(214,'(i9,a,i9,i9,i9,i9,f14.7)') i,"  dihed.",coord_def(i,2),coord_def(i,3), &
                  & coord_def(i,4),coord_def(i,5),internal(i)
            else if (coord_def(i,1) .eq. 4) then
               write(214,'(i9,a,i9,i9,i9,i9,f14.7)') i,"  oop.  ",coord_def(i,2),coord_def(i,3), &
                  & coord_def(i,4),coord_def(i,5),internal(i)
            end if
         end do

         write(215,'(a,i6,a)') "#-------- structure",ref_count,": -------------------------"
         call calc_wilson(coord,internal,B_mat)
         do i=1,nat6
            do j=1,3*natoms
               write(215,'(a,i9,a,i9,a,f16.8)') "    ",i,"         ",j,"           ",B_mat(i,j)
            end do
         end do
!
!     If ordered, calculate and print also the Wilson matrix derivative!
!
         if (wilson_mode .eq. 2) then
            write(216,'(a,i6,a)') "#-------- structure",ref_count,": -------------------------"
            call calc_dwilson(coord,internal,dB_mat)
            do i=1,nat6
               do j=1,3*natoms
                  do k=1,3*natoms
                     write(216,'(a,i9,a,i9,a,i9,a,f16.8)') "    ",i,"        ",j,"         ",k, &
                          & "           ",dB_mat(i,j,k)
                  end do
               end do
            end do
         end if
      end if

!
!     Calculate and write the energies of the single qmdffs!
!
      if ((.not. treq) .and. (.not. pot_ana) .and. (nqmdff .gt. 1) ) then
         call eqmdff(coord*bohr,e_qmdff1,e_qmdff2)
         write(99,*) e_qmdff1,e_qmdff2
      end if
      if (.not. pes_topol) then
         write(44,*) e_evb
         write(45,*) "Structure",ref_count,":"
         do k=1,natoms
            write(45,*) g_evb(:,k)
         end do
      end if
   end do
   if (print_wilson) then
      close(215)
      if (wilson_mode .eq. 2) close(216)
   end if
   if (pes_topol) close(214)

   if (treq) close(48)
   close(172)
!
!    For debugging, print the QMDFF components to two files (one for each QMDFF)
!
   close(166)
   close(192)
   if (.not. pes_topol) then
      close(44)
      close(45)
      close(99)
   end if
   if (treq) close(48)
   if (int_grad_plot) close(192)
   if (do_debug) then
!
!    Components of first QMDFF
!
      open(unit=99,file="debug_qmdff1.dat",status="unknown")
      do i=1,struc_no
         do j=1,comp_no
            write(99,'(e15.7)',advance="no") qmdff_parts(1,i,j)
         end do
         write(99,*) " "
      end do
      close(99)

!
!    Components of second QMDFF
!
      open(unit=99,file="debug_qmdff2.dat",status="unknown")
      do i=1,struc_no
         do j=1,comp_no2
            write(99,'(e15.7)',advance="no") qmdff_parts(2,i,j)
         end do
         write(99,*) " "
      end do
      close(99)

!
!    gnuplot script for first QMDFF
!
      open(unit=99,file="debug1_plot.gnu",status="unknown")
      write(99,'(a,a,a)') "plot 'debug_qmdff1.dat' u 0:1 w l title '",trim(parts_labels(1,1)),"' \"
      do i=2,comp_no
         write(99,'(a,i4,a,a,a)')  ", '' u 0:",i," w l title '",trim(parts_labels(1,i)),"' \"
      end do
      write(99,*) " "
      write(99,*) "pause -1"
      close(99)

!
!    gnuplot script for first QMDFF
!
      open(unit=99,file="debug2_plot.gnu",status="unknown")
      write(99,'(a,a,a)') "plot 'debug_qmdff2.dat' u 0:1 w l title '",trim(parts_labels(2,1)),"' \"
      do i=2,comp_no2
         write(99,'(a,i4,a,a,a)')  ", '' u 0:",i," w l title '",trim(parts_labels(2,i)),"' \"
      end do
      write(99,*) " "
      write(99,*) "pause -1"
      close(99)
!
!    close output files 
!
      close(288)
      close(289)
      close(290)
   end if
   if (.not. pes_topol) then
      write(*,*) "Energies written to 'energies.dat', gradients written to 'gradients.dat'."
      if ((.not. treq) .and. (.not. pot_ana) .and. (nqmdff .gt. 1) ) then
         write(*,*) "Energies of QMDFFs written to 'single_qmdff.dat'."
      end if
   end if
   goto 389
end if



!
!    Translate the read in coordinates (x,y,z arrays) to the coord array
!
allocate(coord(3,natoms))
do i=1,natoms
   coord(1,i)=x(i)
   coord(2,i)=y(i)
   coord(3,i)=z(i)
end do
!
!     Convert coordinates to bohr
!
coord=coord/bohr

write(15,*) "The structure for which the calculations shall be perfomed contains the "
write(15,*) "cartesian coordinates:"
write(15,*)
do i=1,nats
   write(15,*)i, coord(:,indi(i))
end do

!
!   If optimations shall be done, start with them
!   a simple Steepest descent algorithm is used here
!

if (trim(jobtype) .eq. "OPT_MIN" .or. trim(jobtype) .eq. "OPTFREQ") then
   write(*,*) "Initializing geometry optimization..."
   opt_min=.true.
   write(15,*)
   write(15,*) "---GEOMETRY OPTIMIZATION:-------"
   write(15,*)
   write(15,*) "Initializing geometry optimization..."

   call geoopt(coord) 

   write(*,*) "Finished!"
   write(*,*) "Succession of optimized structures written to geoopt.xyz"
   write(*,*) "Final optimized structure written to opt_final.xyz"
   if (coord_vasp) then
      write(*,*) "Optimized structure also written to CONTCAR (VASP format)"
   end if
   write(*,*) "Calculate gradient and energy of optimized structure..."
else  
   write(*,*) "Calculate energy and gradient of given structure..."
end if 
write(15,*)
write(15,*) "---RESULTS:-------"
write(15,*)
!
!   gradient calculations
!
allocate(g_evb(3,natoms))
if (grad) then
   call gradient(coord,e_evb,g_evb,1,1)
   if (opt_min) then
      write(15,*)"*The optimized structure has the coordinates:"
      do i=1,nats
         write(15,*) i,coord(:,indi(i))
      end do
   end if
   write(15,*)
   write(15,*)"*The energy (Hartree) of the structure is:"
   write(15,*) e_evb
   write(15,*) 
   write(15,*)"*The gradient (Hartree/bohr) of the structure is:"
   write(15,*) 
   do i=1,nats
      write(15,*) i, g_evb(:,indi(i))
   end do
   write(15,*)
end if
!
!   hessian calculation and normal mode analysis
!
if (frequency) then
!
!    Default: activate intensities
!
   calc_freq_int = .true.

   write(*,*) "Calculate frequencies of the structure..."
   write(15,*)"*The hessian of the structure is:"
   write(15,*)
   allocate(hess(3*nats,3*nats))
   allocate(h_out(9*nats*nats),freqs(3*nats))
   allocate(eigvecs(3*nats,3*nats))
   call hessevb(coord,hess)
!   write out the hessian
!  hess=hess/bohr/bohr
  do j=1,3*nats
      do m=1,3*nats
         h_out((j-1)*3*nats+m)=hess(m,j)
      end do
   end do
   maxline= int(size(h_out)/5)
   do j=0,maxline-1
      write(15,'((F12.7),(F12.7),(F12.7),(F12.7),(F12.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3),h_out(j*5+4),h_out(j*5+5)
   end do
   if (size(h_out)-maxline*5 .eq.4) then
      write(15,'((F12.7),(F12.7),(F12.7),(F12.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3),h_out(j*5+4)
   else if (size(h_out)-maxline*5 .eq.3) then
      write(15,'((F12.7),(F12.7),(F12.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3)
   else if (size(h_out)-maxline*5 .eq.2) then
      write(15,'((F12.7),(F12.7))') &
          &h_out(j*5+1),h_out(j*5+2)
   else if (size(h_out)-maxline*5 .eq.1) then
      write(15,'((F12.7))') &
          &h_out(j*5+1)
   end if
   write(15,*) 
   
!
!     For QMDFF generations: write out the hessian in a pseudo orca.hess
!     file, if ordered by keyword ORCA_FAKE
!
   if (orca_fake) then
      open(unit=16,file="explore.hess",status="unknown")
      write(16,*)
      write(16,*) "$orca_hessian_file"
      write(16,*)
      write(16,*) "$act_atom"
      write(16,*) natoms
      write(16,*) 
      write(16,*) "$act_coord"
      write(16,*) natoms
      write(16,*) 
      write(16,*) "$act_energy"
      write(16,*) e_evb
      write(16,*)
      write(16,*) "$hessian"
      write(16,*) 3*natoms

      maxline= int((3*natoms)/5)
      do j=0,maxline-1
         write(16,*) "   ",(j*5+k-1,"    ",k=1,5)
         do k=1,3*natoms 
            write(16,'(i5,5e17.7)') k-1,hess(k,j*5+1:j*5+5)
         end do
      end do

      if ((3*natoms-maxline*5) .eq. 4) then 
         write(16,*) "   ",(maxline*5+k-1,"    ",k=1,4)
         do k=1,3*natoms 
            write(16,'(i5,5e17.7)') k-1,hess(k,maxline*5+1:maxline*5+4)
         end do
      else if ((3*natoms-maxline*5) .eq. 3) then
         write(16,*) "   ",(maxline*5+k-1,"    ",k=1,3)
         do k=1,3*natoms
            write(16,'(i5,5e17.7)') k-1,hess(k,maxline*5+1:maxline*5+3)
         end do
      else if ((3*natoms-maxline*5) .eq. 2) then
         write(16,*) "   ",(maxline*5+k-1,"    ",k=1,2)
         do k=1,3*natoms
            write(16,'(i5,5e17.7)') k-1,hess(k,maxline*5+1:maxline*5+2)
         end do
      else if ((3*natoms-maxline*5) .eq. 1) then
         write(16,*) "   ",(maxline*5+k-1,"    ",k=1,1)
         do k=1,3*natoms
            write(16,'(i5,5e17.7)') k-1,hess(k,maxline*5+1:maxline*5+1)
         end do
      end if

      write(16,*) 
      write(16,*) "#"
      write(16,*) "# The atoms: label  mass x y z (in bohrs)"
      write(16,*) "#"
      write(16,*) "$atoms"
      write(16,*) natoms
      do j=1,natoms
         write(16,*) name(j),mass(j),coord(:,j)
      end do
      write(16,*) 
!
!    write line with final single point energy to insert it into orca.out file
!
      write(16,*) "The following line(s) should be included into the orca.out file &
              &(at the end):"
      write(16,*)
      write(16,'(a)') "-------------------------   --------------------"
      write(16,'(a,f20.12)') "FINAL SINGLE POINT ENERGY ", e_evb
      write(16,'(a)') "-------------------------   --------------------"
      close(16)
      write(*,*)
      write(*,*) "Fake orca-output file explore.hess written."
      write(*,*) "Insert last line of explore.hess into orca.out file to get &
              &the correct energy!"
      write(*,*) "This file can be used for QMDFF generations!"
      write(*,*) 
   end if
!
!     estimate atom numbers from element symbols
! 
   allocate(atind(nats))
   do i=1,natoms
      call elem(name(i),atind(i))
   end do
!
!   do the normal coordinate analysis
!
   call calc_freq(hess,freqs,eigvecs,.true.)
!
!     Write out the normal modes in molden formate for better handling
!
!   call g98fake("explore.molden",natoms,atind,coord,freqs,eigvecs,eigvecs)
!   write(15,*) "File explore.molden with normal mode spectum written."

   deallocate(hess)
   deallocate(eigvecs)
   write(*,*) "File explore.molden with normal mode spectrum written."
   write(*,*) "Finished!"
   
end if
!
!     TS and IRC:
!
!     Do all calculations in internal coordinates!
!     Specify them automatically and the full set of 3N-6 coordinates
!     unless read in of internal coordinates is ordered explicitly
!
if (ts_opt .or. calc_irc) then
   if (.not. read_coord) then
      write(*,*) "Automatic generation of coordinates is unstable at the moment."
      write(*,*) "Please add the READ_COORD keyword!"
      call fatal
   end if
   if (.not. dg_evb .and. .not. treq) then
      int_mode=3
      if (read_coord) int_mode=1
      call init_int(xyzfile,1,1,int_mode)
   end if
end if
!
!     TS optimization: call the respective subroutine
!
if (ts_opt) then
   write(*,*) "Starting the TS optimization..."
   if (rank .eq. 0) then
      call opt_ts
   end if
!
!     Energy/gradient calculation of the optimized structure
!
   if (grad) then
      call gradient(coord,e_evb,g_evb,1,1)
      write(15,*)"*The optimized structure has the coordinates:"
      do i=1,nats
         write(15,*) i,coord(:,indi(i))
      end do
      write(15,*)
      write(15,*)"*The energy (Hartree) of the structure is:"
      write(15,*) e_evb
      write(15,*)
      write(15,*)"*The gradient (Hartree/bohr) of the structure is:"
      write(15,*)
      do i=1,nats
         write(15,*) i, g_evb(:,indi(i))
      end do
      write(15,*)
   end if


   write(*,*) "Finished!"
end if
!
!     IRC optimization: call the respective subroutine
!
if (calc_irc) then
   write(*,*) "Starting the Reaction path optimization (IRC)..."
   if (rank .eq. 0) then
      call opt_irc
   end if
   write(*,*) "Finished!"
end if

389 continue

call system_clock(time_int2)
time2_omp = omp_get_wtime()

!duration=real(time_int2)/1000-real(time_int1)/1000
duration=time2_omp-time1_omp

write(15,*)
write(15,*) "Calculation finished at: "
call timestamp ( )


write(*,*)  ".. .----. -- -.. --- -. . "
write(*,*) "Calculation successfully finished!"
write(*,*) "Output was written to explore.log"
write(*,'(A, F10.3, A)') " The calculation needed a time of",duration," seconds."
!
!     perform any final tasks before program exit
!

close(15)
deallocate(coord)
end program explore

