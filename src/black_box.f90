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
!     ##  program black_box  --  automatization utility tool         ##
!     ##                                                             ##
!     #################################################################
!
!     "black_box" is an on-top program for automatization of force 
!     field generation and upfollowing k(T) calculation.
!     It assumes that you want to use the TREQ coupling term.
!     A bunch of variables are predefined on standard values unless
!     you change them explicitly in the caracal.key file.
!

program black_box 
use general  !  general parameters
use qmdff    !  module with QMDFF global variables
use evb_mod  !  EVB coupling parameters

implicit none 
!     Loop variables
integer::i,j,k,l,m,o,p 
!     For read in of commands and settings
character(len=120)::string
character(len=20)::keyword
logical::prefix_key,exist
integer::next
character(len=120)::record
character(len=80)::prefix,line

character(len=40)::commarg ! string for command line argument
!     the MPI rank (here always 0)
integer::rank
!     total number of processors for called applications 
integer::nprocs_tot 
integer::nat3  ! 3*natoms
!     For read in of the IRC/MEP
logical::scan_path
character(len=:),allocatable::irc_prefix
character(len=80)::log_name,chk_name,fchk_name,xyz_name
character(len=2)::irc_direct  ! direction of the IRC: left2right or right2left
integer::formchk_stat ! exit status of formchk program
integer::sys_stat  ! general exit status for sys calls
integer::lastline,readstat ! if the document is over or corrupted
character(len=80)::a80
character(len=60)::adum,adum2  ! dummy character/string
character(len=120),allocatable::fchk_content(:)
integer::en_line,xyz_line,ch_line ! positions of energy/structure/charge blocks
integer::idum  ! dummy integer
real(kind=8)::rdum ! dummy real
integer::irc_num ! number of structures in the irc
integer::irc_lines  ! number of lines in the irc.fchk file
integer::ts_loc_a(1),ts_loc   ! index of TS structure in IRC
integer::charge,multi  ! charge and multiplicity of the system
integer::energylines,chargelines,coordlines ! for better readin of irc.fchk
integer,allocatable::at_charge(:)
real(kind=8),allocatable::irc_all_enxi(:),irc_all_xyz(:) ! IRC information storage
real(kind=8),allocatable::xi_vals(:)
real(kind=8),allocatable::irc_coords(:,:,:)
real(kind=8),allocatable::irc_energies(:)
real(kind=8)::read_dummy(5000) ! dummy array for scan read in
real(kind=8)::tmp1,tmp2     ! sorting of Xi values and energies 
real(kind=8),allocatable::tmp3(:,:)   ! sorting of IRC structures
!     For calculation of QM reference data of QMDFFs of both minima
logical::reactants_done,products_done  ! if this part was already done before
integer::min_procs,min_mem   ! processors and memors for minimum calculations
integer::ens_procs  ! processors for additional energy calculations
character(len=80)::link_gaussian   ! link to gaussian executable 
character(len=80)::link_orca  ! link to orca executable 
character(len=80)::link_qmdffgen  ! link to qmdffgen.x executable 
character(len=80)::link_rpmd  ! link to rpmd.x executable
character(len=80)::link_mpi  ! link to MPI starter (e.q. mpirun)
character(len=80),allocatable::gf_method(:)   ! commands for QM method
character(len=120)::gf_basis ! commands for QM reference basis set and other settings
integer::gf_met_words,gf_bas_words  ! number of single words for method and basis commands
character(len=80),allocatable::e_method(:)   ! commands for QM method
!                for additional energy calculations 
character(len=120)::e_basis ! command for QM methods basis for add. energy calcs
integer::e_met_words,e_bas_words  ! number of commands for additional energy calcs.
integer::geo_reactants,geo_products  ! if geoopts of reactants/products are already finished
real(kind=8),allocatable::reactants_xyz(:),products_xyz(:) ! for Gaussian: geoopt structures 
integer::minlines   ! structure section of minimum geoopt name.fchk files
!     For calculation of QM reference for RP-EVB coupling 
integer::rp_pts_min  ! minimum number of RP-EVB points (+N/2 will be added)
integer::rp_procs ! processors for calculation of RP-EVB reference
integer::pos_ts1,pos_ts2   ! number of TS structure for two calculation methods
integer::n_rp_pts    ! number of RP-EVB G+H reference points
integer,allocatable::rp_points(:)  ! positions of all RP-EVB reference points
real(kind=8),allocatable::irc_slope(:) ! energy derivatives of the IRC for all points
real(kind=8)::sum_slope_l,sum_slope_r  ! sum of all slopes left and right to the TS
real(kind=8),allocatable::irc_weighted(:)  ! weighted IRC path progress for RP points
integer::n_pts_left,n_pts_right  ! number of RP-EVB points left/right to the TS
real(kind=8)::pos_act  ! actual ideal position of RP-EVB point
real(kind=8),allocatable::pos_diff(:)  ! in order to determine the actual RP-EVB point
real(kind=8)::plot_min   ! minimum energy value along the IRC for plots
character(len=3),allocatable::rp_names(:)  ! subfolder names of RP-EVB reference calcs
integer::ncalcs_round  ! number of RP-point calculations to be done in parallel
integer::nrounds  ! number of rounds of RP-point calculations to be done
integer,allocatable::round_done(:)  ! indicator, which calculations of current run are finished
real(kind=8)::energy  ! actual energy to read in
integer::gradline,natsum,sum1,mincol,maxcol,sum2,maxline  ! parameters for RP-EVB read in
real(kind=8),allocatable::grad(:,:),hess(:,:)   ! input arrays for gradients and frequencies
real(kind=8),allocatable::grad1d(:),hess1d(:),hess1d2(:)  ! intermediate low dim arrays
real(kind=8),allocatable::xyz_out(:),g_out(:),h_out(:)  ! arrays for write out
logical::snforca
character(len=40)::names
logical::section
integer::ndays,nhours,nminutes,nseconds,tot_time
character(len=1)::software_irc  ! used software for the calculation of the given IRC
character(len=1)::software_gf  ! used software for QMDFF ref. and gradient/hessian calculations
character(len=1)::software_ens ! used software if the energies shall be corrected 
!   For additional calculation of energies 
integer::njobs_round  ! number of calculations to be done in this round
integer::njobs_done ! number of calculations already done in previous runs
logical::ens_extra  ! if energies shall be corrected or not
real(kind=8)::ens_reactants,ens_products  ! extra energy for reactant and product opt structures
!   For generation of dynamic folders 
character(len=80)::qmdff_en_line  ! line with QMDFF energies 
integer::temps(20)   ! list with simulation temperatures for RPMD
integer::temp_num   ! number of different temperatures
integer,allocatable::fragments(:)   ! for all atoms: to which fragment they belong 
integer,allocatable::fragsize(:),atfrags(:)  ! list of numbers of atoms per fragment 
integer,allocatable::frag1(:),frag2(:),frag3(:),frag4(:)  ! atoms for each fragment
integer::num_eds   ! number of molecules in the reactant structure
integer,allocatable::atind(:)  ! atomic indices for bondlength criterion
real(kind=8)::bondlentol  ! bondlength criterion tolerance
integer,allocatable::doubles(:,:) ! all bonds that are present in both structures
integer::n_form,n_break,n_double  ! number of bonds for different categories
integer,allocatable::breaks_tmp(:,:),forms_tmp(:,:),breaks(:,:),forms(:,:) ! bond arrays
character(len=20),allocatable::print_form(:),print_break(:)   ! print-out arrays for bonds
!   The RPMD setting variables 
real(kind=8)::dt,xi_min,xi_max
integer::equi_step,gen_step,nbins,umbr_step
logical::recross_nocheck,recross_mpi
character(len=10)::pmf_minloc
!   For start of RPMD calculations
character(len=4)::adum_nprocs
real(kind=8),allocatable::rates_mol(:,:),rates_molec(:,:)  ! reaction rates to read in 
real(kind=8),allocatable::log_kt_mol(:,:),log_kt_molec(:,:)  ! logarithms of calculated rates
real(kind=8),allocatable::avg_kt_mol(:),avg_kt_molec(:)  ! averaged reaction rates
real(kind=8),allocatable::var_kt_mol(:),var_kt_molec(:)  ! variances for reaction rates
real(kind=8),allocatable::avg_kt_print(:,:) ! averages of k(T) values for print out
real(kind=8),allocatable::ln_kt(:),onet(:)  ! ln(k(T)) and 1/T for Arrhenius plots
real(kind=8)::arrh_slope,arrh_y,corre  ! Arrhenius data from regression
integer::kt_error   ! if an error occured during a k(T) calculation
integer::kt_avg  ! number of k(T) calculations per temperature to be averaged
!   Manual read in of RPMD internal coordinates 
character(len=100)::coord_file 
logical::manual_ints
real(kind=8)::start,finish
integer::starti,finishi
!
!     No MPI parallelization is needed!
!
rank=0
!
!     print general program informations
!
if (rank .eq. 0) then
   call promo
   write(*,*) "PROGRAM BLACK_BOX: AUTOMATED EVB-QMDFF GENERATION AND KT CALCULATION"
end if

!
!    Measure the needed time for initialization of the system
!
call system_clock(time_int(1))


!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("black_box")
      stop
   else
      write(*,*) 
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'black_box.x -help' or 'black_box.x -h'."
      write(*,*) 
   end if
end if


call getkey(rank) 

write(*,*) "---- TASK 0: read in mandatory informations from file ----------"
write(*,*) 
!
!     Read the program, which was used for IRC calculation

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:15) .eq. 'IRC_SOFTWARE ') then
      read(record,*) prefix,software_irc
!     possible options: O (orca), G (gaussian)
      read_software=.true.
      select case (software_irc)
         case("O")
            write(*,*) "* Reference IRC is read in from ORCA output"
         case("G")
            write(*,*) "* Reference IRC is read in from GAUSSIAN output"
         case default
            read_software=.false.
         end select
   end if
end do
!
!     If the software wasn't defined in the caracal.key file, read it in manually!
!
if (.not. read_software) then
   write(*,*) "Please tell which QM software package was used for IRC calculation!"
   write(*,*) "Add the keyword 'IRC_SOFTWARE', with O (orca) or G (Gaussian)"
   call fatal
end if
!
!    If another method shall be used for the calculation of the energies 
!    along the path, activate it here 
!    Additionally, the method for this calculation needs to be assigned!
!
ens_extra=.false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:20) .eq. 'SEPARATE_ENERGY ') then
      write(*,*) "* The keyword SEPARATE_ENERGY is present, therefore energies"
      write(*,*) "   along the IRC path will be calculated separately."
      ens_extra=.true.
   end if
end do


!
!     Read in links to reference QM or other software to be used from keyfile
!

link_gaussian="g16"
link_orca="orca"
link_qmdffgen="qmdffgen.x"
link_rpmd="rpmd"
link_mpi="mpirun"

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:14))) .eq. 'SYMLINKS {' .or. trim(adjustl(record(1:14))) &
              &  .eq. 'SYMLINKS{') then

      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         record=adjustl(record)
         if (keyword(1:10) .eq. 'GAUSSIAN ') then
            read(record(9:120),'(a)') link_gaussian
            write(*,*) "* The Gaussian executable will be invoked with: ",trim(link_gaussian)
         else if (keyword(1:16) .eq. 'ORCA ') then
            read(record(5:120),'(a)') link_orca
            write(*,*) "* The orca executable will be invoked with: ",trim(link_orca)
         else if (keyword(1:10) .eq. 'QMDFFGEN ') then
            read(record(9:120),'(a)') link_qmdffgen
            write(*,*) "* The qmdffgen executable will be invoked with: ",trim(link_qmdffgen)
         else if (keyword(1:10) .eq. 'CALC_RATE ') then
            read(record(10:120),'(a)') link_rpmd
            write(*,*) "* The rpmd executable will be invoked with: ",trim(link_rpmd)
         else if (keyword(1:5) .eq. 'MPI ') then
            read(record(4:120),'(a)') link_mpi
            write(*,*) "* MPI runs will be started with: ",trim(link_qmdffgen)
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The SYMLINKS section has no second delimiter! (})"
            call fatal
         end if
      end do
      exit
   end if
end do

!
!    Read in the total number of processors to be accessible for all included types 
!    of reference as well as rpmd.x calculations 
!
nprocs_tot=1
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'NPROCS_TOTAL ') then
      read(record,*) prefix,nprocs_tot
   end if
end do
!
!     If not defined in file, read processor number from command line
!
if (nprocs_tot .lt. 2) then
   exist=.false.
   do while (.not. exist)
      write(*,*) "Total number of processors to be used for the whole calculation"
      write(*,*) "(At least two processors are required due to rpmd.x!):"
      read(*,*) nprocs_tot
      if (nprocs_tot .ge. 2) exist=.true.
   end do
end if

write(*,'(a,i2,a)') " * In total, ",nprocs_tot," processors will be available."
!
!    Read in the direction of the IRC!
!    Two possibilities: If the reaction from left to right (first to last structure)
!    shall be calculated, insert LEFT2RIGHT. If the other direction shall be calculated,
!    insert RIGHT2LEFT!
!
irc_direct="NN"
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:15) .eq. 'IRC_DIRECTION ') then
      read(record,*) prefix,adum
      call upcase(adum)
      if (trim(adum) .eq. "LEFT2RIGHT") irc_direct="LR"
      if (trim(adum) .eq. "RIGHT2LEFT") irc_direct="RL"
   end if
end do
!
!    If the IRC direction wasn't defined in the caracal.key file, read it in manually!
!
if (irc_direct .eq. "NN") then
   exist=.false.
   do while (.not. exist) 
      write(*,*) "Which direction of the IRC shall be used as reaction mechanism?"
      write(*,*) "1: Left to right (first to last) or 2: Right to left (last to first)?"
      read(*,*,iostat=readstat) idum
      if (readstat .ne. 0) cycle 
      if (idum .eq. 1) irc_direct="LR"
      if (idum .eq. 2) irc_direct="RL"
      if (irc_direct .eq. "LR" .or. irc_direct .eq. "RL") exist=.true.
   end do
end if

if (irc_direct .eq. "LR") then
   write(*,*) "* The IRC will be followed from left to right (first to last structure) during"
   write(*,*) "   the proposed reaction mechanism."
else if (irc_direct .eq. "RL") then
   write(*,*) "* The IRC will be followed from right to left (last to first structure) during"
   write(*,*) "   the proposed reaction mechanism."
end if

!
!    Read in filenames of MEP/IRC input informations!
!
if (software_irc .eq. "G") then
   exist=.false.
!
!    Try to read IRC filename keyword from caracal.key file!
!   
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'IRC_PREFIX ') then
         read(record,*) prefix,line
         if (allocated(irc_prefix)) deallocate(irc_prefix)
         if (.not. allocated(irc_prefix)) allocate(character(len=LEN(TRIM(line)))  &
               & :: irc_prefix)
         irc_prefix=trim(line)
!
!    Define the output filenames from the given filename prefix
!
         log_name=irc_prefix // ".log"
         fchk_name=irc_prefix // ".fchk"
         chk_name=irc_prefix // ".chk"
!
!    Check if the name.log and the name.chk files are present 
!
         inquire(file=log_name,exist=exist)
         if (.not. exist) then
            write(*,*) "The ",trim(log_name)," file could not been found!"
            call fatal
         end if
         inquire(file=chk_name,exist=exist)
         if (.not. exist) then
            write(*,*) "The ",trim(chk_name)," file could not been found!"
            call fatal
         end if
      end if
   end do
else if (software_irc .eq. "O") then
!
!    Try to read IRC filename keyword from caracal.key file!
!   
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'IRC_PREFIX ') then
         read(record,*) prefix,line
         if (allocated(irc_prefix)) deallocate(irc_prefix)
         if (.not. allocated(irc_prefix)) allocate(character(len=LEN(TRIM(line)))  &
               & :: irc_prefix)
         irc_prefix=trim(line)
!
!    Define the output filename from the given filename prefix
!
         xyz_name=irc_prefix // ".xyz"
!
!    Check if the name.log and the name.chk files are present 
!
         inquire(file=xyz_name,exist=exist)
         if (.not. exist) then
            write(*,*) "The ",trim(xyz_name)," file could not been found!"
            call fatal
         end if
      end if
   end do


end if

!
!    Read in number of processors and memory amount for minimizations 
!    Default: equal to total number of processors
!
ens_procs=nprocs_tot
min_procs=nprocs_tot
rp_procs=nprocs_tot
min_mem=nprocs_tot
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'MIN_NPROCS ') then
      read(record,*,iostat=readstat) prefix,min_procs
      if ((min_procs .le. 0) .or. (min_procs .gt. nprocs_tot)) then
         write(*,*) "ERROR! An invalid number of procesors has been given!"
         call fatal
      end if
   else if (keyword(1:16) .eq. 'RP_REF_NPROCS ') then
      read(record,*,iostat=readstat) prefix,rp_procs
      if ((min_procs .le. 0) .or. (min_procs .gt. nprocs_tot)) then
         write(*,*) "ERROR! An invalid number of procesors has been given!"
         call fatal
      end if
   else if (keyword(1:13) .eq. 'REF_MEMORY ') then
      read(record,*,iostat=readstat) prefix,min_mem
      if (min_mem .le. 0) then
         write(*,*) "ERROR! An invalid amount of memory has been given!"
         call fatal
      end if
   else if (keyword(1:13) .eq. 'ENS_NPROCS ') then
      read(record,*,iostat=readstat) prefix,ens_procs
      if (min_mem .le. 0) then
         write(*,*) "ERROR! An invalid amount of memory has been given!"
         call fatal
      end if   
   end if
end do

write(*,'(a,i2,a)') " * Minimizations will be parallelized with ",min_procs," processors."
write(*,'(a,i2,a)') " * RP-EVB reference calcs. will be parallelized with ",rp_procs," processors."
if (ens_extra) then
   write(*,'(a,i2,a)') " * Additional IRC energy calcs. will be parallelized with ",ens_procs," processors."
end if
write(*,'(a,i3,a)') " * Reference calculations will be done with ",min_mem," GB of memory."
!
!     If another number of minimum RP-EVB reference points along the path shall be used 
!     The default value will be 10
!
rp_pts_min=10
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'MIN_RP_POINTS ') then
      read(record,*,iostat=readstat) prefix,rp_pts_min
      write(*,*) "* Chosen minimum/offset number of RP-EVB points: ",rp_pts_min
   end if
end do
!
!    Read in reference method and basis set for minimzations and gradient
!    frequency reference calculations
!
allocate(gf_method(40))
gf_method="NN"
gf_basis="NN"
gf_met_words=0
gf_bas_words=0
charge=-100
multi=-100

section = .false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:14))) .eq. 'GRAD_FREQ {' .or. trim(adjustl(record(1:14))) &
              &  .eq. 'GRAD_FREQ{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)

         if (keyword(1:13) .eq. 'SOFTWARE ') then 
            read(record,*,iostat=readstat) prefix,software_gf
            if (software_gf .eq. "G") then
               write(*,*) "  -> Gaussian will be used for gradient+frequency calculations."
            else if (software_gf .eq. "O") then
               write(*,*) "  -> orca will be used for gradient+frequency calculations."
            else
               write(*,*) "ERROR! No valid software has been chosen for gradient+frequency calculations!"
               call fatal
            end if
         else if (keyword(1:13) .eq. 'METHOD ') then
            read(record,*,iostat=readstat) prefix,gf_method
            do k=1,40
               if (trim(gf_method(k)) .ne. "NN")  gf_met_words=gf_met_words+1
            end do
         else if (keyword(1:13) .eq. 'BASIS ') then
            read(record(10:120),'(a)') gf_basis
            gf_bas_words=1
!
!    If a orca IRC is given, charge and multiplicity of the system must be given as well!
!
         else if (keyword(1:13) .eq. 'CHARGE ') then
            read(record,*,iostat=readstat) names,charge
            if (readstat .ne. 0) then 
               write(*,*) "Correct format: CHARGE [total charge as integer]"
               call fatal
            end if
         else if (keyword(1:13) .eq. 'MULTI') then
            read(record,*,iostat=readstat) names,multi
            if (readstat .ne. 0) then
               write(*,*) "Correct format: MULTI [spin multiplicity as integer]"
               call fatal
            end if         
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The GRAD_FREQ section has no second delimiter! (})"
            call fatal
         end if
      end do
      exit
   end if
end do
if (.not. section) then
   write(*,*) "Please add a GRAD_FREQ section with settings for the gradient+frequency"
   write(*,*) " reference calculations!"
   call fatal
end if

!
!    If Gaussian is not used as reference method, require the input if charge and multiplicity
!
if (software_irc .ne. "G") then
   if (charge .eq. -100) then
      write(*,*) "Please give the charge of the system in the GRAD_FREQ section!"
      call fatal
   end if
   if (multi .eq. -100) then
      write(*,*) "Please give the spin multiplicity of the system in the GRAD_FREQ section!"
      call fatal
   end if
end if

!
!    Write out info about reference basis set and reference method
!

write(*,'(a)',advance="no") " * The used reference method is: "
do i=1,gf_met_words
   write(*,'(a,a)',advance="no") trim(gf_method(i))," "
end do
write(*,*) " "
write(*,'(a)',advance="no") " * The used reference basis set is: "
write(*,'(a,a)',advance="no") trim(gf_basis)," "
write(*,*) " "
!
!    Read in also the method and basis set for the additional energy
!    calculations if they are desired!
!
if (ens_extra) then
   allocate(e_method(40))
   e_method="NN"
   e_basis="NN"
   e_met_words=0
   e_bas_words=0

   section = .false.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:18))) .eq. 'ENERGY_EXTRA {' .or. trim(adjustl(record(1:18))) &
              &  .eq. 'ENERGY_EXTRA{') then
         section=.true.
         do j=1,nkey-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)

            if (keyword(1:13) .eq. 'SOFTWARE ') then
               read(record,*,iostat=readstat) prefix,software_ens
               if (software_ens .eq. "G") then
                  write(*,*) "  -> Gaussian will be used for extra energy calculations."
               else if (software_ens .eq. "O") then
                  write(*,*) "  -> orca will be used for extra energy calculations."
               else
                  write(*,*) "ERROR! No valid software has been chosen for extra energy calculations!"
                  call fatal
               end if
            else if (keyword(1:13) .eq. 'METHOD ') then
               read(record,*,iostat=readstat) prefix,e_method
               do k=1,40
                  if (trim(gf_method(k)) .ne. "NN")  e_met_words=e_met_words+1
               end do
            else if (keyword(1:13) .eq. 'BASIS ') then
               read(record(10:120),'(a)') e_basis
               e_bas_words=1
            end if
            if (keyword(1:13) .eq. '}') exit
            if (j .eq. nkey-i) then
               write(*,*) "The ENERGY_EXTRA section has no second delimiter! (})"
               call fatal
            end if
         end do
         exit
      end if
   end do
   if (.not. section) then
      write(*,*) "Please add a ENERGY_EXTRA section with settings for the extra energy"
      write(*,*) " reference calculations!"
      call fatal
   end if

   write(*,'(a)',advance="no") " * The add. energy reference method is: "
   do i=1,e_met_words
      write(*,'(a,a)',advance="no") trim(e_method(i))," "
   end do
   write(*,*) " "
   write(*,'(a)',advance="no") " * The add. energy reference basis set is: "
   write(*,'(a)') trim(e_basis)
end if


!
!    If needed, check if executables are really in the stated links 
!
if ((software_irc .eq. "G") .or. (software_gf .eq. "G") .or. (software_ens) .eq. "G") then
   sys_stat=system("command -v " // trim(link_gaussian) // " > sys.log")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Gaussian executable could not been found at ",trim(link_gaussian)," !"
      write(*,*) " Please add or alter the GAUSSIAN keyword (SYMLINKS section)!"
      call fatal
   end if
end if

if ((software_irc .eq. "O") .or. (software_ens .eq. "O") .or. (software_gf) .eq. "O") then
   sys_stat=system("command -v " // trim(link_orca) // " > sys.log")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! orca executable could not been found at ",trim(link_gaussian)," !"
      write(*,*) " Please add or alter the ORCA keyword (SYMLINKS section)!"
      call fatal
   end if
end if


sys_stat=system("command -v " // trim(link_qmdffgen) // " > sys.log")
if (sys_stat .ne. 0) then
   write(*,*) "ERROR! qmdffgen.x executable could not been found at ",trim(link_qmdffgen)," !"
   write(*,*) " Please add or alter the QMDFFGEN keyword (SYMLINKS section)!"
   call fatal
end if

sys_stat=system("command -v " // trim(link_rpmd) // " > sys.log")
if (sys_stat .ne. 0) then
   write(*,*) "ERROR! calc_rate.x executable could not been found at ",trim(link_rpmd)," !"
   write(*,*) " Please add or alter the CALC_RATE keyword (SYMLINKS section)!"
   call fatal
end if


!
!    If the automatically generated coordinate set for RPMD fails, a manual
!    alignment might be needed
!    In this case, this keyword activates, that a given internal coordinate file 
!    will be copied to all RPMD subfolders and the read_coord command will
!    be added to the RPMD keyfiles
!
temps=-10
temp_num=0
manual_ints=.false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'MANUAL_INTS ') then
      read(record,*,iostat=readstat) prefix,coord_file
      manual_ints=.true.
      inquire(file=coord_file, exist=exist)
      if (.not. exist) then
         write(*,*) "ERROR! You have activated MANUAL_INTS (manual read in of"
         write(*,*) " internal coordinates for rpmd.x) but the proposed coordinate"
         write(*,*) " file ",trim(coord_file)," is not present!"
         call fatal
      end if
   end if
end do



!
!    Read in the temperatures for CALC_RATE calculations
!
temps=-10
temp_num=0
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'RATE_TEMPS ') then
      read(record,*,iostat=readstat) prefix,temps 
      do j=1,20
         if (temps(j) .gt. 0.d0) temp_num=temp_num+1
      end do
   end if
end do
!
!     If no reference basis keyword was there, read temperatures from command line 
!
if (temp_num .eq. 0) then
   exist=.false.
   do while (.not. exist)
      write(*,*) "Give a list of temperatures for the RPMD calculations!"
      write(*,*) "At least one, maximum 20 temperatures might be listed."
      read(*,'(a)') adum
      exist=.true.
   end do
   read(adum,*,iostat=readstat) temps
   do j=1,20
      if (temps(j) .gt. 0.d0) temp_num=temp_num+1
   end do
end if
write(*,'(a,20i5)') " * RPMD calculations will be done at temperatures (K): ",temps(1:temp_num)
!
!    If more than one k(T) calculation shall be done per temperature, that
!    will be averaged in order to determine the RPMD convergence, activate it here 
!
kt_avg=1
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'KT_AVERAGE ') then
      read(record,*,iostat=readstat) prefix,kt_avg
   end if
end do

if (kt_avg .eq. 1) then
   write(*,*) "* One k(T) calculation will be done for each temperature."
else if (kt_avg .gt. 1) then
   write(*,'(a,i3,a)') " * Per temperature, ",kt_avg," RPMD k(T) calculations will be done!"
   write(*,*) "    Averages and variances will be calculated in the following."
else 
   write(*,*) "You have activated the KT_AVERAGE option, but no useful number of"
   write(*,*) " k(T) calculations per temperature was given! Check your input!"
   call fatal
end if

!
!    Read in all needed CALC_RATE.x keywords! 
!    Default values will be given for all of them; unless no keyword
!    is given into the caracal.key file, these values will be used 
!
nbeads=1
npaths=1
pre_exp=25.d0
rp_mid_tot=0.6d0
rp_mid_trans=0.1d0
dt=0.2d0
r_inf=10.d0
k_force=0.1d0
umbr_lo=-0.05d0
umbr_hi=1.10d0
umbr_dist=0.01d0
gen_step=10000
equi_step=10000
umbr_step=10000
umbr_traj=10
xi_min=-0.02
xi_max=1.10d0
nbins=5000
pmf_method="integration"
thermo="andersen"
andersen_step=60
nose_q=100.0
recr_equi=10000
child_tot=2000
child_interv=2000
child_point=100
child_evol=1000
energy_tol=500d0
err_max=1000
recross_check=.true.
recross_mpi=.false.
scan_path=.false. ! if a scan instead of an IRC was calculated
!
!    0) General dynamical parameters
!

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:15) .eq. 'RPMD_BEADS ') then
      read(record,*,iostat=readstat) prefix,nbeads
      if (readstat .ne. 0) then
         write(*,*) "Correct format: RPMD_BEADS [No. of beads]"
         call fatal
      end if
   else if (keyword(1:13) .eq. 'DELTAT ') then
      read(record,*,iostat=readstat) prefix,dt
      if (readstat .ne. 0) then
         write(*,*) "Correct format: DELTAT [time step (fs)]"
         call fatal
      end if
   else if (keyword(1:18) .eq. 'RPMD_EN_TOL ') then
      read(record,*,iostat=readstat) prefix,energy_tol
      if (readstat .ne. 0) then
         write(*,*) "Correct format: RPMD_EN_TOL [energy tolerance (kJ/mol)]"
         call fatal
      end if
   else if (keyword(1:13) .eq. 'MAX_ERROR ') then
      read(record,*,iostat=readstat) prefix,err_max
      if (readstat .ne. 0) then
         write(*,*) "Correct format: MAX_ERROR [Max. number of errors]"
         call fatal
      end if
   else if (keyword(1:13) .eq. 'SCAN_PATH ') then
      scan_path=.true.
      if (readstat .ne. 0) then
         write(*,*) "Correct format: SCAN_PATH "
         call fatal
      end if
   end if
end do

!
!   A) The NVT section
!
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'NVT {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'NVT{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!    The thermostat 
         if (keyword(1:20) .eq. 'THERMOSTAT ') then
            read(record,*,iostat=readstat) names,thermo
            if (readstat .ne. 0) then
               write(*,*) "Correct format: THERMOSTAT [thermostat]"
               call fatal
            end if
            call upcase(thermo)
!    How often the Andersen thermostat velocity rescaling is applied
         else if (keyword(1:16) .eq. 'ANDERSEN_STEP ') then
            read(record,*,iostat=readstat) names,andersen_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: ANDERSEN_STEP [No. of steps]"
               call fatal
            end if
!    The damping factor for the Nose-Hoover thermostat
         else if (keyword(1:13) .eq. 'NOSE_DAMP ') then
            read(record,*,iostat=readstat) names,nose_q
            if (readstat .ne. 0) then
               write(*,*) "Correct format: NOSE_DAMP [integer (>1)]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The NVT section has no second delimiter! (})"
            call fatal
         end if
      end do
      exit
   end if
end do


!
!   A) The MECHA section
!

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'MECHA {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'MECHA{') then
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:20) .eq. 'N_PATHS ') then
            read(record,*,iostat=readstat) prefix,npaths
         else if (keyword(1:13) .eq. 'DIST_INF ') then
            read(record,*,iostat=readstat) prefix,r_inf
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The MECHA section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do

!
!    B) The UMBRELLA section
!

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'UMBRELLA {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'UMBRELLA{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!
!     Now read in all other informaton of the MECHA section
!
!     The strength of the umbrella force constant
         if (keyword(1:20) .eq. 'BIAS ') then
            read(record,*,iostat=readstat) names,k_force
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BIAS [force constant (a.u.)]"
               call fatal
            end if
!     The interval of the reaction coordinate Xi where sampling shall happen
         else if (keyword(1:20) .eq. 'BONDS ') then
            read(record,*,iostat=readstat) names,umbr_lo,umbr_hi
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BONDS [lower Xi-val, upper Xi-val]"
               call fatal
            end if
!     Distance between two umbrella windows
         else if (keyword(1:20) .eq. 'DIST ') then
            read(record,*,iostat=readstat) names,umbr_dist
            if (readstat .ne. 0) then
               write(*,*) "Correct format: DIST [Xi-interval between two windows]"
               call fatal
            end if
!     Number of MD steps per window for structure generation
         else if (keyword(1:20) .eq. 'GEN_STEPS ') then
            read(record,*,iostat=readstat) names,gen_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: GEN_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of MD steps per window for umbrella equilibration
         else if (keyword(1:20) .eq. 'EQUI_STEPS ') then
            read(record,*,iostat=readstat) names,equi_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: EQUI_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of MD steps per window for umbrella sampling
         else if (keyword(1:20) .eq. 'SAMPLE_STEPS ') then
            read(record,*,iostat=readstat) names,umbr_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: SAMPLE_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of sampling trajectories per umbrella window
         else if (keyword(1:20) .eq. 'SAMPLE_TRAJS ') then
            read(record,*,iostat=readstat) names,umbr_traj
            if (readstat .ne. 0) then
               write(*,*) "Correct format: SAMPLE_TRAJS [Number of trajectories]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The UMBRELLA section has no second delimiter! (})"
            call fatal
         end if

      end do
   end if
end do

!
!     C) The PMF section
!
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'PMF {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'PMF{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!     Range along the reaction path variable for integration
         if (keyword(1:20) .eq. 'XI_RANGE ') then
            read(record,*,iostat=readstat) names,xi_min,xi_max
            if (readstat .ne. 0) then
               write(*,*) "Correct format: XI_RANGE [lower Xi-val, upper Xi-val]"
               call fatal
            end if
!     Number of bins in which the PMF shall be calculated
         else if (keyword(1:20) .eq. 'BINS ') then
            read(record,*,iostat=readstat) names,nbins
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BINS [number of PMF calculation bins]"
               call fatal
            end if
!     Which PMF calculation method shall be used
         else if (keyword(1:20) .eq. 'METHOD ') then
            read(record,*,iostat=readstat) names,pmf_method
            if (readstat .ne. 0) then
               write(*,*) "Correct format: METHOD [method identifiert]"
               call fatal
            end if
!     How the minimum along the PMF surface shall be located
         else if (keyword(1:20) .eq. 'MINLOC ') then
            read(record,*,iostat=readstat) names,pmf_minloc
            if (readstat .ne. 0) then
               write(*,*) "Correct format: MINLOC [method identifiert]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The PMF section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do

!
!    D) The RECROSS section
!
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'RECROSS {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'RECROSS{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!     Number of MD steps for equilibration trajectory
         if (keyword(1:20) .eq. 'EQUI_STEPS ') then
            read(record,*,iostat=readstat) names,recr_equi
            if (readstat .ne. 0) then
               write(*,*) "Correct format: EQUI_STEPS [Number of MD steps]"
               call fatal
            end if
!     Total number of child trajectories
         else if (keyword(1:20) .eq. 'CHILD_TOTAL ') then
            read(record,*,iostat=readstat) names,child_tot
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_TOTAL [Number of child trajs.]"
               call fatal
            end if
!     Number of parent MD steps between child spawnings
         else if (keyword(1:20) .eq. 'CHILD_INTERVAL ') then
            read(record,*,iostat=readstat) names,child_interv
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_INTERVAL [MD steps between childs]"
               call fatal
            end if
!     Number of children trajectories per spawning
         else if (keyword(1:20) .eq. 'CHILD_PERPOINT ') then
            read(record,*,iostat=readstat) names,child_point
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_PERPOINT [Number of children p.p.]"
               call fatal
            end if
!     Number of MD steps for each child trajectory
         else if (keyword(1:20) .eq. 'CHILD_STEPS ') then
            read(record,*,iostat=readstat) names,child_evol
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_STEPS [Number of MD steps per child]"
               call fatal
            end if
!     If recrossing shall be calculated with MPI
         else if (keyword(1:20) .eq. 'MPI ') then
            recross_mpi = .true.
!     If recrossing part shall not be checked for errors
         else if (keyword(1:20) .eq. 'NO_CHECK ') then
            recross_check = .false.
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The RECROSS section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do

!
!    print further informations about RPMD settings:
!
write(*,*) " -- further settings for RPMD k(T) calculations: --"
write(*,'(a,i3)') " * The number of RPMD beads in the system is: ",nbeads
write(*,'(a,f15.7)') " * The damping coefficient for RP-EVB is: ", pre_exp
write(*,'(a,2f15.7)') " * The interval for direct interpolation and its borders &
             & are: ", rp_mid_tot,rp_mid_trans
write(*,'(a,f15.7)') " * The MD timestep in fs is: ",dt
write(*,'(a,f15.7)') " * The asymptotic distance for reactants in Angstroms is: ",r_inf
write(*,'(a,f15.7)') " * The strength of the umbrella potential (a.u.): ",k_force 
write(*,'(a,2f15.7)') " * Borders of umbrella window distribution: ",umbr_lo,umbr_hi
write(*,'(a,f15.7)') " * Distance between two umbrella windows: ",umbr_dist
write(*,'(a,i8)') " * MD steps for initial structure generation: ",gen_step
write(*,'(a,i8)') " * MD steps for umbrella equilibration: ",equi_step
write(*,'(a,i8)') " * MD steps for umbrella samplings: ",umbr_step
write(*,'(a,i4)') " * Number of umbrella trajectories per window: ",umbr_traj
write(*,'(a,2f15.7)') " * Borders for umbrella integration (xi): ",xi_min,xi_max
write(*,'(a,i6)') " * Number of bins along the path for integration: ",nbins
write(*,'(a,a)') " * Method for generation of PMF profile: ",pmf_method
write(*,'(a,i8)') " * MD steps for recrossing equilibration: ",recr_equi
write(*,'(a,i6)') " * Total number of recrossing child trajectories: ",child_tot
write(*,'(a,i8)') " * Interval of steps between two child spawnings: ",child_interv
write(*,'(a,i8)') " * Number of childs per spawning points: ",child_point
write(*,'(a,i8)') " * MD steps for child trajectory evaluation: ", child_evol
write(*,'(a,f15.7)') " * Energy tolerance for MD errors: ",energy_tol
write(*,'(a,i8)') " * Maximum error number during RPMD calculation: ",err_max
if (scan_path) then
   write(*,*) " * A simple coordinate scan was taken was reference path calculation."
end if

!
!    TASK 1 : EXTRACT IRC / MEP
!
!    The structures and energies of the minimum energy reactionpath
!    will be read in
!    
!
!    Gaussian input: Only results of an IRC calculation will be accepted!
!      the output files of that calculcation (name.log and name.chk) need
!      to be present. 
!
write(*,*)
write(*,*) "---- TASK 1: extract IRC/MEP reference informations ----------"
write(*,*) 
if (software_irc .eq. "G") then
   formchk_stat=system("formchk " // chk_name // " " // fchk_name )
   if (formchk_stat .ne. 0) then
      write(*,*) "ERROR! The formchk command needed for extraction of the"
      write(*,*) "  name.chk file wasn't successful! Check your settings."
      call fatal
   end if    
!
!    Open the converted name.fchk file and read in the structures and energies 
!    of the IRC
!   
   open(unit=16,file=fchk_name,status="old")
   irc_num=0
   natoms=0
   charge=1000
   multi=1000
   irc_lines=0
!
!    A: For usual IRC calculations of the reaction path of interest
!
   if (.not. scan_path) then
      do 
         read(16,'(a)',iostat=lastline) a80
 
         if (lastline .ne. 0) exit
         irc_lines=irc_lines+1
!
!    Determine number of structures in the IRC
!
         if (index(a80,'point       1 Results for each geome') .ne. 0) then
            en_line=irc_lines  ! starting line for energy printout
            read(a80,*) adum,adum,idum,adum,adum,adum,adum,adum,adum,irc_num
!
!    only the half number because gradient norms are shown as well
!   
            irc_num=irc_num/2
            write(*,*) "Read in: numbers of structures in the IRC:",irc_num
         end if
!
!    Determine number of atoms per structure (i.e. in the system)
! 
         if (index(a80,'Atomic numbers                             I   N=') .ne. 0) then
            read(a80,*) adum,adum,adum,adum,natoms
            ch_line=irc_lines  ! starting line of charge printout 
            write(*,*) "Read in: number of atoms in the system:",natoms
         end if
!
!    Read in the charge of the system
!
         if (index(a80,'Charge                                     I') .ne. 0) then
            read(a80,*) adum,adum,charge
            write(*,*) "Read in: net charge of the system:",charge
         end if
!
!    Read in the multiplicity of the system
!
         if (index(a80,'Multiplicity                               I') .ne. 0) then
            read(a80,*) adum,adum,multi
            write(*,*) "Read in: spin multiplicity of the system:",multi
         end if
!
!    Determine the start of the geometry section
!
         if (index(a80,'point       1 Geometries') .ne. 0) then
            xyz_line=irc_lines  ! starting line of geometry printout
         end if

      end do
      close(16)
!
!    Check if all important parameters were read in
!
      if (irc_num .eq. 0) then
         write(*,*) "ERROR! No number of IRC structures has been read in!"
         call fatal
      else if (natoms .eq. 0) then
         write(*,*) "ERROR! No number of atoms was defined!"
         call fatal
      else if (charge .eq. 1000) then
         write(*,*) "ERROR! No net charge of the system has been read in!"
         call fatal
      else if (multi .eq. 1000) then
         write(*,*) "ERROR! No multiplicity of the system has been read in!"
         call fatal  
      end if
!
!     Calculate needed parameters for upfollowing readin
!
      energylines=int((irc_num*2)/5)
      chargelines=int(natoms/6)
      coordlines=int((irc_num*natoms*3)/5)
!
!     Allocate IRC information arrays
!
      allocate(irc_all_enxi(irc_num*2))
      allocate(irc_all_xyz(irc_num*natoms*3))
      allocate(at_charge(natoms))
!
!     Read in the irc.fchk file again in order to store all lines as array
!
      open(unit=16,file=fchk_name,status="old")
      allocate(fchk_content(irc_lines))
      do i=1,irc_lines
         read(16,'(a)') fchk_content(i) 
      end do
      close(16)
!
!     Extract all energy informations from the irc.fchk array
!
      irc_all_enxi=0.d0
      do i=1,energylines
         read(fchk_content(i+en_line),*) irc_all_enxi((i-1)*5+1:i*5) 
      end do 
!
!     If a last partial line remains, read it in
!
      if ((irc_num*2-energylines*5) .ne. 0) then
         read(fchk_content(energylines+en_line+1),*) irc_all_enxi(energylines*5+1:irc_num*2)
      end if

!
!    Extract all structural informations from the irc.fchk array
!
      irc_all_xyz=0.d0
      do i=1,coordlines
         read(fchk_content(i+xyz_line),*) irc_all_xyz((i-1)*5+1:i*5)
      end do
!
!     If a last partial line remains, read it in
!
      if ((irc_num*natoms*3-coordlines*5) .ne. 0) then
         read(fchk_content(i+xyz_line),*) irc_all_xyz(coordlines*5+1:irc_num*natoms*3)
      end if

!
!     Extract atomic charges in order to determine element symbols 
!
      i=0
      at_charge=0
      if (natoms .ge. 6) then
         do i=1,chargelines
            read(fchk_content(i+ch_line),*) at_charge((i-1)*6+1:i*6)
         end do
      end if
!
!     If a last partial line remains, read it in
!
      if ((natoms-chargelines*6) .ne. 0) then
         read(fchk_content(chargelines+ch_line+1),*) at_charge(chargelines*6+1:natoms)
      end if  
!
!     determine the elements for all atoms in the system from the atomic charge 
!
      do i=1,natoms
         call atomname(at_charge(i),name(i)) 
         call upcase(name(i))
      end do
!
!     Finally, store the values of IRC energies and Xi values in separate arrays
!
      allocate(xi_vals(irc_num))
      allocate(irc_energies(irc_num))
      do i=1,irc_num
         irc_energies(i)=irc_all_enxi(i*2-1)
         xi_vals(i)=irc_all_enxi(i*2)
      end do
!
!     Store the geometry in a better suited array
! 
      allocate(irc_coords(irc_num,natoms,3)) 
      do i=1,irc_num
         do j=1,natoms
            irc_coords(i,j,:)=irc_all_xyz((i-1)*(3*natoms)+(j-1)*3+1:(i-1)*(3*natoms)+j*3)
         end do
      end do  
!
!     B: For coordinate scans with Gaussian that generated the reaction path
! 
   else 
      open(unit=16,file=fchk_name,status="old")
      do
         read(16,'(a)',iostat=lastline) a80

         if (lastline .ne. 0) exit
         irc_lines=irc_lines+1
!
!    Determine number of atoms per structure (i.e. in the system)
! 
         if (index(a80,'Atomic numbers                             I   N=') .ne. 0) then
            read(a80,*) adum,adum,adum,adum,natoms
            ch_line=irc_lines  ! starting line of charge printout 
            write(*,*) "Read in: number of atoms in the system:",natoms
         end if
!        
!    Read in the charge of the system
!       
         if (index(a80,'Charge                                     I') .ne. 0) then
            read(a80,*) adum,adum,charge
            write(*,*) "Read in: net charge of the system:",charge
         end if 
!   
!    Read in the multiplicity of the system
!
         if (index(a80,'Multiplicity                               I') .ne. 0) then
            read(a80,*) adum,adum,multi
            write(*,*) "Read in: spin multiplicity of the system:",multi
         end if
!
!    Determine number of structures during the scan
!           
         if (index(a80,'Optimization Number of geometries          I   N= ') .ne. 0) then
            read(a80,*) adum,adum,adum,adum,adum,adum,irc_num
         end if

      end do
      close(16)
!
!    Check if all important parameters were read in
!
      if (irc_num .eq. 0) then
         write(*,*) "ERROR! No number of IRC structures has been read in!"
         call fatal
      else if (natoms .eq. 0) then
         write(*,*) "ERROR! No number of atoms was defined!"
         call fatal
      else if (charge .eq. 1000) then
         write(*,*) "ERROR! No net charge of the system has been read in!"
         call fatal
      else if (multi .eq. 1000) then
         write(*,*) "ERROR! No multiplicity of the system has been read in!"
         call fatal
      end if

!
!    Allocate final storage arrays
!
      allocate(xi_vals(irc_num))
      allocate(irc_energies(irc_num))
      allocate(irc_coords(irc_num,natoms,3))
      allocate(at_charge(natoms))
      open(unit=16,file=fchk_name,status="old")
      k=0
      do

         read(16,'(a)',iostat=lastline) a80

         if (lastline .ne. 0) exit


!
!     Read in the partial charges of the system
! 
         if (index(a80,'Atomic numbers                             I   N=') .ne. 0) then
            chargelines=int((natoms)/6)
            do i=1,chargelines
               read(16,*) at_charge((i-1)*6+1:i*6)
            end do
!
!     If a last partial line remains, read it in
!
            if ((natoms-chargelines*6) .ne. 0) then
               read(16,*) at_charge(chargelines*6+1:natoms)
            end if

!
!     determine the elements for all atoms in the system from the atomic charge 
!
            do i=1,natoms
               call atomname(at_charge(i),name(i))
               call upcase(name(i))
            end do
         end if

!
!     Read in the energies of the path structures
!

         if (index(a80,'Results for each geome ') .ne. 0) then
            k=k+1
            read(a80,*) adum,adum,adum,adum,adum,adum,adum,adum,adum,idum
            energylines=int((idum)/5)
            do i=1,energylines
               read(16,*) read_dummy((i-1)*5+1:i*5)
            end do
!
!     If a last partial line remains, read it in
!
            if ((idum-energylines*5) .ne. 0) then
               read(16,*) read_dummy(energylines*5+1:idum)
            end if
            irc_energies(k)=read_dummy(idum-1)
            xi_vals(k)=k 
         end if
!
!     Read in the structures itself an store them directly into global array
!
         if (index(a80,'Geometries               R   N= ') .ne. 0) then
            read(a80,*) adum,adum,adum,adum,adum,adum,idum
            coordlines=int((idum)/5)
            do i=1,coordlines
               read(16,*) read_dummy((i-1)*5+1:i*5)
            end do
!
!     If a last partial line remains, read it in
!
            if ((idum-coordlines*5) .ne. 0) then
               read(16,*) read_dummy(coordlines*5+1:idum)
            end if

            do i=1,natoms
               irc_coords(k,i,:)=read_dummy((idum-natoms*3)+(i-1)*3+1:(idum-natoms*3)+(i-1)*3+3)
            end do
            

         end if

      end do
      close(16)
   end if
   
!
!     Sort Xi coordinates, energies and structures such that the Xi variable
!     rises monotonically throughout the path
!
   allocate(tmp3(natoms,3))
   do i=1,irc_num
      tmp1=xi_vals(i)
      tmp2=irc_energies(i)
      tmp3(:,:)=irc_coords(i,:,:)
      j=i-1
      do while(j .gt. 0 .and. xi_vals(j) .gt. tmp1) 
         xi_vals(j+1)=xi_vals(j)
         irc_energies(j+1)=irc_energies(j)   
         irc_coords(j+1,:,:)=irc_coords(j,:,:)
         j=j-1
      end do
      xi_vals(j+1)=tmp1
      irc_energies(j+1)=tmp2 
      irc_coords(j+1,:,:)=tmp3
   end do
!
!     If the IRC shall be readed from right to left, invert all entries of the relevant arrays!
!
   if (irc_direct .eq. "RL") then
      write(*,*) "The IRC direction shall be right2left, therefore invert the upfollowing of the"
      write(*,*) " read in IRC data."
      do i=1,irc_num/2
!
!     First, store the rightest element in a temporary variable
!
         tmp1=xi_vals(irc_num-i+1)
         tmp2=irc_energies(irc_num-i+1)
         tmp3(:,:)=irc_coords(irc_num-i+1,:,:)
!
!     Second, Overwrite the rightest element with the leftest
!
         xi_vals(irc_num-i+1)=xi_vals(i)
         irc_energies(irc_num-i+1)=irc_energies(i)
         irc_coords(irc_num-i+1,:,:)=irc_coords(i,:,:)
!
!     Third, overwrite the leftest element with the temporary variable 
!
         xi_vals(i)=tmp1
         irc_energies(i)=tmp2
         irc_coords(i,:,:)=tmp3(:,:)
      end do
!
!     Swith the signs of all Xi values 
!
      do i=1,irc_num
         xi_vals(i)=-xi_vals(i)
      end do
   end if
!
!     Generate output folder for IRC informations
!
   inquire(file="mep_irc", exist=exist)
   if (.not. exist)  call system("mkdir mep_irc")
!
!     Output file for energies
!   
   open(unit=17,file="mep_irc/irc_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) irc_energies(i)
   end do
   close(17)
!
!     Output file for energies and reaction coordinates Xi
!
   open(unit=17,file="mep_irc/irc_xi_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) xi_vals(i),irc_energies(i)
   end do
   close(17)

!
!     Output file for xyz structures of the path
!
   open(unit=17,file="mep_irc/irc.xyz",status="replace")
   do i=1,irc_num
      write(17,*) natoms 
      write(17,*) 
      do j=1,natoms  
         write(17,*) name(j),irc_coords(i,j,:)*bohr
      end do   
   end do
   close(17)

else if (software_irc .eq. "O") then
!
!     Read in orca IRC file: only xyz trajectory with energies
!      as comments!
!
   open(unit=14,file=xyz_name,status="old")
   read(14,*) natoms
   i=1
   do 
      read(14,*,iostat=readstat)
      if (readstat .ne. 0) exit
      i=i+1
   end do
   close(14)
   irc_num=i/(natoms+2)
  
   allocate(xi_vals(irc_num))
   allocate(irc_energies(irc_num))
   allocate(irc_coords(irc_num,natoms,3))
 
   open(unit=14,file=xyz_name,status="old")
   do i=1,irc_num
      read(14,*) idum
      read(14,*) adum,adum,adum,adum,adum,irc_energies(i)
      do j=1,natoms
         read(14,*) name(j),irc_coords(i,j,:)
      end do
   end do
   close(14)
   irc_coords=irc_coords/bohr
!
!     Determine TS position for xi_vals determination: highest energy!
!
   ts_loc_a=maxloc(irc_energies)
   ts_loc=ts_loc_a(1)

   do i=ts_loc,1,-1
      xi_vals(ts_loc-i+1)=-real(i)+1
   enddo
   do i=ts_loc+1,irc_num
      xi_vals(i)=real(i-ts_loc)
   end do


!
!     If the IRC shall be readed from right to left, invert all entries of the relevant arrays!
!
   allocate(tmp3(natoms,3))
   if (irc_direct .eq. "RL") then
      write(*,*) "The IRC direction shall be right2left, therefore invert the upfollowing of the"
      write(*,*) " read in IRC data."
      do i=1,irc_num/2
!
!     First, store the rightest element in a temporary variable
!
         tmp1=xi_vals(irc_num-i+1)
         tmp2=irc_energies(irc_num-i+1)
         tmp3(:,:)=irc_coords(irc_num-i+1,:,:)
!
!     Second, Overwrite the rightest element with the leftest
!
         xi_vals(irc_num-i+1)=xi_vals(i)
         irc_energies(irc_num-i+1)=irc_energies(i)
         irc_coords(irc_num-i+1,:,:)=irc_coords(i,:,:)
!
!     Third, overwrite the leftest element with the temporary variable 
!
         xi_vals(i)=tmp1
         irc_energies(i)=tmp2
         irc_coords(i,:,:)=tmp3(:,:)
      end do
!
!     Swith the signs of all Xi values 
!
      do i=1,irc_num
         xi_vals(i)=-xi_vals(i)
      end do
   end if

!
!     Generate output folder for IRC informations
!
   inquire(file="mep_irc", exist=exist)
   if (.not. exist)  call system("mkdir mep_irc")
!
!     Output file for energies
!   
   open(unit=17,file="mep_irc/irc_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) irc_energies(i)
   end do
   close(17)
!
!     Output file for energies and reaction coordinates Xi
!
   open(unit=17,file="mep_irc/irc_xi_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) xi_vals(i),irc_energies(i)
   end do
   close(17)

!
!     Output file for xyz structures of the path
!
   open(unit=17,file="mep_irc/irc.xyz",status="replace")
   do i=1,irc_num
      write(17,*) natoms
      write(17,*)
      do j=1,natoms
         write(17,*) name(j),irc_coords(i,j,:)*bohr
      end do
   end do
   close(17)


end if

call system_clock(time_int(2))


!
!    TASK 2 : CALCULATE QMDFF REFERENCE DATA
!
!    Two minima will be considered: reactants (min1) and products (min2)
!    The first structure of the IRC will be starting point for the reactants 
!    calculation, and the last structure for the products calculation
!    In case of orca calculations, all tasks can be done in one calculation,
!    in case of Gaussian calculations, the geometry optimization needs 
!    to be done separately of the frequency and CM5 charge calculation
!
write(*,*)
write(*,*) "---- TASK 2: calculate QMDFF reference data of minima ----------"
write(*,*) 
!
!    Check if these calculations were already done in a run before 
!
inquire(file="reactants/calc_done", exist=reactants_done)
inquire(file="products/calc_done", exist=products_done)


if (.not. reactants_done .or. .not. products_done) then

   call system("mkdir reactants")
   call system("mkdir products")

   if (software_gf .eq. "G") then

write(*,*) "Due to different printout formates, therefore geometry optimizations"
write(*,*) "and QMDFF reference calculations need to be done in different runs!"
write(*,*)
  
!
!    PART A : Do geometry optimization for reactants and producs 
! 
!    First, write appropriate input files for geoopt calculations 
!
      do k=1,2
         if (k .eq. 1) then
            sys_stat=chdir("reactants")
         else if (k .eq. 2) then
            sys_stat=chdir("products")
         end if
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants/producs can't be accessed!"
            call fatal
         end if
         if (k .eq. 1) then
            open(unit=16,file="reactants_opt.com") 
            write(16,*) "%chk=reactants_opt.chk"
         else if (k .eq. 2) then
            open(unit=16,file="products_opt.com")
            write(16,*) "%chk=products_opt.chk"
         end if
         if (min_procs .lt. 10) then
            write(16,'(a,i1)') " %nprocshared=",min_procs
         else 
            write(16,'(a,i2)') " %nprocshared=",min_procs
         end if
         if (min_mem*min_procs .lt. 10) then
            write(16,'(a,i1,a)') " %mem=",min_mem,"GB"
         else 
            write(16,'(a,i2,a)') " %mem=",min_mem,"GB"
         end if
         write(16,'(a)',advance="no") " #p opt=maxcycles=200 "
         do i=1,gf_met_words
            write(16,'(a,a)',advance="no") trim(gf_method(i))," "
         end do
         write(16,'(a,a)',advance="no") trim(gf_basis)," "
         write(16,*) " " 
         write(16,*)
         if (k .eq. 1) then
            write(16,*) "Educts geometry optimization"
         else if (k .eq. 2) then
            write(16,*) "Products geometry optimization"
         end if
         write(16,*)
         write(16,'(i4,i4)') charge,multi
         do i=1,natoms
            if (k .eq. 1) then
               write(16,*) name(i),irc_coords(1,i,:)*bohr
            else if (k .eq. 2) then
               write(16,*) name(i),irc_coords(irc_num,i,:)*bohr
            end if
         end do
         
         write(16,*) 
         close(16)
         sys_stat=chdir("..")
      end do
!
!     Now start Gaussian Jobs for geoopts 
!
!     If less than two calculations might be done at once, start them subsequently
      if (nprocs_tot .lt. 2*min_procs) then
!
!     First reactants 
!
         write(*,*) "Start geometry optimization of reactants minimum."
         sys_stat=chdir("reactants") 
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " reactants_opt.com")            
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start geometry optimization of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// "  products_opt.com")          
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")

      else 
!
!     If two calculations can be done at once, start them parallel and test finishing
!
!
!     First reactants 
!
         write(*,*) "Start geometry optimization of reactants minimum."
         sys_stat=chdir("reactants")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " reactants_opt.com &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start geometry optimization of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " products_opt.com &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Now test each second, if the calculations are finished (then, the last line of the 
!     name.log file begins with 'Normal termination of Gaussian')
! 
         geo_reactants=0
         geo_products=0
         call sleep(1)
         do 
!
!     First reactants
!
            if (geo_reactants .ne. 1) then
               sys_stat=chdir("reactants")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory reactants can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -1 reactants_opt.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out") 
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'Normal termination of Gaussian') .ne. 0) then    
                     geo_reactants=1
                     write(*,*) "Geometry optimization of reactants minimum finished!"
                  else if (index(a80,'File lengths (MBytes)') .ne. 0) then
                     write(*,*) "ERROR! The optimization of the reactants minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
            end if
!
!     Then products 
!
            if (geo_products .ne. 1) then
               sys_stat=chdir("products")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory products can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -1 products_opt.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out")
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'Normal termination of Gaussian') .ne. 0) then       
                     geo_products=1
                     write(*,*) "Geometry optimization of products minimum finished!"
                  else if (index(a80,'File lengths (MBytes)') .ne. 0) then
                     write(*,*) "ERROR! The optimization of the products minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
            end if
            if ((geo_reactants .eq. 1) .and. (geo_products .eq. 1)) exit
            call sleep(1)
         end do
         sys_stat=system("rm reactants/tail.out")
         sys_stat=system("rm products/tail.out")
      end if
!
!    SECOND: Evaluate geometry optimizations and start QMDFF reference calculation
! 
!
!    Convert the name.chk output files of both minimizations 
!       
      sys_stat=chdir("reactants")
      formchk_stat=system("formchk reactants_opt.chk reactants_opt.fchk > formchk.log")
      if (formchk_stat .ne. 0) then
         write(*,*) "ERROR! Converting of reactants_opt.chk file failed!"
         call fatal
      end if 
      sys_stat=chdir("..")
      
      sys_stat=chdir("products")
      formchk_stat=system("formchk products_opt.chk products_opt.fchk > formchk.log")
      if (formchk_stat .ne. 0) then
         write(*,*) "ERROR! Converting of products_opt.chk file failed!"
         call fatal
      end if 
      sys_stat=chdir("..")
!
!     Read in both fchk files and get the optimized structures of minima
!
      allocate(reactants_xyz(3*natoms))
      allocate(products_xyz(3*natoms))
      minlines=int((natoms*3)/5)
!
!     First reactants minimum
!
      open(unit=16,file="reactants/reactants_opt.fchk",status="old")
      do
         read(16,'(a)',iostat=lastline) a80
         if (lastline .ne. 0) exit
         if (index(a80,'Current cartesian coordinates              R') .ne. 0) then      
            do i=1,minlines   
               read(16,*)  reactants_xyz((i-1)*5+1:i*5)
            end do
            if (natoms*3-minlines .ne. 0) then
               read(16,*)  reactants_xyz(minlines*5+1:3*natoms)
            end if
            exit
         end if
      end do
      close(16) 
!
!    Second products minimum
!
      open(unit=16,file="products/products_opt.fchk",status="old")
      do
         read(16,'(a)',iostat=lastline) a80
         if (lastline .ne. 0) exit
         if (index(a80,'Current cartesian coordinates              R') .ne. 0) then
            do i=1,minlines
               read(16,*)  products_xyz((i-1)*5+1:i*5)
            end do
            if (natoms*3-minlines .ne. 0) then
               read(16,*)  products_xyz(minlines*5+1:3*natoms)
            end if
            exit
         end if
      end do
      close(16)
!
!    PART B : Do reference data calculations for reactants and producs
! 
!    First, write appropriate input files for calculations 
!

      do k=1,2
         if (k .eq. 1) then
            sys_stat=chdir("reactants")
         else if (k .eq. 2) then
            sys_stat=chdir("products")
         end if
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants/producs can't be accessed!"
            call fatal
         end if
         if (k .eq. 1) then
            open(unit=16,file="reactants_ref.com")
            write(16,*) "%chk=reactants_ref.chk"
         else if (k .eq. 2) then
            open(unit=16,file="products_ref.com")
            write(16,*) "%chk=products_ref.chk"
         end if
         if (min_procs .lt. 10) then
            write(16,'(a,i1)') " %nprocshared=",min_procs
         else 
            write(16,'(a,i2)') " %nprocshared=",min_procs
         end if
         if (min_mem*min_procs .lt. 10) then
            write(16,'(a,i1,a)') " %mem=",min_mem,"GB"
         else 
            write(16,'(a,i2,a)') " %mem=",min_mem,"GB"
         end if
         write(16,'(a)',advance="no") " #p freq "
         do i=1,gf_met_words
            write(16,'(a,a)',advance="no") trim(gf_method(i))," "
         end do
         write(16,'(a,a)',advance="no") trim(gf_basis)," "
         write(16,*) " "
         write(16,*) "pop=(hirshfeld, nboread) density=current iop(7/33=1,6/80=1)" 
         write(16,*)
         if (k .eq. 1) then
            write(16,*) "Educts QMDFF reference calculation"
         else if (k .eq. 2) then
            write(16,*) "Products QMDFF reference calculation"
         end if
         write(16,*)
         write(16,'(i4,i4)') charge,multi
         do i=1,natoms
            if (k .eq. 1) then
               write(16,*) name(i),reactants_xyz((i-1)*3+1:i*3)*bohr
            else if (k .eq. 2) then
               write(16,*) name(i),products_xyz((i-1)*3+1:i*3)*bohr
            end if
         end do
         
         write(16,*) 
         WRITE(16,*) "$nbo bndidx $end"
         close(16)
         sys_stat=chdir("..")
      end do
!
!     Now start Gaussian Jobs for geoopts 
!
!     If less than two calculations might be done at once, start them subsequently
      if (nprocs_tot .lt. 2*min_procs) then
!
!     First reactants
!
         write(*,*) "Start QMDFF reference calculation of reactants minimum."
         sys_stat=chdir("reactants")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " reactants_ref.com")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start QMDFF reference calculation of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " products_ref.com")          
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")

      else 
!
!     If two calculations can be done at once, start them parallel and test finishing
!
!
!     First reactants
!
         write(*,*) "Start QMDFF reference calculation of reactants minimum."
         sys_stat=chdir("reactants")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " reactants_ref.com &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start QMDFF reference calculation of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_gaussian)// " products_ref.com &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The Gaussian16 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Now test each second, if the calculations are finished (then, the last line of the 
!     name.log file begins with 'Normal termination of Gaussian')
! 
         geo_reactants=0
         geo_products=0
         call sleep(1)
         do 
!
!     First reactants
!
            if (geo_reactants .ne. 1) then
               sys_stat=chdir("reactants")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory reactants can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -1 reactants_ref.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out") 
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'Normal termination of Gaussian') .ne. 0) then    
                     geo_reactants=1
                     write(*,*) "QMDFF reference calculation of reactants minimum finished!"
                     call system("touch calc_done")
                  else if (index(a80,'File lengths (MBytes)') .ne. 0) then
                     write(*,*) "ERROR! The QMDFF reference calculation of the reactants minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
            end if
!
!     Then products 
!
            if (geo_products .ne. 1) then
               sys_stat=chdir("products")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory products can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -1 products_ref.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out")
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'Normal termination of Gaussian') .ne. 0) then       
                     geo_products=1
                     write(*,*) "QMDFF reference calculation of products minimum finished!"
                     call system("touch calc_done")
                  else if (index(a80,'File lengths (MBytes)') .ne. 0) then
                     write(*,*) "ERROR! The QMDFF reference calculation of the products minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
               end if
            if ((geo_reactants .eq. 1) .and. (geo_products .eq. 1)) exit
            call sleep(1)
         end do
         sys_stat=system("rm reactants/tail.out")
         sys_stat=system("rm products/tail.out")
      end if
!
!     For later informations, write out xyz files with structures of reactants and products
!     minima!
!  
      open(unit=18,file="reactants/reactants_opt.xyz")
      open(unit=19,file="products/products_opt.xyz")

      write(18,*) natoms; write(18,*) 
      write(19,*) natoms; write(19,*)
      do i=1,natoms
         write(18,*) name(i),reactants_xyz((i-1)*3+1:i*3)*bohr
         write(19,*) name(i),products_xyz((i-1)*3+1:i*3)*bohr
      end do
  
      close(18)
      close(19)
  
   else if (software_gf .eq. "O") then
!
!     First, write input files for combined opt+freq calculations
!
      do k=1,2
         if (k .eq. 1) then
            sys_stat=chdir("reactants")
         else if (k .eq. 2) then
            sys_stat=chdir("products")
         end if
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants/producs can't be accessed!"
            call fatal
         end if
         if (k .eq. 1) then
            open(unit=16,file="reactants.inp")
         else if (k .eq. 2) then
            open(unit=16,file="products.inp")
         end if
         if (min_procs .gt. 1) then
            if (min_procs .lt. 10) then
               write(16,'(a,i1)') "! PAL",min_procs
            else
               write(16,'(a,i2)') "! PAL",min_procs
            end if
         end if
         write(16,'(a)',advance="no") "! opt freq "
         do i=1,gf_met_words
            write(16,'(a,a)',advance="no") trim(gf_method(i))," "
         end do
         write(16,'(a,a)',advance="no") trim(gf_basis)," "
         write(16,*)
         write(16,'(a,i4,i4)') "* xyz ", charge,multi
         do i=1,natoms
            if (k .eq. 1) then
               write(16,*) name(i),irc_coords(1,i,:)*bohr
            else if (k .eq. 2) then
               write(16,*) name(i),irc_coords(irc_num,i,:)*bohr
            end if
         end do

         write(16,*) "*"
         write(16,*) "%output"
         write(16,*) "    Print[ P_Hirshfeld ] 1"
         write(16,*) "end"
         close(16)
         sys_stat=chdir("..")
      end do
!
!     Now start orca Jobs for geoopts 
!
!     If less than two calculations might be done at once, start them subsequently
      if (nprocs_tot .lt. 2*min_procs) then
!
!     First reactants
!
         write(*,*) "Start QMDFF reference calculation of reactants minimum."
         sys_stat=chdir("reactants")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_orca)// " reactants.inp > reactants.log")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca5 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start QMDFF reference calculation of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_orca)// " products.inp > products.log")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca5 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")

      else
!
!     If two calculations can be done at once, start them parallel and test finishing
!
!
!     First reactants
!
         write(*,*) "Start QMDFF reference calculation of reactants minimum."
         sys_stat=chdir("reactants")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory reactants can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_orca)// " reactants.inp > reactants.log &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca5 calculation in 'reactants' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Then products 
!
         write(*,*) "Start QMDFF reference calculation of products minimum."
         sys_stat=chdir("products")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! Directory products can't be accessed!"
            call fatal
         end if
         sys_stat=system(trim(link_orca)// " products.inp > products.log &")
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca5 calculation in 'products' failed!"
            call fatal
         end if
         sys_stat=chdir("..")
!
!     Now test each second, if the calculations are finished (then, the last line of the 
!     name.log file begins with 'Normal termination of Gaussian')
! 
         geo_reactants=0
         geo_products=0
         call sleep(1)
         do
!
!     First reactants
!
            if (geo_reactants .ne. 1) then
               sys_stat=chdir("reactants")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory reactants can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -2 reactants.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out")
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'****ORCA TERMINATED NORMALLY****') .ne. 0) then
                     geo_reactants=1 
                     call system("touch calc_done")
                     write(*,*) "QMDFF reference calculation of reactants minimum finished!"
                  else if ((index(a80,'ERROR') .ne. 0) .or. (index(a80,'error') .ne. 0) .or. &
                      &  (index(a80,'Aborting the run') .ne. 0)) then
                     write(*,*) "ERROR! The reference calculaiton of the reactants minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
            end if
!
!     Then products 
!
            if (geo_products .ne. 1) then
               sys_stat=chdir("products")
               if (sys_stat .ne. 0) then
                  write(*,*) "ERROR! Directory products can't be accessed!"
                  call fatal
               end if
               sys_stat=system("tail -2 products.log > tail.out &")
               sys_stat=system("sleep 0.1")
               open (unit=16,file="tail.out")
               read(16,'(a)',iostat=readstat) a80
               if (readstat .eq. 0) then
                  if (index(a80,'****ORCA TERMINATED NORMALLY****') .ne. 0) then
                     geo_products=1
                     call system("touch calc_done")
                     write(*,*) "QMDFF reference calculation of products minimum finished!"
                  else if ((index(a80,'ERROR') .ne. 0) .or. (index(a80,'error') .ne. 0) .or. & 
                     &   (index(a80,'Aborting the run') .ne. 0)) then
                     write(*,*) "ERROR! The reference calculation of the products minimum failed!"
                     call fatal
                  end if
               end if
               close(16)
               sys_stat=chdir("..")
            end if
            if ((geo_reactants .eq. 1) .and. (geo_products .eq. 1)) exit
            call sleep(1)
         end do
         sys_stat=system("rm reactants/tail.out")
         sys_stat=system("rm products/tail.out")
      end if
   end if
else 
   write(*,*) "QMDFF reference information on both minima was already calculated!"
 
end if

!
!    TASK 2-B : CALCULATE QMDFF REFERENCE ENERGY WITH EXTRA REFERENCE
!
!    The optimized structures of both minima will be taken and energies 
!    of the energy correction reference will be calculated with them
!    These new energies will be written into a file in order to say 
!    qmdffgen that the energy shifts shall be corrected 
!

if (ens_extra) then
   write(*,*)
   write(*,*) "---- TASK 2-B: calculate QMDFF reference energy with extra reference -------"
   write(*,*)
   write(*,*) 
   ens_reactants=0.d0
   ens_products=0.d0
   do k=1,2
      if (k .eq. 1) then
         sys_stat=chdir("reactants")
      else if (k .eq. 2) then
         sys_stat=chdir("products")
      end if
      inquire(file="add_ens", exist=exist)
      if (.not. exist) call system("mkdir add_ens") 
     
      call chdir("add_ens")

!
!     Write orca input file for calculation
!     Use the external xyz file formate such that the optimized 
!     structures can be used directly!
!
      open(unit=15,file="orca.inp",status="unknown")
      write(15,'(a)',advance="no") "! "
      do j=1,e_met_words
         write(15,'(a,a)',advance="no") adjustl(trim(e_method(j)))," "
      end do
      write(15,'(a)') trim(e_basis)
      if (ens_procs .gt. 1) then
         write(15,'(a,i1)') "! PAL",ens_procs
      end if
      if (min_mem .le. 9) then
         write(15,'(a,i1,a)') "%maxcore ",min_mem,"000"
      else
         write(15,'(a,i2,a)') "%maxcore ",min_mem,"000"
      end if
      write(15,'(a,i3,a,i3,a)') "*xyzfile ",charge," ",multi," struc.xyz"
      close(15)
      if (software_gf .eq. "G") then
         if (k .eq. 1) then
            call system("cp ../reactants_opt.xyz struc.xyz")
         else if (k .eq. 2) then
            call system("cp ../products_opt.xyz struc.xyz")
      end if
      else if (software_gf .eq. "O") then
         if (k .eq. 1) then
            call system("cp ../reactants.xyz struc.xyz")
         else if (k .eq. 2) then
            call system("cp ../products.xyz struc.xyz")
         end if
      end if
!
!     Start the calculations if they are not already finished
!   
      inquire(file="done", exist=exist)
      if (.not. exist) then 
         if (k .eq. 1) then
            write(*,*) "Start extra energy calculation for optimized reactants.."
         else if (k .eq. 2) then
            write(*,*) "Start extra energy calculation for optimized products.."
         end if 
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca calculation in this folder failed!"
            call fatal
         end if
         sys_stat=system(trim(link_orca)// " orca.inp > orca.out")
         write(*,*) " --> done!"
         call system("touch done")
      else 
         if (k .eq. 1) then
            write(*,*) " Energy calculation already done for reactants!"
         else if (k .eq. 2) then
            write(*,*) " Energy calculation already done for products!"
         end if
      end if     
!
!     Evaluate the calculations
!
      open(unit=16,file="orca.out")
9     read(16,'(a)',end=28)a80
         if (index(a80,'FINAL SINGLE POINT ENERGY ') .ne. 0) then
            if (k .eq. 1) then
               read(a80,*) adum,adum,adum,adum,ens_reactants
            else if (k .eq. 2) then
               read(a80,*) adum,adum,adum,adum,ens_products 
            end if
         end if 
      goto 9
28    continue
      close(16,status='keep')

      call chdir("..")
      call chdir("..")
   end do
end if 



call system_clock(time_int(3))

!
!    TASK 3 : GENERATE THE QMDFFs
!
!    After the QM reference data needed for the QMDFF generation was calculated, generate 
!    both QMDFF force fields as first part of the full EVB-QMDFF force field
!
!    
!    Copy all needed input files to a new evb_qmdff folder 
!
write(*,*)
write(*,*) "---- TASK 3: Generate the QMDFF force fields from reference ----------"
write(*,*)

write(*,*) "Start the QMDFF generation by invoking qmdffgen.x ..."
inquire(file="evb_qmdff", exist=exist)
if (exist) call system("rm -r evb_qmdff")
call system("mkdir evb_qmdff")

if (software_gf .eq. "G") then
   sys_stat=system("cp reactants/reactants_ref.log evb_qmdff/reactants.log")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of reactant reference failed!"
      call fatal
   end if
   sys_stat=system("cp reactants/reactants_ref.chk evb_qmdff/reactants.chk")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of reactant reference failed!"
      call fatal
   end if
   sys_stat=system("cp products/products_ref.log evb_qmdff/products.log") 
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of products reference failed!"
   end if
   sys_stat=system("cp products/products_ref.chk evb_qmdff/products.chk")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of products reference failed!"
      call fatal
   end if
else if (software_gf .eq. "O") then
   sys_stat=system("cp reactants/reactants.log evb_qmdff/reactants.out")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of reactant reference failed!"
      call fatal
   end if
   sys_stat=system("cp reactants/reactants.hess evb_qmdff/reactants.hess")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of reactant reference failed!"
      call fatal
   end if
   sys_stat=system("cp products/products.log evb_qmdff/products.out")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of products reference failed!"
   end if
   sys_stat=system("cp products/products.hess evb_qmdff/products.hess")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Copying of products reference failed!"
      call fatal
   end if
end if

open(unit=15,file="evb_qmdff/qmdffgen.key",status="replace")
write(15,*) "software ",software_gf
write(15,*) "2qmdff reactants products"
close(15)
!
!     Invoke qmdffgen.x to generate the QMDFFs
!
sys_stat=chdir("evb_qmdff")
!
!     If an additional level of theory was used for calculation of energies, write the 
!     energies to the file that tells QMDFFGEN to activate the correction and read 
!     them in thereafter         
! 
!     Remove the file before, for the case that no correction shall be applied
!
inquire(file="min_energies.dat", exist=exist)
if (exist) call system("rm min_energies.dat")

if (ens_extra) then
   open(unit=16,file="min_energies.dat",status="new")
   write(16,*) ens_reactants
   write(16,*) ens_products
   close(16) 
end if 

sys_stat=system(trim(link_qmdffgen) // " qmdffgen.key > qmdffgen.log")
if (sys_stat .ne. 0) then
   write(*,*) "ERROR! Application of qmdffgen.x failed! Look into qmdff.log!"
   call fatal
end if
write(*,*) "QMDFF generation successful!"
sys_stat=chdir("..")

!
!    TASK 4 : CALCULATE TREQ GRADIENT+HESSIAN REFERENCE
!
!    For the parametrization of TREQ, gradients and frequencies need to be calculated
!    at a number of reference points along the path 
!    
!    First, the position of these points will be determined from the path topology
!    and the size of the system (more atoms=more reference points)
! 
!    Then, QM input data will be generated for each of this points, the calculations 
!    will be started and the results will be collected 
!

call system_clock(time_int(4))

write(*,*)
write(*,*) "---- TASK 4: Generate the reference data for TREQ  ----------"
write(*,*)
write(*,*) "First: determine positions of all RP-EVB points at which gradients"
write(*,*) " of the QM reference will be calculated."
!
!    Determine the position of the transition state (TS)
!    Take the highest energy and the lowest Xi coordinate value!
!    Only OK, if both criteria give the same structure of the IRC
!
pos_ts1=maxloc(irc_energies(:),dim=1)
pos_ts2=minloc(abs(xi_vals(:)),dim=1)

if (.not. scan_path) then
   if (pos_ts1 .ne. pos_ts2) then
      write(*,*) "ERROR! The structure of highest energy on the IRC has not the Xi"
      write(*,*) " value nearest to 0! Something seems to be broken with the IRC!"
      call fatal
   end if
   if ((pos_ts1 .eq. 1) .or. (pos_ts1 .eq. irc_num)) then
      write(*,*) "ERROR! The TS of your path is located at the right or left end"
      write(*,*) " of this path! This is not useful!"
      call fatal
   end if
end if
write(*,'(a,i3,a,i4,a)') " Position of the TS: IRC structure ",pos_ts1," of ",irc_num,"."
!
!    Determine number of RP-EVB reference points on the path
!    If the number of RP points is larger than half the IRC points, set is 
!    to this number
!
n_rp_pts=rp_pts_min+int(natoms/2)
if (n_rp_pts .gt. irc_num/2) then
   n_rp_pts=irc_num/2
end if

allocate(rp_points(n_rp_pts))
write(*,*) "Calculate slope of IRC energy profile..."
!
!    Determine numerical energy slope for each point along the IRC
!    The first and last point will be neglected because they only have one 
!    neighbor
!
allocate(irc_slope(irc_num-2))
do i=1,irc_num-2
   irc_slope(i)=0.5d0*(abs(irc_energies(i)-irc_energies(i+1))+abs(irc_energies(i+1)-&
               & irc_energies(i+2)))
end do

open(unit=15,file="mep_irc/irc_slopes.dat",status="replace")
do i=1,irc_num-2
   write(15,*) xi_vals(i+1),irc_slope(i)
end do
close(15)
write(*,*) "--> done! Values written to mep_irc/irc_slopes.dat."
write(*,*) "Weight IRC progress coordinate Xi with magnitude of energy slope..."
!
!    Sum up slopes right and left to the ts in order to get proper weighing parameters 
!
sum_slope_l=0.d0
sum_slope_r=0.d0
do i=1,pos_ts1-1 
   sum_slope_l=sum_slope_l+irc_slope(i)
end do
do i=pos_ts1-1,irc_num-2
   sum_slope_r=sum_slope_r+irc_slope(i)
end do

!
!    Fill array with weighted values to include the effect of energy slope 
!    along the IRC
!
allocate(irc_weighted(irc_num))
rdum=0.d0
do i=pos_ts1,1,-1
   if (i .eq. pos_ts1) then
      rdum=xi_vals(i)
   else 
      rdum=rdum+(xi_vals(i)-xi_vals(i+1))+0.5d0*(irc_slope(i)/(sum_slope_l)*xi_vals(1))
   end if 
   irc_weighted(i)=rdum
end do
rdum=0.d0
do i=pos_ts1,irc_num,1
   if (i .eq. pos_ts1) then
      rdum=xi_vals(i)
   else 
      rdum=rdum+(xi_vals(i)-xi_vals(i-1))+0.5d0*(irc_slope(i-1)/(sum_slope_r)*xi_vals(irc_num))
   end if
   irc_weighted(i)=rdum
end do
!
!    From the relative lengths of weighted paths, determine how many of the total
!    RP-EVB point number need to be right or left the TS
!
write(*,*) "Calculate integer positions of RP-EVB points along the TS."
n_pts_left=(n_rp_pts-3)*abs(xi_vals(1))/(abs(xi_vals(irc_num)- &
                & xi_vals(1)))
n_pts_right=(n_rp_pts-3)*abs(xi_vals(irc_num))/(abs(xi_vals(irc_num)- & 
                & xi_vals(1)))
!
!     If both numbers are too small due to rounding errors, increment them
!

if ((n_pts_left+n_pts_right) .lt. (n_rp_pts-4)) then
   n_pts_left=n_pts_left+1
   n_pts_right=n_pts_right+1
end if

!
!     If the sum of both is one too small, increase the smaller number by one 
!     (if both are equal, increase the left number: there will be more sampling)
!
if ((n_pts_left+n_pts_right) .lt. (n_rp_pts-3)) then
   if (n_pts_left .lt. n_pts_right) then
      n_pts_left=n_pts_left+1
   else if (n_pts_right .lt. n_pts_left) then
      n_pts_right=n_pts_right+1
   else  
      n_pts_left=n_pts_left+1
   end if
end if

!
!     Determine ideal positions of reference points on both sides of the TS,
!     referring to the scaled range 
!     Then look directly what point in the irc_weighted array is next to this 
!     ideal position; chose its index for the final RP-EVB point!
!
allocate(pos_diff(irc_num))
!
!     Left to the TS
!
rp_points(1)=1
do i=1,n_pts_left
   pos_act=i*irc_weighted(1)/(n_pts_left+1)
   pos_diff(:)=abs(irc_weighted(:)-pos_act)
   rp_points(2+n_pts_left-i)=minloc(pos_diff(:),dim=1)
end do
rp_points(2+n_pts_left)=pos_ts1
!
!     Right to the TS
!
do i=1,n_pts_right
   pos_act=i*irc_weighted(irc_num)/(n_pts_right+1)
   pos_diff(:)=abs(irc_weighted(:)-pos_act)
   rp_points(i+n_pts_left+2)=minloc(pos_diff(:),dim=1)
end do
rp_points(n_rp_pts)=irc_num

do i=1,n_rp_pts
   do j=1,n_rp_pts
      if (i .ne. j) then
         if (rp_points(i) .eq. rp_points(j)) then
            write(*,*) "ERROR! Two RP-EVB reference points seem to be located at the same"
            write(*,*) " structure! This is not useful!"
            write(*,*) " Decrease the value at RP_POINT_MIN or give more path structures!"
            call fatal
         end if
      end if
   end do
end do

write(*,'(a,i4,a)') " --> done! There are in total ",n_rp_pts," RP-EVB points (10+natoms/2)"
write(*,'(a,i4,a,i4,a)') " There are ",n_pts_left+1," points left and ",n_pts_right+1,&
           & " points right the TS."
write(*,*) " Positions written to rp_ref/rp_points.dat, plot rp_points.gnu to visualize!"

!
!     New folder for upfollowing RP-EVB reference data 
!
inquire(file="rp_ref", exist=exist)
if (.not. exist) call system("mkdir rp_ref")
!
!     Write positions of RP-EVB points to file for illustration
!     Additional file: energies of the IRC 
!

open(unit=15,file="rp_ref/rp_points.dat",status="replace")
write(15,*) "# No.   IRC-point   Xi-value   energy"
do i=1,n_rp_pts
   write(15,*) i,rp_points(i),xi_vals(rp_points(i)),irc_energies(rp_points(i))
end do
close(15)
open(unit=15,file="rp_ref/irc_ens.dat",status="replace")
do i=1,irc_num
   write(15,*) xi_vals(i),irc_energies(i)
end do
close(15)
!
!     Write gnuplot file for possible plot of RP-EVB points on the path
!
plot_min=minval(irc_energies)
open(unit=15,file="rp_ref/rp_points.gnu",status="replace")
write(15,*) "set title 'Positions of the RP-EVB reference points along the IRC'"
write(15,*) "set xlabel 'Reaction coordinate (Xi) / a.u.'"
write(15,*) "set ylabel 'energy (kJ/mol)'"
write(15,*) "set linestyle 1 linetype 2 linewidth 1 pointtype 5 pointsize 0.4 lc rgb '#006400'"
write(15,*) "plot 'irc_ens.dat' u 1:(($2-",plot_min," )*2650.50) with linespoints pointtype 4 & 
                & pointsize 0.7 lt 1 lc rgb 'red' title 'IRC energies', \"
write(15,*) " 'rp_points.dat' u 3:(($4-",plot_min," )*2650.50) with points pointsize 1.5 & 
                & pointtype 19 lc rgb 'blue' title 'RP-EVB points'"
write(15,*) "pause -1"
close(15)

write(*,*) "Generate folders and write input files for RP-EVB QM references.."
!
!     Generate QM input for all RP-EVB points, don't start them directly!
!
!
!     Fill array with names of subfolders for RP-EVB points
!
allocate(rp_names(n_rp_pts))
do i=1,n_rp_pts
   if (i .le. 9) then
      write(rp_names(i),'(i1)') i      
   else if (i .le. 99) then
      write(rp_names(i),'(i2)') i
   else if (i .le. 999) then
      write(rp_names(i),'(i3)') i
   else 
      write(*,*) "ERROR! No more than 999 RP-EVB points can be used!"
      call fatal
   end if
end do
!
!     If Gaussian shall be used for QM reference calculations
!
if (software_gf .eq. "G") then
   sys_stat=chdir("rp_ref")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! The folder rp_ref seems to be nonexistent!"
      call fatal
   end if 
   do i=1,n_rp_pts
      inquire(file=adjustl(trim(rp_names(i))), exist=exist)
      if (.not. exist) call system("mkdir " // adjustl(trim(rp_names(i))))
      open(unit=16,file=adjustl(trim(rp_names(i))) // "/gauss.com",status="replace")
      write(16,*) "%chk=gauss.chk"
      if (min_procs .lt. 10) then
         write(16,'(a,i1)') " %nprocshared=",rp_procs
      else
         write(16,'(a,i2)') " %nprocshared=",rp_procs
      end if
      if (min_mem*min_procs .lt. 10) then
         write(16,'(a,i1,a)') " %mem=",min_mem,"GB"
      else
         write(16,'(a,i2,a)') " %mem=",min_mem,"GB"
      end if
      write(16,'(a)',advance="no") " #p freq "
      do j=1,gf_met_words
         write(16,'(a,a)',advance="no") adjustl(trim(gf_method(j)))," "
      end do
      write(16,'(a,a)',advance="no") adjustl(trim(gf_basis))," "
      write(16,*) " "
      write(16,*)
      write(16,*) " RP-EVB reference point No. ",adjustl(trim(rp_names(i)))
      write(16,*)
      write(16,'(i4,i4)') charge,multi
      do j=1,natoms
         write(16,*) name(j),irc_coords(rp_points(i),j,:)*bohr
      end do
      write(16,*)
      close(16)
   end do
else if (software_gf .eq. "O") then
   sys_stat=chdir("rp_ref")
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! The folder rp_ref seems to be nonexistent!"
      call fatal
   end if
   do i=1,n_rp_pts
      inquire(file=adjustl(trim(rp_names(i))), exist=exist)
      if (.not. exist) call system("mkdir " // adjustl(trim(rp_names(i))))
      open(unit=16,file=adjustl(trim(rp_names(i))) // "/orca.inp",status="replace")
      if (min_procs .gt. 1) then
         if (min_procs .lt. 10) then
            write(16,'(a,i1)') " ! PAL",rp_procs
         else
            write(16,'(a,i2)') " ! PAL",rp_procs
         end if
      end if
      write(16,'(a)',advance="no") " ! freq engrad "
      do j=1,gf_met_words
         write(16,'(a,a)',advance="no") adjustl(trim(gf_method(j)))," "
      end do
      write(16,'(a,a)',advance="no") adjustl(trim(gf_basis))," "
      write(16,*) " "
      write(16,'(a,i4,i4)') "* xyz ", charge,multi
      do j=1,natoms
         write(16,*) name(j),irc_coords(rp_points(i),j,:)*bohr
      end do
      write(16,*) "*"
      close(16)
   end do
end if
write(*,*) " --> done!"
!
!     Start RP-EVB point calculations! If more processors are available than needed 
!     for a single calculation, start some of them in parallel.
!     Check if calculations are done, if its so, start new ones 
!

write(*,*) "Start QM reference calculations for all RP-EVB gradient+hessian points!"
write(*,*) " If enough processors are available, this might be done in parallel!"

ncalcs_round=nprocs_tot/rp_procs
if (mod(n_rp_pts,ncalcs_round) .eq. 0) then
   nrounds=n_rp_pts/ncalcs_round
else 
   nrounds=n_rp_pts/ncalcs_round+1
end if

write(*,'(a,i3)') " Number of needed reference calculations: ",n_rp_pts
write(*,'(a,i3)') " Total processors available: ",nprocs_tot
write(*,'(a,i2)') " Processors per calculation: ",rp_procs
write(*,'(a,i2)') " Calculations to be done simultaneously: ",ncalcs_round
write(*,'(a,i2)') " Rounds of calculations to be done: ",nrounds
!
!    Go into large loop for calculation rounds: First, start the calculations simultaneously
!    Then, wait until all started calculations are finished 
!
!    A) GAUSSIAN calculations
!
if (software_gf .eq. "G") then
   do i=1,nrounds
!
!    directly include possible "rest calculations"
!   
      if (i .lt. nrounds) then
         ncalcs_round=nprocs_tot/rp_procs
      else 
         if (mod(n_rp_pts,ncalcs_round) .ne. 0) then
            ncalcs_round=mod(n_rp_pts,ncalcs_round)
         else 
            ncalcs_round=nprocs_tot/rp_procs
         end if
      end if
!
!    Start the calculations 
!
   
      do j=1,ncalcs_round 
         sys_stat=chdir(adjustl(trim(rp_names((i-1)*nprocs_tot/rp_procs+j))))
         inquire(file="done", exist=exist)
         if (.not. exist) then
            sys_stat=system(trim(link_gaussian)// " gauss.com &")
            write(*,'(a,i3)') " Reference calculation stated for RP-EVB point ",&
                      & (i-1)*nprocs_tot/rp_procs+j
            if (sys_stat .ne. 0) then
               write(*,*) "ERROR! The Gaussian16 calculation in RP folder ",&
                      & trim(rp_names((i-1)*nprocs_tot/rp_procs+j))," failed!"
               call fatal
            end if 
         else 
            write(*,'(a,i3)') " Reference calculation already done for RP-EVB point ",&
                 & (i-1)*nprocs_tot/rp_procs+j
         end if
         sys_stat=chdir("..")
      end do
!
!     Check if all are finished 
!
      if (.not. allocated(round_done)) allocate(round_done(ncalcs_round))
      round_done=0
      do 
         do j=1,ncalcs_round 
            sys_stat=chdir(adjustl(trim(rp_names((i-1)*nprocs_tot/rp_procs+j))))
            sys_stat=system("tail -1 gauss.log > tail.out &")
            sys_stat=system("sleep 0.1")
            open (unit=16,file="tail.out")
            read(16,'(a)',iostat=readstat) a80
            if (readstat .eq. 0) then
               if (index(a80,'Normal termination of Gaussian') .ne. 0) then
                  if (round_done(j) .eq. 0) then
                     write(*,*) " Finished reference calculation for RP-EVB point ",&
                           & (i-1)*nprocs_tot/rp_procs+j
                     call system(" touch done")
                  end if
                  round_done(j)=1
                  write(*,*) "ERROR! The QMDFF reference calculation for RP-EVB point ",&
                           & (i-1)*nprocs_tot/rp_procs+j,"failed!"
                  call fatal
               end if
            end if
            close(16)
            sys_stat=chdir("..")
         end do
!
!     If all actual calculations are finished, all stata are 1 (sum=N)
!
         idum=sum(round_done)
         if (idum .eq. ncalcs_round) exit 
      end do 
      if (allocated(round_done)) deallocate(round_done)
   end do
!
!    B) ORCA calculations
!
else if (software_gf .eq. "O") then
   do i=1,nrounds
!
!    directly include possible "rest calculations"
!   
      if (i .lt. nrounds) then
         ncalcs_round=nprocs_tot/rp_procs
      else
         if (mod(n_rp_pts,ncalcs_round) .ne. 0) then
            ncalcs_round=mod(n_rp_pts,ncalcs_round)
         else
            ncalcs_round=nprocs_tot/rp_procs
         end if
      end if
!
!    Start the calculations 
!
   
      do j=1,ncalcs_round
         sys_stat=chdir(adjustl(trim(rp_names((i-1)*nprocs_tot/rp_procs+j))))
         inquire(file="done", exist=exist)
         if (.not. exist) then
            sys_stat=system(trim(link_orca)// " orca.inp > orca.log &")
            write(*,'(a,i3)') " Reference calculation stated for RP-EVB point ",&
                      & (i-1)*nprocs_tot/rp_procs+j
            if (sys_stat .ne. 0) then
               write(*,*) "ERROR! The orca5 calculation in RP folder ",&
                      & trim(rp_names((i-1)*nprocs_tot/rp_procs+j))," failed!"
               call fatal
            end if
         else
            write(*,'(a,i3)') " Reference calculation already done for RP-EVB point ",&
                 & (i-1)*nprocs_tot/rp_procs+j
         end if
         sys_stat=chdir("..")
      end do
!
!     Check if all are finished 
!
      if (.not. allocated(round_done)) allocate(round_done(ncalcs_round))
      round_done=0
      do
         do j=1,ncalcs_round
            sys_stat=chdir(adjustl(trim(rp_names((i-1)*nprocs_tot/rp_procs+j))))
            sys_stat=system("tail -2 orca.log > tail.out &")
            sys_stat=system("sleep 0.1")
            open (unit=16,file="tail.out")
            read(16,'(a)',iostat=readstat) a80
            if (readstat .eq. 0) then
               if (index(a80,'****ORCA TERMINATED NORMALLY****') .ne. 0) then
                  if (round_done(j) .eq. 0) then
                     write(*,*) " Finished reference calculation for RP-EVB point ",&
                           & (i-1)*nprocs_tot/rp_procs+j
                     call system(" touch done")
                  end if
                  round_done(j)=1
               else if ((index(a80,'ERROR') .ne. 0) .or. (index(a80,'error') .ne. 0) .or. &
                     &   (index(a80,'Aborting the run') .ne. 0)) then

                  write(*,*) "ERROR! The QMDFF reference calculation for RP-EVB point ",&
                           & (i-1)*nprocs_tot/rp_procs+j,"failed!"
                  call fatal
               end if
            end if
            close(16)
            sys_stat=chdir("..")
         end do
!
!     If all actual calculations are finished, all stata are 1 (sum=N)
!
         idum=sum(round_done)
         if (idum .eq. ncalcs_round) exit
      end do
      if (allocated(round_done)) deallocate(round_done)
   end do

end if

call system_clock(time_int(5))


!
!     TASK 4-B : CALCULATE SEPARATE TREQ ENERGY REFERENCE
!
!     If a higher level of theory shall be used for energy reference 
!     along the path, generate the needed input data for these 
!     calculations in this folder 
!
if (ens_extra) then
   call chdir("..")
   inquire(file="add_ens", exist=exist)
   if (.not. exist) call system("mkdir add_ens")
   call chdir("add_ens")
 
   write(*,*)
   write(*,*) "---- TASK 4-B: Calculate separate energy data for TREQ  ----------"
   write(*,*)
   ! irc_num
!
!    Determine how many calculations might be done in parallel: 
!    multiplicity of ens_procs that fits into total number of processors
!
   if (2*ens_procs .le. nprocs_tot) then
      nrounds=nprocs_tot/ens_procs
   else 
      nrounds=1
   end if
   njobs_done=1
   write(*,*) "The additional energy calculations will be done in parallel with"
   write(*,'(a,i2,a)') " ",nrounds," calculations, each with a part of the IRC."
   do i=1,nrounds 
      if (i .le. 9) then
         write(adum,'(a,i1)') "round",i
      else 
         write(adum,'(a,i2)') "round",i
      end if
!
!     Determine range of structures to be calculated in this round
!
      if (mod(irc_num,nrounds) .ge. i) then
         njobs_round=irc_num/nrounds+1
      else 
         njobs_round=irc_num/nrounds
      end if
      inquire(file=trim(adum), exist=exist)
      if (.not. exist) call system("mkdir "//trim(adum))
      call chdir(trim(adum))
      write(*,'(a,i2,a,i5,a,i5,a)') " Write input for calculation part ",i, &
                 " (structures ",njobs_done," to ",njobs_done+njobs_round-1,")"
!
!     Write orca input file for calculation
!
      open(unit=15,file="orca.inp",status="unknown") 
      write(15,'(a)',advance="no") "! "
      do j=1,e_met_words
         write(15,'(a,a)',advance="no") adjustl(trim(e_method(j)))," "
      end do
      write(15,'(a)') trim(e_basis)
      if (ens_procs .gt. 1) then
         write(15,'(a,i1)') "! PAL",ens_procs
      end if
      if (min_mem .le. 9) then
         write(15,'(a,i1,a)') "%maxcore ",min_mem,"000"
      else 
         write(15,'(a,i2,a)') "%maxcore ",min_mem,"000"
      end if
      write(15,'(a,i3,a,i3,a)') "*xyzfile ",charge," ",multi," traj.xyz"
      close(15)
!
!     Write external trajectory file for calculation
! 
      open(unit=16,file="traj.xyz",status="unknown")
      do j=njobs_done,njobs_done+njobs_round-1
         if (j .ne. njobs_done) then
            write(16,'(a)') ">"
         end if
         write(16,*) natoms 
         write(16,*) 
         do k=1,natoms
            write(16,*) name(k),irc_coords(j,k,:)*bohr
         end do
      end do
      write(16,*) 
      close(16)
!
!     Go into next loop and increment upper border...
!
      call chdir("..")
      njobs_done=njobs_done+njobs_round

      write(*,*) " --> done..."
   end do   
!
!    Start the calculations 
!

   do i=1,nrounds
      if (i .le. 9) then
         write(adum,'(a,i1)') "round",i
      else
         write(adum,'(a,i2)') "round",i
      end if
      sys_stat=chdir(adjustl(trim(adum)))
      inquire(file="done", exist=exist)
      if (.not. exist) then
         sys_stat=system(trim(link_orca)// " orca.inp > orca.out &")
         write(*,'(a,i3)') " Extra energy reference calculation started for part ",i
         if (sys_stat .ne. 0) then
            write(*,*) "ERROR! The orca calculation in folder ",trim(adum)," failed!"
            call fatal
         end if
      else
         write(*,'(a,i3)') " Extra energy reference calculation already done for part ",i
      end if
      sys_stat=chdir("..")
   end do

!
!     Check if all are finished 
!
   if (allocated(round_done)) deallocate(round_done)
   if (.not. allocated(round_done)) allocate(round_done(nrounds))
   round_done=0
   do
      do j=1,nrounds
         if (j .le. 9) then
            write(adum,'(a,i1)') "round",j
         else
            write(adum,'(a,i2)') "round",j
         end if
         sys_stat=chdir(adjustl(trim(adum)))
         sys_stat=system("tail -1 orca.out > tail.out &")
         sys_stat=system("sleep 0.1")
         open (unit=16,file="tail.out")
         read(16,'(a)',iostat=readstat) a80
         if (readstat .eq. 0) then
            if (index(a80,'TOTAL RUN TIME: ') .ne. 0) then
               if (round_done(j) .eq. 0) then
                  write(*,'(a,i3)') " Finished extra energy calculation for part ",j
                  call system(" touch done")
               end if
               round_done(j)=1
            else if (index(a80,'File lengths (MBytes)') .ne. 0) then
               write(*,*) "ERROR! The extra energy calculation for part ",j,"failed!"
               call fatal
            end if
         end if
         close(16)
         sys_stat=chdir("..")
      end do
!
!     If all actual calculations are finished, all stata are 1 (sum=N)
!
      idum=sum(round_done)
      if (idum .eq. nrounds) exit
   end do
   if (allocated(round_done)) deallocate(round_done)

   call chdir("..")
   call chdir("rp_ref")
end if
!
!     TASK 5-B : READ IN EXTRA ENERGY REFERENCES
!
if (ens_extra) then
   sys_stat=chdir("..")
   write(*,*)
   write(*,*) "---- TASK 5-B: Read in the extra energy data for the IRC  ----------"
   write(*,*)
 
   sys_stat=chdir('add_ens')
   njobs_done=1

   do i=1,nrounds
      if (i .le. 9) then
         write(adum,'(a,i1)') "round",i
      else
         write(adum,'(a,i2)') "round",i
      end if
!
!     Determine range of structures to be calculated in this round
!
      if (mod(irc_num,nrounds) .ge. i) then
         njobs_round=irc_num/nrounds+1
      else
         njobs_round=irc_num/nrounds
      end if
!
!     If only one job will be done in this round, read in directly from orca.out file,
!     else, read in from the orca.xyzact file 
!
      call chdir(trim(adum))
      write(*,'(a,i2,a,i5,a,i5,a)') " Read in extra reference energies for part ",i, &
                " (structures ",njobs_done," to ",njobs_done+njobs_round-1,")"
      call system("pwd")
      if (njobs_round .gt. 1) then
         open(unit=15,file="orca.xyzact.dat",status="old",iostat=readstat)
         if (readstat .ne. 0) then
            write(*,*) "ERROR! The file with the resulting energies (orca.xyzact.dat) is not there!"
            call fatal
         end if
! 
!     Read in the energies: DIRECTLY overwrite the irc_energies!
!
         do j=njobs_done,njobs_done+njobs_round-1
            read(15,*) rdum,irc_energies(j)
         end do
         close(15)
      else 
         open(unit=15,file="orca.out")
17       read(15,'(a)',end=29)a80
         if (index(a80,'FINAL SINGLE POINT ENERGY ') .ne. 0) then
            read(a80,*) adum,adum,adum,adum,irc_energies(njobs_done)
         end if
         goto 17
29       continue
         close(15,status='keep')
      end if
!
!     Go into next loop and increment upper border...
!
      call chdir("..")
      njobs_done=njobs_done+njobs_round

      write(*,*) " --> done..."
   end do
!
!     Rewrite the file(s) for the IRC energies!
!
   call chdir('..')
   open(unit=17,file="mep_irc/irc_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) irc_energies(i)
   end do
   close(17)

   open(unit=17,file="mep_irc/irc_xi_ens.dat",status="replace")
   do i=1,irc_num
      write(17,*) xi_vals(i),irc_energies(i)
   end do
   close(17)

   sys_stat=chdir("rp_ref")
end if

!
!     TASK 5 : READ IN RP-EVB GRADIENT+HESSIAN REFERENCE
!
write(*,*)
write(*,*) "---- TASK 5: Read in the reference data for RP-EVB  ----------"
write(*,*)
!
!     Read in the results of the RP-EVB reference calculations!
!     They will be given for the EVB-QMDFF optimization in compressed formate 
!     Write simultaneously the big reference file "grad_hess.dat"
!   
nat3=3*natoms
allocate(grad(3,natoms))
allocate(grad1d(3*natoms))
allocate(hess(nat3,nat3))
allocate(hess1d(((3*natoms)*(3*natoms+1))/2))
allocate(hess1d2((3*natoms)*(3*natoms)))
allocate(g_out(nat3),xyz_out(nat3),h_out(nat3*nat3))
dg_ref_file="grad_hess.dat"
open(unit=15,file=dg_ref_file,status="unknown")
write(15,'(A)') "# This is an input file for TREQ calculations, "
write(15,'(A)') "# generated by black_box.x."
write(15,*)
write(15,'((A),(I3))') "NPOINTS",n_rp_pts
write(15,*)
write(15,'((A),(I4))') "NATOM",natoms
write(15,*)
!
!     Generate file header and headers for single point sections 
!
do i=1,n_rp_pts
   sys_stat=chdir(adjustl(trim(rp_names(i))))
   if (sys_stat .ne. 0) then
      write(*,*) "ERROR! Folder ",adjustl(trim(rp_names(i))),"seems to be nonexistent!" 
      call fatal
   end if
   write(15,'((a),(i3))') "*POINT",i
   write(*,'((a),(i3),(a))') " Read in RP-EVB point No.",i," ... "
   write(15,*)
!
!      Now read in reference data 
!  
!      A: for orca calculations 
!
   if (software_gf .eq. "O") then
      open(unit=11,file="orca.engrad")

!
!     Read in actual energy
!
10    read(11,'(a)',end=30)a80
      if(index(a80,'# The current total energy in Eh').ne.0) then
         read(11,'(a)',end=30)a80
         read(11,*) energy
      end if
!
!     Read in the gradient 
!
      if(index(a80,'# The current gradient in Eh/bohr').ne.0)then
         read(11,'(a)',end=30)a80
         do k=1,natoms
            do m=1,3
               read(11,*) grad(m,k)
            enddo
         enddo
      endif
      goto 10
30    close(11)
!
!     Read in the hessian
!
      open(unit=13,file="orca.hess")
      m=nat3/5
      if(mod(nat3,5).gt.0)m=m+1

11    read(13,'(a)',end=31)a80
      if(index(a80,'$act_energy').ne.0)snforca=.false.
      if(index(a80,'$hessian').ne.0)then
         read(13,'(a)',end=31)a80
         maxcol = 0
         do k=1,m
            read(13,'(a)')a80
            mincol = maxcol + 1
            maxcol = min(maxcol+5,nat3)
            if(snforca)then
               do l=1,nat3
                  read(13,*)adum,(hess(l,j),j=mincol,maxcol)
               enddo
            else
               do l=1,nat3
                  read(13,*)idum,(hess(l,j),j=mincol,maxcol)
               enddo
            endif
         enddo
      endif
      goto 11
31    close(13,status='keep')
!
!     B: for Gaussian calculations
!
   else if (software_gf .eq. "G") then
!
!    First execute the formchk command in order to get the formatted output
! 
      formchk_stat=system("formchk gauss.chk gauss.fchk > formchk.log")
      if (formchk_stat .ne. 0) then
         write(*,*) "ERROR! The formchk command was not successful!"
         call fatal
      end if
      open(unit=11,file="gauss.fchk")
!
!   Read in the energy
!
12    read(11,'(a)',end=32)a80
      if(index(a80,'Total Energy').ne.0) then
         read(a80,*) adum,adum,adum,energy
      end if
      maxcol=0

! 
!   Read in gradient (1d Array)
!
      if(index(a80,'Cartesian Gradient').ne.0)then
!
!   if all lines are full, read one less
!
         if (int(nat3/5)*5 .eq. int(nat3)) then
            gradline=int(nat3/5)
         else
            gradline=int(nat3/5)+1
         end if

         do l=1,gradline
            mincol = maxcol + 1
            maxcol = min(maxcol+5,nat3)

            read(11,*) (grad1d(j),j=mincol,maxcol)
         enddo
         do l=1,natoms
            do j=1,3
               grad(j,l)=grad1d((l-1)*3+j)
            end do
         end do
      end if
!
!   Read in the hessian (upper triangular matrix)
!
      
      if(index(a80,'Cartesian Force Constants').ne.0)then
         do k=1,int(((nat3*(nat3+1))/10)+1)
            mincol = maxcol + 1
            maxcol = min(maxcol+5,(nat3*(nat3+1))/2)
            read(11,*) (hess1d(j),j=mincol,maxcol)
         enddo
         do k=1,nat3
            do j=1,nat3
               if (j .le. k) then
                  natsum=((k-1)*k)/2
                  sum1=(k-1)*nat3+j
                  sum2=(j-1)*nat3+k
                  hess1d2(sum1)=hess1d(natsum+j)
                  hess1d2(sum2)=hess1d2(sum1)
               end if
            end do
         end do
      end if

      goto 12
      stop
32    close(11)

!
!     Write better two dimensional hessian array
!

      do k=1,nat3
         do j=1,nat3
            hess(k,j)=hess1d2((k-1)*nat3+j)
         end do
      end do
   end if

!
!     Now write reference informations for actual RP-EVB point to file 
!     grad_hess.dat !
!
   write(*,*) "Write reference informations to file grad_hess.dat ..."
!
!     1) The energy (either from the extra theory level or the actual level)
!
   if (ens_extra) then
      write(15,'((a),(f18.10))') "ENERGY",irc_energies(rp_points(i))
   else 
      write(15,'((a),(f18.10))') "ENERGY",energy
   end if
   write(15,*)
!
!     2) The cartesian coordinates of the structure 
!
   write(15,'(a)') "GEOMETRY"
   do j=1,natoms
      do m=1,3
         xyz_out((j-1)*3+m)=irc_coords(rp_points(i),j,m)
      end do
   end do
   maxline= int(size(xyz_out)/5)
   do j=0,maxline-1
      write(15,'((e16.7),(e16.7),(e16.7),(e16.7),(e16.7))') &
          &xyz_out(j*5+1),xyz_out(j*5+2),xyz_out(j*5+3),xyz_out(j*5+4),xyz_out(j*5+5)
   end do
   if (size(xyz_out)-maxline*5 .eq.4) then
       write(15,'((e16.7),(e16.7),(e16.7),(e16.7))') &
          &xyz_out(j*5+1),xyz_out(j*5+2),xyz_out(j*5+3),xyz_out(j*5+4)
   else if (size(xyz_out)-maxline*5 .eq.3) then
       write(15,'((e16.7),(e16.7),(e16.7))') &
          &xyz_out(j*5+1),xyz_out(j*5+2),xyz_out(j*5+3)
   else if (size(xyz_out)-maxline*5 .eq.2) then
       write(15,'((e16.7),(e16.7))') &
          &xyz_out(j*5+1),xyz_out(j*5+2)
   else if (size(xyz_out)-maxline*5 .eq.1) then
       write(15,'((e16.7))') &
          &xyz_out(j*5+1)
   end if
   write(15,'(A)') "END"
   write(15,*)
!
!     3) The cartesian gradient of the structure 
!
   write(15,'(A)') "GRADIENT"
   do j=1,natoms
      do m=1,3
         g_out((j-1)*3+m)=grad(m,j)
      end do
   end do
   do j=0,maxline-1
      write(15,'((e16.7),(e16.7),(e16.7),(e16.7),(e16.7))') &
          &g_out(j*5+1),g_out(j*5+2),g_out(j*5+3),g_out(j*5+4),g_out(j*5+5)
   end do
   if (size(xyz_out)-maxline*5 .eq.4) then
      write(15,'((e16.7),(e16.7),(e16.7),(e16.7))') &
          &g_out(j*5+1),g_out(j*5+2),g_out(j*5+3),g_out(j*5+4)
   else if (size(xyz_out)-maxline*5 .eq.3) then
      write(15,'((e16.7),(e16.7),(e16.7))') &
          &g_out(j*5+1),g_out(j*5+2),g_out(j*5+3)
   else if (size(xyz_out)-maxline*5 .eq.2) then
      write(15,'((e16.7),(e16.7))') &
          &g_out(j*5+1),g_out(j*5+2)
   else if (size(xyz_out)-maxline*5 .eq.1) then
      write(15,'((e16.7))') &
          &g_out(j*5+1)
   end if
!
!     4) The cartesian hessian of the structure 
!
   write(15,'(A)') "END"
   write(15,*)
   write(15,'(A)') "HESSIAN"
   do j=1,3*natoms
      do m=1,3*natoms
         h_out((j-1)*3*natoms+m)=hess(m,j)
      end do
   end do
   maxline= int(size(h_out)/5)
   do j=0,maxline-1
      write(15,'((e16.7),(e16.7),(e16.7),(e16.7),(e16.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3),h_out(j*5+4),h_out(j*5+5)
   end do
   if (size(h_out)-maxline*5 .eq.4) then
      write(15,'((e16.7),(e16.7),(e16.7),(e16.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3),h_out(j*5+4)
   else if (size(h_out)-maxline*5 .eq.3) then
      write(15,'((e16.7),(e16.7),(e16.7))') &
          &h_out(j*5+1),h_out(j*5+2),h_out(j*5+3)
   else if (size(h_out)-maxline*5 .eq.2) then
      write(15,'((e16.7),(e16.7))') &
          &h_out(j*5+1),h_out(j*5+2)
   else if (size(h_out)-maxline*5 .eq.1) then
      write(15,'((e16.7))') &
          &h_out(j*5+1)
   end if

   write(15,'(a)') "END"
   write(15,*)




   write(*,*) " --> done..."
   sys_stat=chdir("..")
end do




call system_clock(time_int(6))

!
!     TASK 6 : Generate input folder for dynamic calculations
!
write(*,*)
write(*,*) "---- TASK 6: Generate input folder for dynamic calculations  ----------"
write(*,*)
sys_stat=chdir("..")
inquire(file="rpmd", exist=exist)
if (.not. exist) sys_stat=system("mkdir calc_rate")

!
!     First, read in the QMDFFs, from their arrays informations about 
!     fragments and bonds will be extracted 
!
sys_stat=chdir("evb_qmdff")
call prepare ("reactants.qmdff","products.qmdff","dummy",2)
! 
!     For later usage: read in the qmdff energies!
!
open(unit=15,file="caracal.key",status="old")
do 
   read(15,'(a)',iostat=lastline) a80
   if (lastline .ne. 0) exit
   if (index(a80,' ESHIFT ') .ne. 0) then
      qmdff_en_line=a80
   end if
end do
close(15)
!
!     Read in fragments from qmdff logfile!
!
allocate(fragments(natoms))
open(unit=15,file="reactants_qmdff.log",status="old")
do
   read(15,'(a)',iostat=lastline) a80

   if (lastline .ne. 0) exit
!
!     Read in fragment definition fragment 
!
   if (index(a80,'#        belongs to fragment') .ne. 0) then
      read(15,'(a)') a80
      do i=1,natoms 
         read(15,*) idum,fragments(i)     
      end do
      exit 
   end if
end do

close(15)

!
!     Determine number of fragments and belonging atoms for the reactant
!     structure 
!

allocate(fragsize(4))
fragsize=0
do i=1,natoms
   if (fragments(i) .gt. 4) then
      write(*,*) "ERROR! More than 4 fragments found in reactant structure!"
      call fatal
   end if 
   do j=1,4
      if (fragments(i) .eq. j) then
         fragsize(j)=fragsize(j)+1
      end if 
   end do
end do
!
!     Allocate atom index arrays for all existent fragments
!
if (fragsize(1) .gt. 0) then
   allocate(frag1(fragsize(1)))
end if
if (fragsize(2) .gt. 0) then
   allocate(frag2(fragsize(2)))
end if
if (fragsize(3) .gt. 0) then
   allocate(frag3(fragsize(3)))
end if
if (fragsize(4) .gt. 0) then
   allocate(frag4(fragsize(4)))
end if
!
!     Determine members of fragments 
!

allocate(atfrags(4))
atfrags=0
do i=1,natoms 
   if (fragments(i) .eq. 1) then 
      atfrags(1)=atfrags(1)+1
      frag1(atfrags(1))=i
   else if (fragments(i) .eq. 2) then
      atfrags(2)=atfrags(2)+1
      frag2(atfrags(2))=i
   else if (fragments(i) .eq. 3) then
      atfrags(3)=atfrags(3)+1
      frag3(atfrags(3))=i
   else if (fragments(i) .eq. 4) then
      atfrags(4)=atfrags(4)+1
      frag4(atfrags(4))=i 
   end if
end do

!
!     Determine number of molecules in reactant system
!

if (atfrags(4) .ne. 0)  then
   num_eds=4
else if (atfrags(3) .ne. 0)  then
   num_eds=3
else if (atfrags(2) .ne. 0)  then
   num_eds=2
else 
   num_eds=1
end if
!
!     Write info message about atoms in reactant structures
!
write(*,'(a,i1)') " Number of molecules in reactant system: ",num_eds
write(*,'(a,40i4)') " Atoms of first molecule: ",frag1
if (num_eds .ge. 2) then
   write(*,'(a,40i4)') " Atoms of second molecule: ",frag2
else if (num_eds .ge. 3) then
   write(*,'(a,40i4)') " Atoms of third molecule: ",frag3
else 
   write(*,'(a,40i4)') " Atoms of fourth molecule: ",frag4
end if

!
!     Determine which bonds are the forming or breaking bonds during 
!     the reaction: Determine bond pools of both reactants and products
!     and determine differences between them
!
allocate(atind(natoms))
!
!     Determine element numbers of the systems atoms
!
do i=1,natoms
   atind(i)=at(i)
end do
!
!     COVALENT RADII of all elements
!     based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!     in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!     edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!     corrected Nov. 17, 2010 for the 92nd edition.
!
rad = (/ &
 & 0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0, &
 & 0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0, &
 & 2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0, &
 & 1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0, &
 & 2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0, &
 & 1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0, &
 & 2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0, &
 & 1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0, &
 & 1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0, &
 & 1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0, &
 & 1.84D0,1.83D0,1.80D0,1.80D0 /)
rad=rad/bohr
!
!     set tolerance factor for bondlengths (if a bond shall still be defined..)
!
bondlentol=1.2d0
!
!     Loop over all combinations of bonds that are short enough to be 
!     covalent. Store all bonds that are present in both structures to
!     array with double bonds
!
allocate(doubles(2,max(nbond,nbond_two)))
n_double=0
do i=1,nbond
   if (vbond(1,i) .lt. (rad(atind(bond(1,i)))+rad(atind(bond(2,i))))*bondlentol) then
      do j=1,nbond_two
         if (vbond_two(1,j) .lt. (rad(atind(bond_two(1,j)))+rad(atind(bond_two(2,j))))*&
                    & bondlentol) then
            if (((bond_two(1,j) .eq. bond(1,i)) .and. (bond_two(2,j) .eq. bond(2,i)) .or. &
              &  ((bond_two(1,j) .eq. bond(2,i)) .and. (bond_two(2,j) .eq. bond(1,i))))) then
           !    write(*,*) "both!",bond(:,i),bond_two(:,j)
               n_double=n_double+1
               doubles(:,n_double)=bond(:,i)
            else 
            end if
         end if
      end do
   end if
end do
!
!     Determine the number and indices of breaking bonds: All bonds that 
!     are in reactants but not in the double array
!
allocate(breaks_tmp(2,max(nbond,nbond_two)))
n_break=0
do i=1,nbond
   idum=0   
   if (vbond(1,i) .lt. (rad(atind(bond(1,i)))+rad(atind(bond(2,i))))*bondlentol) then
      do j=1,n_double
         if (((doubles(1,j) .eq. bond(1,i)) .and. (doubles(2,j) .eq. bond(2,i)) .or. &
            &  ((doubles(1,j) .eq. bond(2,i)) .and. (doubles(2,j) .eq. bond(1,i))))) then
            idum=1 
         end if

      end do
   else 
      idum=1
   end if
   if (idum .eq. 0) then
      n_break=n_break+1
      breaks_tmp(:,n_break)=bond(:,i)
   end if
end do
allocate(breaks(2,n_break))
allocate(print_break(n_break))
do i=1,n_break
   breaks(:,i)=breaks_tmp(:,i)
   if (breaks(1,i) .lt. 10) then
      if (breaks(2,i) .lt. 10) then
         write(print_break(i),'(i1,a,i1)') breaks(1,i),"-",breaks(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_break(i),'(i1,a,i2)') breaks(1,i),"-",breaks(2,i)
      else
         write(print_break(i),'(i1,a,i3)') breaks(1,i),"-",breaks(2,i)
      end if
   else if (breaks(1,i) .lt. 100) then
      if (breaks(2,i) .lt. 10) then
         write(print_break(i),'(i2,a,i1)') breaks(1,i),"-",breaks(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_break(i),'(i2,a,i2)') breaks(1,i),"-",breaks(2,i)
      else 
         write(print_break(i),'(i2,a,i3)') breaks(1,i),"-",breaks(2,i)
      end if
   else 
      if (breaks(2,i) .lt. 10) then
         write(print_break(i),'(i3,a,i1)') breaks(1,i),"-",breaks(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_break(i),'(i3,a,i2)') breaks(1,i),"-",breaks(2,i)
      else 
         write(print_break(i),'(i3,a,i3)') breaks(1,i),"-",breaks(2,i)
      end if
   end if
end do
deallocate(breaks_tmp)
!
!     print info about breaking bonds 
!
write(*,'(a)',advance="no") " List of breaking bonds: "
do i=1,n_break
   write(*,'(a,a)',advance="no") " ",trim(print_break(i))
end do
write(*,'(a)') " "

!
!     Determine the number and indices of forming bonds: All bonds that 
!     are in products but not in the double array
!
allocate(forms_tmp(2,max(nbond,nbond_two)))
n_form=0
do i=1,nbond_two
   idum=0
   if (vbond_two(1,i) .lt. (rad(atind(bond_two(1,i)))+rad(atind(bond_two(2,i))))*bondlentol) then
      do j=1,n_double
         if (((doubles(1,j) .eq. bond_two(1,i)) .and. (doubles(2,j) .eq. bond_two(2,i)) .or. &
            &  ((doubles(1,j) .eq. bond_two(2,i)) .and. (doubles(2,j) .eq. bond_two(1,i))))) then
            idum=1
         end if

      end do
   else
      idum=1
   end if
   if (idum .eq. 0) then
      n_form=n_form+1
      forms_tmp(:,n_form)=bond_two(:,i)
   end if
end do
allocate(forms(2,n_form))
allocate(print_form(n_form))
do i=1,n_form
   forms(:,i)=forms_tmp(:,i)
   if (forms(1,i) .lt. 10) then
      if (forms(2,i) .lt. 10) then
         write(print_form(i),'(i1,a,i1)') forms(1,i),"-",forms(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_form(i),'(i1,a,i2)') forms(1,i),"-",forms(2,i)
      else
         write(print_form(i),'(i1,a,i3)') forms(1,i),"-",forms(2,i)
      end if
   else if (forms(1,i) .lt. 100) then
      if (forms(2,i) .lt. 10) then
         write(print_form(i),'(i2,a,i1)') forms(1,i),"-",forms(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_form(i),'(i2,a,i2)') forms(1,i),"-",forms(2,i)
      else
         write(print_form(i),'(i2,a,i3)') forms(1,i),"-",forms(2,i)
      end if
   else
      if (forms(2,i) .lt. 10) then
         write(print_form(i),'(i3,a,i1)') forms(1,i),"-",forms(2,i)
      else if (breaks(2,i) .lt. 100) then
         write(print_form(i),'(i3,a,i2)') forms(1,i),"-",forms(2,i)
      else
         write(print_form(i),'(i3,a,i3)') forms(1,i),"-",forms(2,i)
      end if
   end if
end do
deallocate(forms_tmp)

write(*,'(a)',advance="no") " List of forming bonds: "
do i=1,n_form
   write(*,'(a,a)',advance="no") " ",trim(print_form(i))
end do
write(*,'(a)') " "

!
!     Determine which type of reaction the actual path describes!
!     --> This defines the type of dividing surfaces used in the samplings 
!     It will be determined from the number of reactants and the number 
!     of forming and breaking bonds
!

if ((n_form .eq. 1) .and. (n_break .eq. 1) .and. (num_eds .eq. 1)) then
   write(*,*) "There are one breaking and one forming bonds as well as one "
   write(*,*) " reactant molecule! Therefore the reaction type is REARRANGE!"
   umbr_type="REARRANGE"
else if ((n_form .eq. 1) .and. (n_break .eq. 1) .and. (num_eds .eq. 2)) then
   write(*,*) "There are one breaking and one forming bonds as well as two "
   write(*,*) " reactant molecules! Therefore the reaction type is BIMOLECULAR!"
   umbr_type="BIMOLEC"
else if ((n_form .eq. 2) .and. (n_break .eq. 0) .and. (num_eds .eq. 2)) then
   write(*,*) "There are two breaking and no forming bonds as well as two "
   write(*,*) " reactant molecules! Therefore the reaction is a CYCLOADDIOTION!"
   umbr_type="CYCLOADD"
else if ((n_form .eq. 1) .and. (n_break .eq. 0) .and. (num_eds .eq. 2)) then
   write(*,*) "There are no breaking and one forming bonds as well as two "
   write(*,*) " reactant molecules! Therefore the reaction is an MERGING!"
   umbr_type="MERGING"
else if ((n_form .eq. 2) .and. (n_break .eq. 1) .and. (num_eds .eq. 2)) then
   write(*,*) "There are one breaking and two forming bonds as well as two "
   write(*,*) " reactant molecules! Therefore the reaction is an ADDITION!"
   umbr_type="ADDITION"
else if ((n_form .eq. 3) .and. (n_break .eq. 2) .and. (num_eds .eq. 3)) then
   write(*,*) "There are three breaking and two forming bonds as well as three "
   write(*,*) " reactant molecules! Therefore the reaction is a cyclic addition with"
   write(*,*) " three molecules (ADDITION3)!"
   umbr_type="ADDITION3"
else if ((n_form .eq. 4) .and. (n_break .eq. 3) .and. (num_eds .eq. 4)) then
   write(*,*) "There are four breaking and three forming bonds as well as four "
   write(*,*) " reactant molecules! Therefore the reaction is a cyclic addition with"
   write(*,*) " four molecules (ADDITION4)!"
   umbr_type="ADDITION4"
else if ((n_form .eq. 0) .and. (n_break .eq. 2) .and. (num_eds .eq. 1)) then
   write(*,*) "There are two breaking and no forming bonds for one reactant "
   write(*,*) " molecule. The reaction type is a CYCLOREVER(SION)!"
   umbr_type="CYCLOREVER"
else if ((n_form .eq. 0) .and. (n_break .eq. 1) .and. (num_eds .eq. 1)) then
   write(*,*) "There are one breaking and no forming bonds for one reactant "
   write(*,*) " molecule. The reaction type is a 1bond-decomposition"
   write(*,*) " (DECOM_1BOND)"
   umbr_type="DECOM_1BOND"
else if ((n_form .eq. 1) .and. (n_break .eq. 2) .and. (num_eds .eq. 1)) then
   write(*,*) "There are two breaking and one forming bonds for one reactant "
   write(*,*) " molecule. The reaction type is an usual elimination"
   write(*,*) " (ELIMINATION)"
   umbr_type="ELIMINATION"
else 
   write(*,*) "ERROR! No reaction type with ",n_form," forming, ",n_break," breaking"
   write(*,*) " bonds and ",num_eds," reactant molecules is supported so far!"
   call fatal
end if

sys_stat=chdir("..")
sys_stat=chdir("calc_rate")
!
!     Generate input folders for RPMD calculations and generate needed 
!     input files (copy structure/energy info and write caracal.key file)
!

do i=1,temp_num
   do k=1,kt_avg 
      if (temps(i) .lt. 10) then
         write(adum,'(i1)') temps(i)
      else if (temps(i) .lt. 100) then
         write(adum,'(i2)') temps(i)
      else if (temps(i) .lt. 1000) then
         write(adum,'(i3)') temps(i)
      else 
         write(adum,'(i4)') temps(i)
      end if
!
!     Generate the new folder if not already existent
!
      if (k .eq. 1) then
         write(adum,'(a,a)') adjustl(trim(adum)),"K"
         inquire(file=adjustl(trim(adum)), exist=exist)
         if (.not. exist) sys_stat=system("mkdir " // adjustl(trim(adum))) 
         call chdir(adjustl(trim(adum)))
      end if
!
!     If more than one run shall be started, generate subfolders 
!
      if (kt_avg .gt. 1) then
         if (k .gt. 9) then
            write(adum2,'(a,i2)') "run",k
         else 
            write(adum2,'(a,i1)') "run",k
         end if
         inquire(file=adjustl(trim(adum2)), exist=exist)
         if (.not. exist) sys_stat=system("mkdir " // adjustl(trim(adum2)))
         call chdir(adjustl(trim(adum2)))
      end if
!
!     Copy input data to file 
!
      if (kt_avg .eq. 1) then
         call system("cp ../../rp_ref/grad_hess.dat .")
         call system("cp ../../evb_qmdff/reactants.qmdff .")
         call system("cp ../../evb_qmdff/products.qmdff .")
         call system("cp ../../mep_irc/irc.xyz .")
         call system("cp ../../mep_irc/irc_ens.dat .")
         if (manual_ints) then
            call system("cp ../../" // trim(coord_file) // " coord_def.inp")
         end if
      else 
         call system("cp ../../../rp_ref/grad_hess.dat .")
         call system("cp ../../../evb_qmdff/reactants.qmdff .")
         call system("cp ../../../evb_qmdff/products.qmdff .")
         call system("cp ../../../mep_irc/irc.xyz .")
         call system("cp ../../../mep_irc/irc_ens.dat .")
         if (manual_ints) then
            call system("cp ../../../" // trim(coord_file) // " coord_def.inp")
         end if
      end if
!
!     Write new file for TS starting structure 
!    
      open(unit=15,file="ts.xyz",status="unknown")
      write(15,*) natoms 
      write(15,*)
      do j=1,natoms 
         write(15,*) name(j),irc_coords(pos_ts1,j,:)*bohr
      end do
      close(15) 
!
!     For unimolecular reactions: Write reactants structure to extra file 
!
      open(unit=15,file="reactants.xyz",status="unknown")
      write(15,*) natoms
      write(15,*)
      do j=1,natoms
         write(15,*) name(j),irc_coords(1,j,:)*bohr
      end do
      close(15)

!
!     Write the caracal.key file 
!
      open(unit=15,file="calc_rate.key",status="unknown")
      write(15,'(a)') "####################################"
      write(15,'(a)') "#  This is the keyfile for a k(T)  #"
      write(15,'(a)') "#  calculation with calc_rate.x    #"
      write(15,'(a)') "#   Generated by black_box.x       #"
      write(15,'(a)') "####################################"
      write(15,*)    
      write(15,'(a)') "##################################"
      write(15,'(a)') "# general settings"
      write(15,'(a)') "##################################"
      write(15,'(a)') "# xyz file with start structure of TS:"
      write(15,'(a)') "    ts_struc ts.xyz"
      write(15,'(a)') "# number of ring polymer beads in the system:"
      write(15,'(a,i4)') "    rpmd_beads ",nbeads 
      write(15,'(a)') "##################################"
      write(15,'(a)') "# potential energy surface settings"
      write(15,'(a)') "##################################"
      write(15,'(a)') "# names of the QMDFF files 1 and 2:"
      write(15,'(a)') "    qmdffnames reactants.qmdff products.qmdff"
      write(15,'(a)') "# QMDFF energy shifts for both minima:"
      write(15,'(a,a)') "   ",qmdff_en_line 
      write(15,'(a)') "# Speficy the TREQ method as PES description."
      write(15,'(a)') " pes treq "
      write(15,*)
      write(15,'(a)') " treq { "
      write(15,'(a)') "# file with xyz structures of the IRC:"
      write(15,'(a)') "    irc_struc irc.xyz"
      write(15,'(a)') "# file with reference energies of the IRC:"
      write(15,'(a)') "    irc_ens irc_ens.dat"
      write(15,'(a)') "# number of TREQ G+H reference points:"
      write(15,'(a,i5)') "    points ",n_rp_pts
      write(15,'(a)') "# RP-EVB damping coefficient for the coupling:"
      write(15,'(a,f14.7)') "    rp_exp_coeff ",pre_exp
      write(15,'(a)') "# borders of the RP direct interpolation region:"
      write(15,'(a,2f12.7)') "    rp_evb_mid ",rp_mid_tot,rp_mid_trans
      if (manual_ints) then
         write(15,'(a)') "# The set of internal coordinates will be read in from file!:"
         write(15,'(a)') "    read_coord"
      end if 
      write(15,'(a)') " } "
      write(15,'(a)') "##################################"
      write(15,'(a)') "# global dynamics settings"
      write(15,'(a)') "##################################" 
      write(15,'(a)') "# MD timestep (fs):"  
      write(15,'(a,f12.7)') "    deltat ",dt
      write(15,*)  
      write(15,'(a)') " nvt {"
      write(15,'(a)') "# temperature of the simulations (K):"
      write(15,'(a,i10)') "    temp ",temps(i)
      write(15,'(a)') "# the termostat to be used:"
      write(15,'(a,a)') "    thermostat ",thermo
      write(15,'(a)') "# the andersen frequency (if used):"
      write(15,'(a,i8)') "     andersen_step", andersen_step
      write(15,'(a)') "# the Nose-hoover damping factor (if used):"
      write(15,'(a,f20.8)') "     nose_damp",nose_q
      write(15,'(a)') " }"
      write(15,'(a)') "##################################"
      write(15,'(a)') "# reactive system settings"
      write(15,'(a)') "##################################" 
      write(15,'(a)') " mecha {"
      write(15,'(a)') "# species of the reaction mechanism:"  
      write(15,'(a,a)') "    type ", umbr_type  
      write(15,'(a)') "# list of all reactants with atom numbers:"
      do j=1,num_eds
         if (j .eq. 1) then
            write(15,'(a,i1,40i4)') "    reactant",j,frag1(:)
         else if (j .eq. 2) then
            write(15,'(a,i1,40i4)') "    reactant",j,frag2(:)
         else if (j .eq. 3) then
            write(15,'(a,i1,40i4)') "    reactant",j,frag3(:)
         else if (j .eq. 4) then
            write(15,'(a,i1,40i4)') "    reactant",j,frag4(:)
         end if   
         
      end do
      write(15,'(a)') "# The asymptotic distance in the reaction (Angstroms):"
      write(15,'(a,f15.7)') "    dist_inf ",r_inf
      write(15,'(a)') "# list of forming bonds in the reaction:"
      write(15,'(a)',advance="no") "    bond_form "
      do j=1,n_form
         write(15,'(a,a)',advance="no") " ",trim(print_form(j))
      end do
      write(15,'(a)') " "
      write(15,'(a)') "# list of breaking bonds in the reaction:"
      write(15,'(a)',advance="no") "    bond_break "
      do j=1,n_break
         write(15,'(a,a)',advance="no") " ",trim(print_break(j))
      end do
      write(15,'(a)') " "
      write(15,'(a)') "# number of equivalent reaction paths:"
      write(15,'(a,i3)') "    n_paths ",npaths
      if ((umbr_type .eq. "CYCLOREVER") .or. (umbr_type .eq. "REARRANGE") .or.&
            &  (umbr_type .eq. "DECOM_1BOND")) then
         write(15,'(a)') " "
         write(15,'(a)') "# file with reactants structure for unimolecular reactions"
         write(15,'(a)') "      reactants_struc reactants.xyz "
      end if
      write(15,'(a)') " }"
      write(15,'(a)') " "
      write(15,'(a)') "###################################"
      write(15,'(a)') "# umbrella sampling settings"
      write(15,'(a)') "###################################"
      write(15,'(a)') " umbrella {"
      write(15,'(a)') "# umbrella force constant (a.u.):"
      write(15,'(a,f14.7)') "    bias",k_force
      write(15,'(a)') "# borders of the umbrella samplings (xi):"
      write(15,'(a,2f13.7)') "    bonds",umbr_lo,umbr_hi
      write(15,'(a)') "# distance between two umbrella windows(xi):"
      write(15,'(a,f13.7)') "    dist",umbr_dist
      write(15,'(a)') "# number of MD timesteps for structure generation:"
      write(15,'(a,i9)') "    gen_steps",gen_step
      write(15,'(a)') "# number of MD timesteps for umbrella equilibration:"
      write(15,'(a,i9)') "    equi_steps",equi_step
      write(15,'(a)') "# number of MD timesteps per umbrella sampling:"
      write(15,'(a,i9)') "    sample_steps",umbr_step
      write(15,'(a)') "# number of umbrella trajectories per window:"
      write(15,'(a,i9)') "    sample_trajs",umbr_traj
      write(15,'(a)') " }"
      write(15,'(a)') "###################################"
      write(15,'(a)') "# potential of mean force calculation settings"
      write(15,'(a)') "###################################"
      write(15,'(a)') " pmf { "
      write(15,'(a)') "# borders of integragion for PMF (xi):"
      write(15,'(a,2f13.7)') "    xi_range ",xi_min,xi_max
      write(15,'(a)') "# number of integration gridpoints for PMF:"
      write(15,'(a,i9)') "    bins ",nbins
      write(15,'(a)') "# reaction mechanism category:"
      write(15,'(a,a)') "    method ",pmf_method
      if ((pmf_minloc .eq. 'ZERO') .or. (pmf_minloc .eq. 'PMF_MIN')) then
         write(15,'(a)') "# method for calculation of lower free energy:"
         write(15,'(a,a)') "    minloc ",pmf_minloc
      end if
      write(15,'(a)') " } "
      write(15,'(a)') "###################################"
      write(15,'(a)') "# recrossing calculation settings"
      write(15,'(a)') "###################################"
      write(15,'(a)') " recross { "
      write(15,'(a)') "# number of MD timesteps for parent sampling:"
      write(15,'(a,i9)') "    equi_steps",recr_equi
      write(15,'(a)') "# total number of child trajectories:"
      write(15,'(a,i10)') "    child_total",child_tot
      write(15,'(a)') "# number of MD timesteps between two child spawnings:"
      write(15,'(a,i9)') "    child_interval",child_interv
      write(15,'(a)') "# number of child trajectories per spawning point:"
      write(15,'(a,i9)') "    child_perpoint" ,child_point
      write(15,'(a)') "# number of evolution timesteps for each child:"
      write(15,'(a,i9)') "    child_steps",child_evol
      if (.not. recross_check) then
         write(15,'(a)') "# no error checking will be done for the recrossing part.."
         write(15,'(a)') "    no_check"
      end if
      if (recross_mpi) then
         write(15,'(a)') "# The recrossing calculation will be parallelized"
         write(15,'(a)') "    mpi"
      end if
      write(15,'(a)') " }"
      write(15,'(a)') "###################################"
      write(15,'(a)') "# miscellaneous"
      write(15,'(a)') "###################################"
      write(15,'(a)') "# maximum number of allowed errors per MD calculation:" 
      write(15,'(a,i9)') "    max_error",err_max
      write(15,'(a)') "# energy tolerance over which an error is stated:"
      write(15,'(a,f14.7)') "    rpmd_en_tol",energy_tol
      close(15)
      if (kt_avg .gt. 1) then
         call chdir("..")
      end if  
   end do
   call chdir("..")
end do
!
!     TASK 7 : START k(T) CALCULATIONS FOR ALL TEMPERATURES
!
!     Start all calculations in the folders that were generated 
!     above! Only one calculation will be started at a time, 
!     the program will wait until the actual calculation is 
!     finished
!  
!
write(*,*)
write(*,*) "---- TASK 7: start k(T) calculations for all temperatures -------"
write(*,*)
write(*,*) "calc_rate.x calculations are parallelized quite well, therefore only "
write(*,*) " one will be started simultenously with all cores available!"


!
!     Take the number of processors as character
!
write(adum_nprocs,'(i3)') nprocs_tot
write(*,'(a,i3)') " Number of needed reference calculations: ",n_rp_pts
write(*,'(a,i3)') " Total processors available: ",nprocs_tot
write(*,'(a,i2)') " Processors per calculation: ",nprocs_tot

do i=1,temp_num
!
!    Generate folder name
!
   if (temps(i) .lt. 10) then
      write(adum,'(i1)') temps(i)
   else if (temps(i) .lt. 100) then
      write(adum,'(i2)') temps(i)
   else if (temps(i) .lt. 1000) then
      write(adum,'(i3)') temps(i)
   else
      write(adum,'(i4)') temps(i)
   end if
   write(adum,'(a,a)') adjustl(trim(adum)),"K"
   call chdir(adjustl(trim(adum)))
!
!    Loop over all subfolders if existent 
!
   do k=1,kt_avg
      if (kt_avg .gt. 1) then
         if (k .gt. 9) then
            write(adum2,'(a,i2)') "run",k
         else
            write(adum2,'(a,i1)') "run",k
         end if
         call chdir(adjustl(trim(adum2)))
      end if

      inquire(file="done", exist=exist)
      if (.not. exist) then
         if (kt_avg .eq. 1) then
            write(*,'(a,i4,a)') " Black box k(T) calculation started for T= ",temps(i)," Kelvin!"
         else 
            write(*,'(a,i4,a,i2,a,i2)') " Black box k(T) calculation started for T= ",temps(i), &
                    &" Kelvin! (run ",k," of ",kt_avg,")"
         end if
         sys_stat=system(trim(link_mpi) // trim(adum_nprocs) // " " &
                         & //trim(link_rpmd)// " calc_rate.key > calc_rate.log")
         write(*,'(a)') " ---> calculation finished!"
         call system("touch done")
         if (sys_stat .ne. 0) then
            write(*,'(a,i4,a)') " ERROR! The Black box k(T) calculation at T= ",temps(i)," K failed!"
            call fatal
         end if
      else
         if (kt_avg .eq. 1) then
            write(*,'(a,i4,a)') " Black box k(T) calculation already done for T= ",temps(i), &
                    & " K ..."
         else 
            write(*,'(a,i4,a,a,a)') " Black box k(T) calculation already done for T= ",temps(i), &
                    & " K ,",trim(adum2)," ..."            
         end if
      end if
      if (kt_avg .gt. 1) then
         sys_stat=chdir("..")
      end if
   end do
   sys_stat=chdir("..")
end do
!
!     TASK 8 : EVALUATE THE K(T) CALCULATIONS 
!
!     Evaluate the RPMD k(T) calculations that were done above!
!     Read in the rate constants in both units.
!     if more than one calculation was done, calculate also
!     Arrhenius parameters from a linear regression
!
kt_error=0
write(*,*)
write(*,*) "---- TASK 8: evaluate k(T) calculation results -------"
write(*,*)
write(*,*) "Read in calculated rate constants for all temperatures."
write(*,*) "If two or more different temperatures were calculated, perform"
write(*,*) "also a linear regression for calculation of Arrhenius"
write(*,*) "parameters!"

allocate(rates_mol(temp_num,kt_avg))
allocate(rates_molec(temp_num,kt_avg))
allocate(log_kt_mol(temp_num,kt_avg))
allocate(log_kt_molec(temp_num,kt_avg))
allocate(avg_kt_mol(temp_num))
allocate(avg_kt_molec(temp_num))
allocate(var_kt_mol(temp_num))
allocate(var_kt_molec(temp_num))
allocate(ln_kt(temp_num))
allocate(onet(temp_num))
write(*,*)
write(*,*) "The calculated reaction rate coefficients are:"

avg_kt_mol=0.d0
avg_kt_molec=0.d0
do i=1,temp_num
!
!    Generate folder name
!
   if (temps(i) .lt. 10) then
      write(adum,'(i1)') temps(i)
   else if (temps(i) .lt. 100) then
      write(adum,'(i2)') temps(i)
   else if (temps(i) .lt. 1000) then
      write(adum,'(i3)') temps(i)
   else
      write(adum,'(i4)') temps(i)
   end if
   write(adum,'(a,a)') adjustl(trim(adum)),"K"
   call chdir(adjustl(trim(adum)))
   rates_mol(i,:)=0.d0
   rates_molec(i,:)=0.d0
!
!     If more than one rate was calculated per temperature, enter the inner 
!     loop to go into all these folders 
!
   do k=1,kt_avg
      if (kt_avg .gt. 1) then
         if (k .gt. 9) then
            write(adum2,'(a,i2)') "run",k
         else
            write(adum2,'(a,i1)') "run",k
         end if
         call chdir(adjustl(trim(adum2)))
      end if
!
!     For the first calculation, copy the written file qmdff_ref.dat
!     with the QMDFF energies for the reactionpath to the main folder 
!     mep_irc in order to make a plot possible there 
!
      if ((i .eq. 1) .and. (k .eq. 1)) then
         write(*,*) "INFO: file qmdff_ref.dat was copied to mep_irc folder."
         write(*,*) " this file contains plots of QMDFF single energies!"
         call system("cp qmdff_ref.dat ../../../mep_irc")
      end if
!
!     Open the qmdff.log output and search for the obtained 
!     rate constants       
!  

      open(unit=16,file="calc_rate.log",status="old")
      do 
         read(16,'(a)',iostat=lastline) a80
         if (lastline .ne. 0) exit
         if (index(a80,'The value of the reaction rate constant is:') .ne. 0) then
            read(16,*) rates_mol(i,k)
            read(16,*) rates_molec(i,k)    
            exit
         end if
      end do
      if (rates_mol(i,k) .eq. 0d0) then 
         if (kt_avg .eq. 1) then
            write(*,*) "ERROR! No useful rate constant could be obtained for ",adum,"!"
         else 
            write(*,*) "ERROR! No useful rate constant could be obtained for ",adum," ",adum2,"!"
         end if
         write(*,*) "Check your input!"
         kt_error=1
      else 
         if (kt_avg .eq. 1) then
            write(*,'(a,a,a,es16.8,a,es16.8,a)') " ",trim(adum),":  ",rates_mol(i,k),& 
                &  " cm^3/(mol*s) , ", rates_molec(i,k)," cm^3/(molec*s)"
         else 
            write(*,'(a,a,a,a,a,es16.8,a,es16.8,a)') " ",trim(adum),",",trim(adum2),":  ",rates_mol(i,k),&
                &  " cm^3/(mol*s) , ", rates_molec(i,k)," cm^3/(molec*s)"
         end if
      end if
      close(16)
      if (kt_avg .gt. 1) then
         call chdir("..")
      end if
   end do
   call chdir("..")
end do
if (kt_avg .gt. 1) then
!
!     Determine averages in usual units for printout if two or more calculations 
!     are done for each temperature
!
   write(*,*)
   write(*,*) "Averaged k(T) values for all temperatures:"
   allocate(avg_kt_print(temp_num,2))
   avg_kt_print=0d0
   do i=1,temp_num
      if (temps(i) .lt. 10) then
         write(adum,'(i1)') temps(i)
      else if (temps(i) .lt. 100) then
         write(adum,'(i2)') temps(i)
      else if (temps(i) .lt. 1000) then
         write(adum,'(i3)') temps(i)
      else
         write(adum,'(i4)') temps(i)
      end if
      write(adum,'(a,a)') adjustl(trim(adum)),"K"
      do j=1,kt_avg
         avg_kt_print(i,1)=avg_kt_print(i,1)+rates_mol(i,j)
         avg_kt_print(i,2)=avg_kt_print(i,2)+rates_molec(i,j)      
      end do
      avg_kt_print(i,1)=avg_kt_print(i,1)/kt_avg
      avg_kt_print(i,2)=avg_kt_print(i,2)/kt_avg
      write(*,'(a,a,a,es16.8,a,es16.8,a)') " ",trim(adum),":  ",avg_kt_print(i,1),&
                &  " cm^3/(mol*s) , ", avg_kt_print(i,2)," cm^3/(molec*s)"
   end do
end if
!
!     Determine logarithms of all rate constants 
!
do i=1,temp_num
   do j=1,kt_avg
      log_kt_mol(i,j)=log(rates_mol(i,j))
      log_kt_molec(i,j)=log(rates_molec(i,j))
   end do
end do
!
!     calculate total loratihmaverage by division with number of runs 
!
do i=1,temp_num
   do j=1,kt_avg
      avg_kt_mol(i)=avg_kt_mol(i)+log_kt_mol(i,j)
      avg_kt_molec(i)=avg_kt_molec(i)+log_kt_molec(i,j)
   end do
end do
avg_kt_mol=avg_kt_mol/kt_avg
avg_kt_molec=avg_kt_molec/kt_avg
!     
!     calculate variances for all temperaturs by applying the usual formula
!
var_kt_mol=0.d0
var_kt_molec=0.d0
do i=1,temp_num
   do j=1,kt_avg
      var_kt_mol(i)=var_kt_mol(i)+(log_kt_mol(i,j)-avg_kt_mol(i))*&
                & (log_kt_mol(i,j)-avg_kt_mol(i))
      var_kt_molec(i)=var_kt_molec(i)+(log_kt_molec(i,j)-avg_kt_molec(i))*&
                & (log_kt_molec(i,j)-avg_kt_molec(i)) 
   end do
end do
var_kt_mol=sqrt(var_kt_mol)/kt_avg
var_kt_molec=sqrt(var_kt_molec)/kt_avg
!
!     determine ln(k(T)) and 1/T values for Arrhenius plot
!        
do i=1,temp_num 
   ln_kt(i)=avg_kt_mol(i)
   onet(i)=real(1.d0/real(temps(i)))
end do
!
!     Now do the linear regression for Arrhenius parameters if 
!     enough rates were calculated!
!     Write also files for visualization of plots 
!
if (temp_num .ge. 2) then
   if (kt_error .eq. 0) then
      open(unit=15,file="arrhenius.dat",status="unknown")
      do i=1,temp_num
         write(15,*) onet(i),ln_kt(i),var_kt_mol(i)
      end do
      close(15)

      open(unit=15,file="all_kt.dat",status="unknown")
      do i=1,temp_num
         do j=1,kt_avg
            write(15,*) onet(i),log_kt_mol(i,j)
         end do
      end do
      close(15)      

      write(*,*)
      write(*,*) "Do Arrhenius fit for calculated rates:"

      call lin_reg(onet,ln_kt,temp_num,arrh_slope,arrh_y,corre) 
 
      open(unit=15,file="arrhen_plot.gnu",status="unknown")
      write(15,*) "set title 'Arrhenius-Plot with calculated data:'"
      write(15,*) "set xlabel '1/T (K^-1)'"
      write(15,*) "set ylabel 'ln(k(T))'"
      write(15,*) "set xrange [",onet(temp_num)*(0.9d0),":",onet(1)*(1.1d0),"]"
      write(15,*) "plot 'all_kt.dat' u 1:2 with points pointtype 2 pointsize 1.3 title &
             & 'calculated k(T) data','arrhenius.dat' u 1:2:3 with yerrorbars pointtype 5 &
             & pointsize 1.2 title 'averaged k(T) data', ",arrh_y,"+",&
             & arrh_slope,"*x title 'Arrhenius fit'"
      write(15,*) "pause -1"
      close(15)


      arrh_y=exp(arrh_y)/6.023E23
      arrh_slope=arrh_slope!*8.31447/1000d0
 
      write(*,*) " --> Done! The results are:"
      write(*,*)
      write(*,'(a,f16.7,a,es16.7,a)') " slope: ",-arrh_slope," K, y_intercept: ",&
               & arrh_y," cm^3/(mol*s)"
      arrh_y=arrh_y!exp(arrh_y)/6.023E23
      arrh_slope=arrh_slope*8.31447/1000d0
      write(*,*)
      write(*,*) "or in other units (NIST):"
      write(*,*) 
      write(*,'(a,f16.7,a,es16.7,a)') " Ea (slope): ",-arrh_slope," kJ/mol, A(y_intercept): ",&
             & arrh_y," cm^3/(mol*s)"
      write(*,*) 
      write(*,'(a,f14.8)') " The correlation coefficient:",corre
      write(*,*) 
      write(*,*) "In order to visualize the Arrhenius plot, open: calc_rate/arrhen_plot.gnu"
   else 
      write(*,*) "No Arrhenius fit could be done because an error occured during "
      write(*,*) " a k(T) calculation."
   end if
end if


call system_clock(time_int(7))
!
!    Print out the calculation time partitioning and the final messages
!    Print out time measuring for better informations 
!

do i=1,10
   time(i)=real(time_int(i))/1000
end do

tot_time=int(time(7)-time(1))
ndays=tot_time/86400
nhours=(tot_time-86400*ndays)/3600
nminutes=(tot_time-86400*ndays-3600*nhours)/60
nseconds=tot_time-86400*ndays-3600*nhours-60*nminutes

write(*,*)
write(*,*) " Timings: "
write(*,*) " ----------"
write(*,'(A, F12.3, A)') " Read in of settings and initialization:  ",time(2)-time(1)," s."
write(*,'(A, F12.3, A)') " QMDFF reference calculation:             ",time(3)-time(2)," s."
write(*,'(A, F12.3, A)') " QMDFF generation:                        ",time(4)-time(3)," s."
write(*,'(A, F12.3, A)') " TREQ reference data calculation:         ",time(5)-time(4)," s."
write(*,'(A, F12.3, A)') " Additional energy calculations:          ",time(6)-time(5)," s."
write(*,'(A, F12.3, A)') " Rate constant calculations + Arrhenius:  ",time(7)-time(6)," s."
write(*,*)
write(*,'(A, F12.3, A)') " The calculation needed a total time of   ",time(7)-time(1)," seconds."
write(*,'(A, I6,A,I2,A,I2,A,I2,A)') " These are: ",ndays," days, ",nhours," hours, ",nminutes,&
          & " minutes and ",nseconds," seconds."
write(*,*)



write(*,*) 
write(*,*) " -.-. .- .-. .- -.-. .- .-.."
write(*,*) "CARACAL has successfully finished all tasks!"
write(*,*) "Look into the results folders for more details and plots."
write(*,*) 

sys_stat=chdir("..")
end program black_box

    
