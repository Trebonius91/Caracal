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
!     ##  program evbopt  --  optimizing the evb-coupling-element    ##
!     ##                                                             ##
!     #################################################################
!
!
!     "evbopt" performs the automatical optimization of evb
!     off-diagonal coupling elements. 
!     For this purpose Levenberg-Marquardt optimizations of different
!     basisfunctions and the DG-EVB
!     couplingterms are availiable
!     2x2 and 3x3 EVBs can be handled (for dE and dQ)
!     Now, this subroutine can also executed in parallel using MPI
!     

program evbopt
use general
use evb_mod
use pbc_mod

implicit none
!
!     include MPI library
!
include 'mpif.h'
integer::num_arg,input_unit,i,qmdff_energies_unit,asciinum
integer::qmdffnumber
integer::dg_evb_mode
integer::reference_counter,parameters  !zahl der Strukturen im xyz-file
real(kind=8)::energy,e_qmdff1,e_qmdff2,e_qmdff3
real(kind=8),dimension(:,:),allocatable::coord,energies_qmdff !array für energien der einzelnen FFs
real(kind=8),dimension(:,:),allocatable::energies_tmp  ! für temporäres zwischenspeichern
!    for correction of QMDFF energies in order to reproduce minimum energies exactly
real(kind=8),dimension(:),allocatable::refens
real(kind=8),dimension(:,:,:),allocatable::struc_tmp 
real(kind=8)::diff1,diff2,e_one,e_two
real(kind=8),dimension(:),allocatable::energies_result
character(len=70)::fffile1,fffile2,fffile3,fileinfo,filegeo,fileenergy,filets,filets2
character(len=60)::test,names
character(len=50)::method  ! the EVB coupling method
character(len=100)::syscall
real(kind=8),dimension(:,:),allocatable::ts_coordinates_a,ts_coordinates2_a
character(len=1)::qmdffnum
character(len=20) keyword
character(len=120) record
character(len=120) string
character(len=:),allocatable::char_num ! for the number of reference points in dg_evb
character(len=120) coupling
character(len=30) write_format
character(len=40)::commarg ! string for command line argument
integer mode,next,readstatus,j,k,nat
logical::exist,exists,has_next
logical::evb2,evb3,ffname2,ffname3,defqmdff
logical::path_struc,path_energy,coupl,ts_xyz
!     for MPI parallelization
integer::ierr,ID,psize
integer::rank
integer::readstat
!
!     Start MPI parallel computation:
!     Since the whole program is executed by all processes, all "serial"
!     parts will be executed only by processor #0
!
call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,psize,ierr) ! total number of procs
call mpi_comm_rank(mpi_comm_world,rank,ierr) ! index of current proc

!
!     The explore program is not used
!
use_explore = .false.
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
      call help("evbopt") 
      stop
   else 
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'evbopt.x -help' or 'evbopt.x -h'."
      write(*,*)
   end if
end if
!
!     Read keywords from generic key-file "tinker.key"
!
call getkey(rank)
evb2=.false.
evb3=.false.
ffname2=.false.
ffname3=.false.
path_struc=.false.
path_energy=.false.
coupl=.false.
defqmdff=.false.
evb_dq=.false.
!
!     Open the global output-file for EVB-term optimizations
!
if (rank .eq. 0) then
   open(unit=15,file="evbopt.log",status="unknown")
   call promo_log(15)
end if

!
!     Try to read the needed parameters for the Optimization
!
!     The general Method keyword
!
method=""
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'PES ') then
      read(record,*) names,method
      call upcase(method)
      exit
   end if
end do

!
!      Choose the EVB coupling term to be optimized
!
evb_dq=.false.
evb_de=.false.
dg_evb=.false.
treq=.false.
if (method .eq. "DE_EVB")  then
   write(*,*) "The dE-EVB (energy-gap) coupling term will be used!"
   write(*,*) 
   evb_de=.true.
else if (method .eq. "DQ_EVB")  then
   evb_dq=.true.
   write(*,*) "The dQ-EVB (coordinate) coupling term will be used!"
   write(*,*)
else if (method .eq. "DG_EVB") then
   dg_evb=.true.
   write(*,*) "The DG-EVB (distributed Gaussian) coupling term will be used!"
   write(*,*)
else if (method .eq. "TREQ") then
   treq=.true.
   write(*,*) "The TREQ (transition region corrected reaction path) PES" 
   write(*,*) " description will be used!"
   write(*,*)
else 
   write(*,*) "No valid EVB coupling term was given!"
   write(*,*) "Add the PES keyword with either 'de_evb', 'dq_evb', 'dg_evb'"
   write(*,*) "  or 'treq' as parameter!"
   call fatal
end if
!
!    Read in parameters for the diabatic QMDFF surfaces
!
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)

   if (trim(adjustl(record(1:11))) .eq. 'QMDFF {' .or. trim(adjustl(record(1:11))) &
        &  .eq. 'QMDFF{') then

      do j=1,nkey_lines-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
!
!     Read in the names of the QMDFF files
!
         if (keyword(1:11) .eq. 'FFNAMES ') then
            if (qmdffnumber .eq. 1) then
               read(record,*,iostat=readstat) names,fffile1
               if (readstat .ne. 0) then
                  write(*,*) "Please check the QMDFFNAMES keyword!"
                  call fatal
               end if
            end if
            read(record,*,iostat=readstat) names,fffile1,fffile2,fffile3
            ffname3=.true.
            defqmdff=.true.
            if (readstat .ne. 0) then
               read(record,*,iostat=readstat) names,fffile1,fffile2
               ffname2=.true.
               defqmdff=.true.
               if (readstat .ne. 0) then
                  qmdffnumber=1
               else if (readstat .eq. 0) then
                  evb2=.true.
                  qmdffnumber=2
               end if
            else if (readstat .eq. 0) then
               evb3=.true.
               qmdffnumber=3
            end if
!
!     Read in the relative energy shift of the QMDFFs 
!
         else if (keyword(1:11) .eq. 'ESHIFT ') then
            if (qmdffnumber .eq. 1) then
               read(record,*,iostat=readstat) names,E_zero1
               exist=.true.
               if (readstat .ne. 0) then
                  write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                  call fatal
               end if
            end if
            if (evb2) then
               read(record,*,iostat=readstat) names,E_zero1,E_zero2
               exist=.true.
               if (readstat .ne. 0) then
                  write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                  call fatal
               end if
            end if
            if (evb3) then
               read(record,*,iostat=readstat) names,E_zero1,E_zero2,E_zero3
               exist=.true.
               if (readstat .ne. 0) then
                  write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                  call fatal
               end if
            end if
!           
!     Determine if the QMDFF energy shifts shall be corrected automatically
!     to exactly reproduce the first/last energy of the path or not
!    
         else if (keyword(1:16) .eq. 'SHIFT_MANUAL') then
            shift_man=.true.
         end if
         if (keyword(1:11) .eq. '}') exit
         if (j .eq. nkey_lines-i) then
            write(*,*) "The QMDFF section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do

! 
!     In the case that no QMDFF´s are defined in key file
!     read in dem manually 
!

if (.not. ffname2 .and. .not.  ffname3) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if
   exist=.false.
   do while (.not. exist)
      write(iout,24)
  24  format (/,' Number of used QMDFF´s:  ',$)
      read (*,34)  qmdffnum
  34  format (a1)
      asciinum = ICHAR(qmdffnum)
!
!     Convert vom ASCII to the real number (1=49,2=50,3=51)
!
      select case (asciinum)
      case (50)
         qmdffnumber=2
      case (51)
         qmdffnumber=3
      end select
      if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
         exist=.true.
      end if
   end do
end if

if (.not.defqmdff) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if
   if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,21)
 21      format (/,' Filename of the first QMDFF:  ',$)
         read (*,31)  fffile1
 31      format (a120)
         inquire(file=fffile1,exist=exist)
      end do
 
      exist=.false.
      do while (.not. exist)
         write(iout,22)
 22      format (/,' Filename of the second QMDFF:  ',$)
         read (*,32)  fffile2
 32      format (a120)
         inquire(file=fffile2,exist=exist)
      end do
   end if
   if (qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,23)
 23      format (/,' Filename of the third QMDFF:  ',$)
         read (*,33)  fffile1
 33      format (a120)
         inquire(file=fffile1,exist=exist)
      end do
   end if
end if
!
!     Read in the QMDFF-energys if not defined via keyfile
!
if ((.not. evb2) .and. (.not. evb3)) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if

   if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
      write(iout,27)
 27   format (/,' Shift-Energy of the first QMDFF :  ',$)
      read (*,*) E_zero1
      write(iout,28)
 28   format (/,' Shift-Energy of the second QMDFF :  ',$)
      read (*,*) E_zero2
   end if
   if (qmdffnumber.eq.3) then
      write(iout,29)
 29   format (/,' Sjift-Energy of the third QMDFF :  ',$)
      read (*,*) E_zero3
   end if
end if
!
!    Read in the structures of a reactionpath
!
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'COORDS_REF ') then
      read(record,*) names,filegeo
      inquire(file=filegeo,exist=path_struc)
   end if
end do
if (.not. path_struc) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if

   exist=.false.
   do while (.not. exist)
      write(iout,25)
 25   format (/,' Name of the xyz-File with Structures:  ',$)
      read (*,35)  filegeo
 35   format (a120)
      inquire(file=filegeo,exist=exist)
   end do
end if
!
!    Read in the energies of a reactionpath
!

do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'ENERGIES_REF ') then
      read(record,*) names,fileenergy
      inquire(file=fileenergy,exist=path_energy)
   end if
end do
if (.not. path_energy) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if

   exist=.false.
   do while (.not. exist)
      write(iout,26)
 26   format (/,' Name of the File with energies of reactionpath:  ',$)
      read (*,36)  fileenergy
 36   format (a120)
      inquire(file=fileenergy,exist=exist)
   end do
end if
!
!    For the levenberg-marquardt-algorithms: read in several parameters (see evb_mod.f90)
!
optsteps=100
maxstart=100
lower_bond=1
upper_bond=30
lm_par=0.01d0
lm_threshold=1E-8
lm_par_change=2d0
diff_step=0.0001d0
num_coord=5 ! number of internal coordinates
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'MAXSTEP ') then
      read(record,*) names,optsteps
   else if (keyword(1:16) .eq. 'LM_THRESHOLD ') then
      read(record,*) names,lm_threshold
   else if (keyword(1:16) .eq. 'RANDOM_BONDS ') then
      read(record,*) names,lower_bond,upper_bond
   else if (keyword(1:16) .eq. 'DIFF_STEP ') then
      read(record,*) names,diff_step
   else if (keyword(1:16) .eq. 'START_POINTS ') then
      read(record,*) names,maxstart
   else if (keyword(1:16) .eq. 'READ_COORD ') then
      read_coord=.true.
   else if (keyword(1:16) .eq. 'NUM_COORD ') then
      read(record,*) names,num_coord
   end if
end do

!
!      Measure the needed time for optimization
!
call system_clock(time_int(1))
!
!      Read in the detailed settings for the dE-EVB coupling term
!      Loop over structure with brackets
!
if (evb_de) then
!
!      Set default values 
!
   off_basis="1g"
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:15))) .eq. 'DE_EVB { ' .or. &
                  & trim(adjustl(record(1:15))) .eq. 'DE_EVB{ ') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11) .eq. 'COUPLING ') then
               read(record,*) names,off_basis
               if (off_basis .eq. "1g" .or. off_basis.eq."2g" & 
                   & .or. off_basis.eq."3g" &
                   & .or.off_basis.eq."sp" .or. off_basis.eq."sd" &
                   & .or. off_basis.eq."sd2" .or. off_basis.eq."sp2d3" &
                   & .or. off_basis.eq."sp2d"  ) then
               else
                  if (rank .eq. 0) then
                     write(*,*) "No valid coupling function was defined."
                     write(*,*) "We will use the simple 1g-dE coupling instead."
                  end if
                  off_basis="1g"
               end if
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The DE-EVB section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do
end if
!
!      Read in the detailed settings for the dE-EVB coupling term
!      Loop over structure with brackets
!
if (evb_dq) then
!
!      Set default values 
!
   off_basis="1g"
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:15))) .eq. 'DQ_EVB { ' .or. &
                  & trim(adjustl(record(1:15))) .eq. 'DQ_EVB{ ') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:11) .eq. 'COUPLING ') then
               read(record,*) names,off_basis
               if (off_basis .eq. "1g" .or. off_basis.eq."2g" &
                   & .or. off_basis.eq."sd2") then
               else
                  if (rank .eq. 0) then
                     write(*,*) "No valid coupling function was defined."
                     write(*,*) "We will use the simple 1g-dQ coupling instead."
                  end if
                  off_basis="1g"
               end if
            else if (keyword(1:11) .eq. 'TS_FILE ') then
               read(record,*) names,filets
               inquire(file=filets,exist=ts_xyz)
            end if


            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The DE-EVB section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do
end if

!
!      Read in the detailed settings for the DG-EVB coupling term
!      Loop over structure with brackets
!

if (dg_evb) then
!
!      Set default values 
!
   dg_evb_mode=0
   dg_evb_points=0
   dg_ref_file="grad_hess.dat"
   g_thres=1E-10
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:15))) .eq. 'DG_EVB { ' .or. &
                  &  trim(adjustl(record(1:15))) .eq. 'DG_EVB{ ') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            if (keyword(1:16) .eq. 'MODE ') then
               read(record,*) names,dg_evb_mode
            else if (keyword(1:16) .eq. 'POINTS ') then
               read(record,*) names,dg_evb_points
            else if (keyword(1:16) .eq. 'DOUBLE_ALPHA ') then
               double_alpha=.true.
            else if (keyword(1:16) .eq. 'READ_COORD ') then
               read_coord=.true.
            else if (keyword(1:16) .eq. 'POINTS_REF ') then
               read(record,*) names,dg_ref_file
            else if (keyword(1:16) .eq. 'DG_MAT_NUM ') then
               dg_mat_num=.true.
            else if (keyword(1:18) .eq. 'GAUSS_THRESHOLD ') then
               read(record,*) names,g_thres 
            end if
            if (record .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The DG-EVB section has no second delimiter! (})"
               call fatal
            end if   
         end do
      end if
   end do
   if (double_alpha) then
      add_alph=dg_evb_points
   end if
!
!     Catch missing information
!
    inquire(file=dg_ref_file,exist=exist)
    if (.not. exist) then
       if (rank .eq. 0) then
          dg_evb=.false.
          write(*,*) "The file ", trim(dg_ref_file)," containing the DG-ref. informations"
          write(*,*) "is not availiable!"
          call fatal
       end if
    end if

   if (dg_evb_mode .eq. 0) then
      write(*,*) "Please add the keyword 'MODE' to the DG_EVB section!"
      call fatal
   end if
   if (dg_evb_points .eq. 0) then
      write(*,*) "Please add the keywod 'POINTS' to the DG-EVB section!"
      call fatal
   end if
end if
!
!     In case of a reaction path EVB (RP-EVB), check if the method shall be 
!     activated
!
!do i = 1, nkey_lines
!   next = 1
!   record = keyline(i)
!   call gettext (record,keyword,next)
!   call upcase (keyword)
!   string = record(next:120)
!   if (keyword(1:11) .eq. 'TREQ ') then
!      read(record,*) names,rp_evb_points
!      treq=.true.
!      names="grad_hess.dat"
!      inquire(file=names,exist=exist)
!      if (rank .eq. 0) then
!         if (.not. exist) then
!            treq=.false.
!            write(*,*) "You have ordered a TREQ calculation,"
!            write(*,*) "but the file ", trim(names)," containing the reference informations"
!            write(*,*) "is not availiable!"
!            write(*,*) "Look into the manual for further informations."
!            call fatal
!         end if
!      end if
!   end if
!end do
if (evb_dq .and. .not. ts_xyz) then
   if (psize .gt. 1) then
      if (rank .eq. 0) then
         write(*,*) "Please read in the reference data from file if you want"
         write(*,*) "to do a parallel calculation! Look into the manual for details."
         call fatal
      end if
   end if

   exist=.false.
   do while (.not. exist)
      write(iout,56)
 56   format (/,' No File with TS-structure for evb_dq-calculation! &
              Please enter the name::  ',$)
      read (*,66)  filets
 66   format (a120)
      inquire(file=filets,exist=exist)
   end do
end if
!
!     initialize the evb-qmdff-force-field
!
call prepare (fffile1,fffile2,fffile3,qmdffnumber)
allocate(coord(3,natoms))
parameters=20
nat=natoms

!
!     first produce the QMDFF-energies for all points
!
if (qmdffnumber.eq.2) then
   input_unit=7
   reference_counter=-1
   allocate(energies_tmp(10000,2))
   allocate(struc_tmp(10000,3,natoms))
   open(unit=input_unit,file=filegeo,status='old')
   do
      reference_counter=reference_counter+1
      call next_geo(coord,natoms,input_unit,has_next)
      struc_tmp(reference_counter+1,:,:)=coord
      if (.not.has_next) exit
      call eqmdff(coord,e_qmdff1,e_qmdff2)
      energies_tmp(reference_counter+1,1)=e_qmdff1
      energies_tmp(reference_counter+1,2)=e_qmdff2
   end do
   allocate(energies_qmdff(reference_counter,2))
!
!     New feature (10.07.2018): shift QMDFF energies in order to reproduce exactly 
!     the first and and last energy of the path with the respective QMDFF
!
   if (.not. shift_man) then
      allocate(refens(reference_counter))
      open(unit=56,file=fileenergy,status="unknown")
      do i=1,reference_counter
         read(56,*) refens(i)
      end do
      diff1=abs(refens(1)- energies_tmp(1,1))
      diff2=abs(refens(reference_counter)-energies_tmp(reference_counter,1))
!
!     if the first QMDFF describes the first minimum
!
      if (diff1 .lt. diff2) then
         e_one=energies_tmp(1,1)-E_zero1
         E_zero1=refens(1)-e_one

         e_two=energies_tmp(reference_counter,2)-E_zero2
         E_zero2=refens(reference_counter)-e_two
!
!    if the first QMDFF describes the second minimum
!   
      else
         e_two=energies_tmp(1,2)-E_zero2
         E_zero2=refens(1)-e_two

         e_one=energies_tmp(reference_counter,1)-E_zero1
         E_zero1=refens(reference_counter)-e_one
      end if
      close(56)
      if (rank .eq. 0) then
         write(*,*) "New QMDFF shifts were calculated to exacly reproduce asymptotics of path."
         write(*,*) "Command 2EVB in keyfile was updated.."
      end if
      write(syscall,'(a,f24.14,f24.14,a,a)') 'sed -i "s/ 2EVB  .*/ 2EVB  ',E_zero1,&
                  & E_zero2,'/g"  ',trim(keyfile) 
      call system(syscall)
      deallocate(refens)
   else 
      if (rank .eq. 0) then
         write(*,*) "No automatic corrections on the QMDFF shifts will be applied."
         write(*,*) "To activate these corrections, remove the keyword SHIFT_MANUAL."
      end if
   end if
   

   if (rank .eq. 0) then 
      open (unit=99,file="single_qmdff.dat",status='unknown')
   end if
!
!  write the single-QMDFF-energies
!
   do i=1,reference_counter 
      call eqmdff(struc_tmp(i,:,:),e_qmdff1,e_qmdff2)
      energies_qmdff(i,1) = e_qmdff1
      energies_qmdff(i,2) = e_qmdff2
      if (rank .eq. 0) then
         write(99,*) energies_qmdff(i,1), energies_qmdff(i,2)
      end if
   end do   
   if (rank .eq. 0) then
      close(99)
   end if
   close(7)
   deallocate(energies_tmp)
   deallocate(struc_tmp)
else if (qmdffnumber.eq.3) then
   input_unit=7
   reference_counter=-1
   allocate(energies_tmp(400,3))
   open(unit=input_unit,file=filegeo,status='old')
   do
      reference_counter=reference_counter+1
      call next_geo(coord,natoms,input_unit,has_next)
      if (.not.has_next) exit
      call eqmdff3(coord,e_qmdff1,e_qmdff2,e_qmdff3)
      energies_tmp(reference_counter+1,1)=e_qmdff1
      energies_tmp(reference_counter+1,2)=e_qmdff2
      energies_tmp(reference_counter+1,3)=e_qmdff3
   end do
   allocate(energies_qmdff(reference_counter,3))
   if (rank .eq. 0) then
      open (unit=99,file="single_qmdff.dat",status='unknown')
   end if
!
!  write the single-QMDFF-energies
!
   do i=1,reference_counter
      energies_qmdff(i,1) = energies_tmp(i,1)
      energies_qmdff(i,2) = energies_tmp(i,2)
      energies_qmdff(i,3) = energies_tmp(i,3)
      if (rank .eq. 0) then
         write(99,*) energies_qmdff(i,1), energies_qmdff(i,2), energies_qmdff(i,3)
      end if
   end do
   if (rank .eq. 0) then
      close(99)
   end if
   close(7)
   deallocate(energies_tmp)
end if
!
!     If you want to give start values for the
!     coupling-term parameters manually
!
coupl=.false.
if (qmdffnumber.eq.2) then
!
!     Optimize the dE-off-diagonal term with reference-data
!     by multi start local search with Levenberg-Marquardt
!     --> only serial calculation
!
   if (rank .eq. 0) then
      if ((.not.evb_dq) .and. (.not.dg_evb) .and. (.not.treq)) then
         allocate(energies_result(reference_counter))
         call optimize_dE(energies_qmdff,fileenergy,reference_counter,energies_result)

!   write the EVB-QMDFF-energies for the optimized parameters
!
         open(unit=17,file="energies_final.dat",status="unknown")
         do i=1,reference_counter 
            write(17,*) energies_result(i)
         end do      
         close(17)
      end if
   end if

!
!     This part for the dq-couplingterm
!
!TEST manual ref_counter
   if (rank .eq. 0) then
      if (evb_dq) then
         allocate(energies_result(reference_counter))
         call optimize_dQ(energies_qmdff,natoms,fileenergy,filegeo,filets,&
              reference_counter,energies_result)

         open(unit=17,file="energies_final.dat",status="unknown")
         do i=1,reference_counter
            write(17,*) energies_result(i)
         end do
         close(17)
      end if
   end if
 
!
!    Distributed gaussian EVB (DG-EVB)
!
   if (dg_evb) then 
!
!   default values
!
     maxstart=100
     diff_ana=.false.
     lower_bond=1
     upper_bond=30
     read_coord=.false.
     num_coord=5
     double_alpha=.false.
     more_info=.false.
     do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:16) .eq. 'START_POINTS ') then
            read(record,*) names,maxstart
         end if
         if (keyword(1:16) .eq. 'DIFF_ANA ') then
            diff_ana=.true.
         end if
         if (keyword(1:16) .eq. 'RANDOM_BONDS ') then
            read(record,*) names,lower_bond,upper_bond
         end if
         if (keyword(1:16) .eq. 'READ_COORD ') then
            read_coord=.true.
         end if
         if (keyword(1:16) .eq. 'NUM_COORD ') then
            read(record,*) names,num_coord
         end if
         if (keyword(1:20) .eq. 'MORE_INFO') then
             more_info=.true.
         end if


      end do
      if ((dg_evb_mode .gt. 3) .or. (dg_evb_mode .lt. 1)) then
         if (rank .eq. 0) then
            write(*,*) "You have chosen the DG-EVB-Mode",dg_evb_mode
            write(*,*) "Only the modes 1 to 3 are availiable!"
            call fatal
         end if
      end if 
!
!     Call the DG-EVB optimization routine, parallel or serial version
!
      call mpi_barrier(mpi_comm_world,ierr)
      if (psize .eq. 1) then
         call dg_evb_init_ser(dg_evb_mode,filegeo,fileenergy,reference_counter)
      else 
         call dg_evb_init_par(dg_evb_mode,filegeo,fileenergy,reference_counter,psize,rank)
      end if
   end if
!
!    Transition corrected reaction path EVB (TREQ)  
!
   if (treq) then
       if (rank .eq. 0) then
          stop "No optimization of parameters needed for TREQ!"
       end if
   end if
   
else if (qmdffnumber.eq.3) then
   if (off_basis .ne. "const" .and. off_basis .ne. "1g" .and. off_basis .ne. "sd2") then
      write(*,*) "The coupling-basis ",off_basis,"isn´t implemented yet for 3x3-EVB!"
      call fatal
   end if
!
!   default values
!
   maxstart=100
   lower_bond=1
   upper_bond=30

   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:16) .eq. 'START_POINTS ') then
         read(record,*) names,maxstart
      end if
      if (keyword(1:16) .eq. 'DIFF_ANA ') then
         diff_ana=.true.
      end if
      if (keyword(1:16) .eq. 'RANDOM_BONDS ') then
         read(record,*) names,lower_bond,upper_bond
      end if
   end do

!
!    Decide if the 1-3 EVB-coupling-parameters should be neglect
!
   full=.false.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'FULL_MATRIX ') then
         full=.true.
      end if
   end do
!
!    If the DQ-coupling is desired, read in the coordinates of the both
!    transition-states and the reactionpath

   do i = 1, nkey_lines
       next = 1
       record = keyline(i)
       call gettext (record,keyword,next)
       call upcase (keyword)
       string = record(next:120)
       if (keyword(1:11) .eq. 'EVB_DQ ') then
          read(record,*) names,filets,filets2
          evb_dq = .true.
       end if
   end do
   if (evb_dq .eqv. .true.) then
      allocate(ts_coordinates_a(3,natoms))
      allocate(ts_coordinates2_a(3,natoms))
      allocate(ts_coordinates(3*natoms-6))
      allocate(ts_coordinates2(3*natoms-6))
      open(unit=33,file=filets,status='old')
      call next_geo(ts_coordinates_a,natoms,33,has_next)
      call xyz_2int(ts_coordinates_a,ts_coordinates,nat)
      close(33)
      open(unit=34,file=filets2,status='old')
      call next_geo(ts_coordinates2_a,natoms,34,has_next)
      call xyz_2int(ts_coordinates2_a,ts_coordinates2,nat)
      close(34)
      if (rank .eq. 0) then
         write(*,*) "You´re using the dQ-couplingterm with the ",off_basis,"-basis"
      end if
   end if
      
!
!   Do the levenberg-marquardt optimization: In the 3x3 case, all settings
!   are handled by a single subroutine.
!

    
   allocate(energies_result(reference_counter))
   if (rank .eq. 0) then
      call optimize_3evb(energies_qmdff,reference_counter,&
              &energies_result,fileenergy,filegeo)
   end if
! 
!   Write the path-energies of the optimized parameters
!
   if (rank .eq. 0) then
      open(unit=17,file="energies_final.dat",status="unknown")
      do i=1,reference_counter
         write(17,*) energies_result(i)
      end do
      close(17)
   end if
end if
!
!     perform any final tasks before program exit
!
!     remove artifact file
!call system("rm fort.10");
!
!    calculate the needed time for optimization
!

if (rank .eq. 0) then
   call system_clock(time_int(2))

   do i=1,10
      time(i)=real(time_int(i))/1000
   end do


   duration=time(2)-time(1)

   write(*,'(A, F12.3, A)') " The calculation needed a time of",duration," seconds."
   write(15,'(A, F12.3, A)') " The calculation needed a time of",duration," seconds."
   write(15,*)
end if
deallocate(coord)
close(15)
end program evbopt

