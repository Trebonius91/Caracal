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
! for manual internal coordinates
integer,dimension(:,:),allocatable::coord_tmp
integer,dimension(:),allocatable::coord_types
character(len=40)::coord_line
!
integer mode,next,j,k,readstatus,dg_evb_mode,mat_size
integer::int_mode  ! method for defining internal coordinates (if used)
logical::exist,exists,has_next,coupl1,par_soschl

logical::grad,frequency,opt_min,orca_fake,ts_opt,calc_irc
! for EVB-QMDFF hessian calculation
real(kind=8),dimension(:,:),allocatable::hess,mass_mat
real(kind=8),dimension(:),allocatable::h_out,freqs
integer::maxline
character(len=40)::commarg ! string for command line argument
! for fragment calculations
integer::indi_tmp(100)
character(len=80) a80
integer::maxcycle,rest
!     the MPI rank (here always 0)
integer::rank
!
!     Set MPI rank to zero for this program
!
rank=0
!
!     no MPI is used
!
use_mpi=.false.
!
!     no RPMDrate is used
!
use_rpmdrate = 0
!
!     set up the structure and mechanics calculation
!
call initial(rank)
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
write(15,*) "--- EXPLORE CALCULATION ON SINGLE STRUCTURE ---"
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
call getxyz
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
frequency=.false.  
opt_min=.false.
ts_opt=.false.
calc_irc=.false.
! 
!    the default values
!
maxiter=500
stepmax=0.05d0
ethr=1D-7
gthr=1D-7
dthr=1E-7
gmaxthr=1D-7
dmaxthr=1D-7
newton_raphson=.false.
irc_maxstep=200
irc_steplen=0.1d0
irc_eulerlen=0.005d0
irc_ethr=1D-7
irc_gthr=1D-6

!
!     Which job shall be done?
!
jobtype="NONE"
do i = 1, nkey
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
!
!    Read in sections for different job types 
!
!    A: The OPT_MIN and OPT_TS jobs
!

if (trim(jobtype) .eq. "OPT_MIN" .or. trim(jobtype) .eq. "OPT_TS" &
            &  .or. trim(jobtype) .eq. "OPTFREQ") then
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'OPT {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'OPT{') then
         do j=1,nkey-i+1
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
!     maximum length of a single optimization step
            else if (keyword(1:20) .eq. 'STEPSIZE ') then
               read(record,*,iostat=readstat) names,stepmax
               if (readstat .ne. 0) then
                  write(*,*) "Correct format STEPSIZE [Opt. stepsize (A)]"
                  call fatal
               end if
!     energy change convergence criterion
            else if (keyword(1:20) .eq. 'ETHR ') then
               read(record,*) names,ethr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ETHR [energy (Hartree)]"
                  call fatal
               end if
!     gradient norm convergence criterion
            else if (keyword(1:20) .eq. 'GTHR ') then
               read(record,*) names,gthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format ETHR [gradient comp.]"
                  call fatal
               end if
!     geometry step norm convergence criterion
            else if (keyword(1:20) .eq. 'DTHR ') then
               read(record,*) names,dthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format DTHR [geom. step]"
                  call fatal
               end if
!     largest component in gradient vector
            else if (keyword(1:20) .eq. 'GMAXTHR ') then
               read(record,*) names,gmaxthr
               if (readstat .ne. 0) then
                  write(*,*) "Correct format GMAXTHR [gradient comp.]"
                  call fatal
               end if
!     largest component in geometry change vector
            else if (keyword(1:20) .eq. 'DMAXTHR ') then
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
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'IRC {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'IRC{') then
         do j=1,nkey-i+1
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
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'FREQ {' .or. trim(adjustl(record(1:11))) &
                  &  .eq. 'FREQ{') then
         do j=1,nkey-i+1
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


call cpu_time(time1)
!
!    Open file with start structure and read it in
!
allocate(coord(3,natoms))
write(15,*) "Read in the structure file: ",xyzfile
open (unit=98,file=xyzfile,status="unknown")
   call next_geo(coord,natoms,98,has_next)               
close(98)
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
   call gradient(coord,e_evb,g_evb,1)
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
   call g98fake("explore.molden",natoms,atind,coord,freqs,eigvecs,eigvecs)
   write(15,*) "File explore.molden with normal mode spectum written."

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
      call gradient(coord,e_evb,g_evb,1)
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


call cpu_time(time2)
duration=time2-time1

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

