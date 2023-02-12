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
!     #################################################################
!     ##                                                             ##
!     ##  program irc  --  TS and reaction path optimization         ##
!     ##                                                             ##
!     #################################################################
!
!     "irc" does a transition state optimization followed by an IRC
!     optimization for one of the methods availiable with the whole 
!     EVB-QMDFF program package (EVB-couplingterms as well as analytical
!     potential functions)
!     As initial point, an approximative TS structure will be needed 
!

program irc
use general
use evb_mod

implicit none 

!     the MPI rank (here always 0)
integer::rank
integer::i,j  ! loop indices
!     for keyword read in
character(len=20)::keyword
character(len=60)::names
character(len=120)::record
character(len=120)::string
integer::next
integer::int_mode  ! method for defining internal coordinates
character(len=40)::commarg ! string for command line argument
!
!     Set MPI rank to zero for this program
!
rank=0
!
!     no MPI is used
!
use_mpi=.false.
!
!     set up the structure and molecular mechanics calculation
!
call initial(rank)
!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("irc")
      stop
   else
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'irc.x -help' or 'irc.x -h'."
   end if
end if

!
!     Read keywords from the key-file
!
call getkey(rank)
!
!     Print info message for actual program specification
!
if (rank .eq. 0) then
   write(*,*) "PROGRAM IRC:"
   write(*,*) "Starting from an approximated TS structure, optimize the TS"
   write(*,*) "and upfollowing the IRC (minimum energy path) for the reaction."
   write(*,*) 
end if
!
!     Open logfile for TS and IRC optimization
!     print header and time 
!
open (unit=15,file="irc.log",status="unknown")
write(15,*) "--- EVB-QMDFF CALCULATION FOr TS/IRC-OPTIMIZATION ---"
write(15,*)
write(15,*) "Calculation initiated at: "
call timestamp ( )


call cpu_time(time1)

!
!     Read in needed parameters for IRC optimization
! 
!    the default values
!
maxiter=200
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
!     default method so far..
irc_method="euler"
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
!
!     maximum number of optimization steps
!
   if (keyword(1:20) .eq. 'OPT_MAXITER ') then
      read(record,*) names,maxiter
!
!     maximum length of a single optimization step
!
   else if (keyword(1:20) .eq. 'OPT_STEPSIZE ') then
      read(record,*) names,stepmax
!
!     energy change convergence criterion
!
   else if (keyword(1:20) .eq. 'OPT_ETHR ') then
      read(record,*) names,ethr
!
!     gradient norm convergence criterion
!
   else if (keyword(1:20) .eq. 'OPT_GTHR ') then
      read(record,*) names,gthr
!
!     geometry step norm convergence criterion
!
   else if (keyword(1:20) .eq. 'OPT_DTHR ') then
      read(record,*) names,dthr
!
!     largest component in gradient vector
!
   else if (keyword(1:20) .eq. 'OPT_GMAXTHR ') then
      read(record,*) names,gmaxthr
! 
!     largest component in geometry change vector
!     
   else if (keyword(1:20) .eq. 'OPT_DMAXTHR ') then
      read(record,*) names,dmaxthr
! 
!     If newton raphson shall be used (else: the better P-RFO method)
!     
   else if (keyword(1:20) .eq. 'NEWTON_RAPHSON ') then
      newton_raphson=.true.
      write(*,*) "The simple Newton-Raphson method for PES searching will be used."
      write(15,*) "The simple Newton-Raphson method for PES searching will be used."
!
!     maximum number of IRC steps in each direction
!
   else if (keyword(1:20) .eq. 'IRC_MAXSTEP ') then
      read(record,*) names,irc_maxstep
!
!     Desired length of an IRC step ((amu)^1/2*bohr)
!
   else if (keyword(1:20) .eq. 'IRC_STEPLEN ') then
      read(record,*) names,irc_steplen
!
!     Desired length of an elementary Euler step ((amu)^1/2*bohr)
!     --> should be much smaller than the IRC step itself!
!
   else if (keyword(1:20) .eq. 'IRC_EULERLEN ') then
      read(record,*) names,irc_eulerlen
!
!     IRC energy change convergence criterion
!
   else if (keyword(1:20) .eq. 'IRC_ETHR ') then
      read(record,*) names,irc_ethr
!
!     IRC gradient norm convergence criterion
!
   else if (keyword(1:20) .eq. 'IRC_GTHR ') then
      read(record,*) names,irc_gthr
   end if

end do
write(15,*)
write(15,*) "------------------SETTINGS------------------"
write(15,'(a)') " (A) TS-OPTIMIZATION:"
write(15,'(a,i10)') " * Maximum number of optimization steps: ",maxiter
write(15,'(a,e14.5)') " * Maximum length of a single geometry step: ",stepmax
write(15,'(a,e14.5)') " * Energy change convergence criterion (ETHR): ",ethr
write(15,'(a,e14.5)') " * Gradient norm convergence criterion (GTHR): ",gthr
write(15,'(a,e14.5)') " * Step norm convergence criterion (DTHR): ",dthr
write(15,'(a,e14.5)') " * Largest gradient component for convergence (GMAXTHR): ",gmaxthr
write(15,'(a,e14.5)') " * Largest step component for convergence (DMAXTHR): ",dmaxthr
write(15,'(a)') " (B) IRC-OPTIMIZATION:"
write(15,'(a,i10)') " * Maximum number of IRC steps per side: ",irc_maxstep
write(15,'(a,e14.5)') " * Length of a single IRC step ((amu)^1/2*bohr): ",irc_steplen
write(15,'(a,e14.5)') " * Length of a single Euler IRC step ((amu)^1/2*bohr): ",irc_eulerlen
write(15,'(a,e14.5)') " * IRC energy change convergence criterion (IRC_ETHR): ",irc_ethr
write(15,'(a,e14.5)') " * IRC gradient norm convergence criterion (IRC_GTHR): ",irc_gthr
write(15,*) "--------------------------------------------"
write(15,*)

!
!     Read in the approximated TS structure!
!
call getxyz
!
!     Read in the QMDFF and EVB terms
!
call read_pes(rank)
!
!     Do all calculations in internal coordinates!
!     Specify them automatically and the full set of 3N-6 coordinates
!     unless read in of internal coordinates is ordered explicitly
!
if (.not. dg_evb .and. .not. treq) then
   read_coord=.false.
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:16) .eq. 'READ_COORD ') then
         read_coord=.true.
      end if
   end do
   int_mode=3
   if (read_coord) int_mode=1
   call init_int("ts.xyz",1,1,int_mode)
end if
!
!     FIRST PART: TS optimization: call the respective subroutine
!
if (rank .eq. 0) then
   call opt_ts
end if
!
!     SECOND PART: IRC optimization: call the respective subroutine
!
if (rank .eq. 0) then
   call opt_irc
end if

call cpu_time(time2)
duration=time2-time1

write(*,*)  ".. .----. -- -.. --- -. . "
write(*,*) "Calculation successfully finished!"
write(*,*) "Output was written to evb_qmdff.log"
write(*,'(A, F10.3, A)') " The calculation needed a time of",duration," seconds."

end program irc
