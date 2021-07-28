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
!     ##  program egrad  --  (EVB-)QMDFF energy and gradient         ##
!     ##                                                             ##
!     #################################################################
!
!
!     "egrad" calculates for a series of concatenated xyz-Structures 
!     gradient and energy for a single QMDFF or 2/3 QMDFF´s with
!     EVB-coupling
!     

program egrad
use general
use evb_mod
use debug

implicit none
integer::nat,nat3,ref_count
integer::readstat  ! for error handling
real(kind=8)::e_evb
real(kind=8),dimension(:,:),allocatable::coord
real(kind=8),dimension(:,:),allocatable::g_evb
real(kind=8)::e_qmdff1,e_qmdff2 ! energies of the single QMDFFs
character(len=70)::filegeo,names
character(len=120)::string
character(len=20)::keyword
character(len=120)::record
integer::i,j,k
real(kind=8)::s_var,z_var   ! variables for RP-EVB
logical::has_next,path_struc
logical::exist
integer::next
character(len=40)::commarg ! string for command line argument
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
      call help("egrad")
      stop
   else
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'egrad.x -help' or 'egrad.x -h'."
   end if
end if

!
!     Read keywords from the key-file
!
call getkey(rank)
!
!     Read in the QMDFF and EVB terms
!
call read_evb(rank)
!
!    Read in the structures of a reactionpath
!
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'PATH_STRUCTURE ') then
      read(record,*) names,filegeo
      inquire(file=filegeo,exist=path_struc)
   end if
end do
if (.not. path_struc) then
   exist=.false.
   do while (.not. exist)
      write(iout,'(/," Name of the xyz-File with Structures:  ",$)')
      read (*,'(a120)')  filegeo
      inquire(file=filegeo,exist=exist)
   end do
end if

!
!     If the debugging mode shall be started to print out more details!
!
do_debug=.false.
do i = 1, nkey
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
!     Print debug information message
!     Also allocate debug arrays
!
if (do_debug) then
   write(*,*) "THE DEBUGGING MODE WAS REQUESTED; ADDITIONAL INFOS WILL & 
       & BE PRINTED"

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
!     produce the EVB-QMDFF-energies for all points
!
ref_count=0
open(unit=166,file=filegeo,status='old')
open(unit=44,file="energies.qmdff",status='unknown')
open(unit=45,file="gradients.qmdff",status='unknown')
write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
write(*,*) "Calculating energy and gradient ..."
allocate(g_evb(3,natoms))
allocate(coord(3,natoms))
call cpu_time(time1)  ! measure the needed time for comparization
!
!     open file with single QMDFF energies
!
open (unit=99,file="single_qmdff.dat",status='unknown')  
open(unit=172,file="grad_square.dat",status="unknown")  ! TEST for gradients
if (rp_evb) open (unit=48,file="rp_evb.out",status="unknown")
if (int_grad_plot) open (unit=192,file="int_grad.out",status="unknown")
do
   ref_count=ref_count+1
   call next_geo(coord,natoms,166,has_next)
   if (.not.has_next) exit
! 
!     Distunguish between debug and normal calculation
!
   coord=coord/bohr
   call gradient(coord,e_evb,g_evb,1)  ! else calculate the usual gradient
!
!     Calculate and write the energies of the single qmdffs!
!
   if ((.not. rp_evb) .and. (.not. pot_ana) .and. (nqmdff .gt. 1) ) then
      call eqmdff(coord,e_qmdff1,e_qmdff2)
      write(99,*) e_qmdff1,e_qmdff2
   end if
   if (mueller_brown) then
      write(44,*) coord(1,1),coord(2,1),e_evb
      if ((abs(g_evb(1,1)/100000.d0) .le. 0.5) .and. (abs(g_evb(2,1)/100000.d0).le.0.5)) then
         write(45,*) coord(1,1),coord(2,1),g_evb(1,1)/100000.d0,g_evb(2,1)/100000.d0
      else 
         write(45,*) coord(1,1),coord(2,1),0.d0,0.d0
      end if
   else 
      write(44,*) e_evb
      write(45,*) "Structure",ref_count,":"
      do k=1,natoms
         write(45,*) g_evb(:,k)
      end do
   end if
 !  stop "test run, gradient loop" 
end do
if (rp_evb) close(48)
close(172)
!
!    For debugging, print the QMDFF components to two files (one for each QMDFF)
!
close(166)
close(192)
close(44)
close(45)
close(99)
if (rp_evb) close(48)
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

call cpu_time(time2)
duration=time2-time1

write(*,*)  ".. .----. -- -.. --- -. . "
write(*,*) "Calculation successfully finished!"
write(*,*) "Energies: energies.qmdff, Gradients: gradients.qmdff"
write(*,*) 
write(*,'(A, F10.3, A)') " The calculation needed a time of ",duration," seconds."
!
!     perform any final tasks before program exit
!

deallocate(coord)
end

