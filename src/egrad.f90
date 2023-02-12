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
!  for Wilson matrix printout
real(kind=8),allocatable::internal(:),B_mat(:,:),dB_mat(:,:,:)
logical::has_next,path_struc,print_wilson
integer::wilson_mode
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
call read_pes(rank)
!
!    Read in the structures of a reactionpath
!
filegeo="xyasdfdud"
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:16) .eq. 'COORDS_FILE ') then
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
if (filegeo .eq. "xyasdfdud") then
   write(*,*) "Please add the COORDS_FILE keyword!"
   call fatal
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
!     If, for coordinate analysis etc., the Wilson matrix for each structure shall
!     be written to file 
!
print_wilson=.false.
do i = 1, nkey
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
!     Read in the internal coordinates from file if the Wilson matrix 
!     shall be calculated for each structure
!
if (print_wilson) then
   call init_int("dummy",0,0,1)
   allocate(internal(nat6),B_mat(nat6,3*natoms))
   open(unit=215,file="wilson_mat.dat",status="replace")
   write(215,*) "# No. int. coord.(i)   No. cart. coord(j)        B(i,j)"
   if (wilson_mode .eq. 2) then
      allocate(dB_mat(nat6,3*natoms,3*natoms))
      open(unit=216,file="wilson_deriv.dat",status="replace")
      write(216,*) "# No. int. coord.(i)   No. cart. coord(j)    No.cart.coord(k)     B'(i,j,k)"
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
open(unit=166,file=filegeo,status='old',iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file ",filegeo," with the geometries to be calculated is missing!"
   write(*,*) "Please edit the COORDS_FILE keyword!"
   call fatal
end if

open(unit=44,file="energies.dat",status='unknown')
open(unit=45,file="gradients.dat",status='unknown')
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
if (treq) open (unit=48,file="treq.out",status="unknown")
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
!     If desired, calculate the Wilson matrix the structure and print it to file
!
   if (print_wilson) then
      write(215,*) "Structure",ref_count,":"
      call xyz_2int(coord,internal,natoms)
      call calc_wilson(coord,internal,B_mat)
      do i=1,nat6
         do j=1,3*n_one
            write(215,*) i,"      ",j,"      ",B_mat(i,j)
         end do
      end do
!
!     If ordered, calculate and print also the Wilson matrix derivative!
!
      if (wilson_mode .eq. 2) then
         call calc_dwilson(coord,internal,dB_mat)
         do i=1,nat6
            do j=1,3*n_one
               do k=1,3*n_one
                  write(216,*) i,"      ",j,"      ",k,"      ",dB_mat(i,j,k)
               end do
            end do
         end do
      end if

   end if
!
!     Calculate and write the energies of the single qmdffs!
!
   if ((.not. treq) .and. (.not. pot_ana) .and. (nqmdff .gt. 1) ) then
      call eqmdff(coord,e_qmdff1,e_qmdff2)
      write(99,*) e_qmdff1,e_qmdff2
   end if
   write(44,*) e_evb
   write(45,*) "Structure",ref_count,":"
   do k=1,natoms
      write(45,*) g_evb(:,k)
   end do
 !  stop "test run, gradient loop" 
end do
if (print_wilson) then
   close(215)
   if (wilson_mode .eq. 2) close(216)
end if

if (treq) close(48)
close(172)
!
!    For debugging, print the QMDFF components to two files (one for each QMDFF)
!
close(166)
close(192)
close(44)
close(45)
close(99)
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

call cpu_time(time2)
duration=time2-time1

write(*,*)  ".. .----. -- -.. --- -. . "
write(*,*) "Calculation successfully finished!"
write(*,*) "Energies: energies.dat, Gradients: gradients.dat"
write(*,*) 
write(*,'(A, F10.3, A)') " The calculation needed a time of ",duration," seconds."
!
!     perform any final tasks before program exit
!

deallocate(coord)
end

