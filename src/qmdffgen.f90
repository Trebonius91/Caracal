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
!     ##  program qmdffgen  --  Automatically generate QMDFF´s       ##
!     ##                                                             ##
!     #################################################################
!
!
!     "qmdffgen" simiplifies the generation of a QMDFF from given 
!     orca/gaussian-output files as references. 
!     Two or three QMDFF´s are generated directly. 
!     From the orca-calculation the orca.hess and orca.out
!     files (with the same prefix) needs to be present.
!     Alternatively, the gaussian.chk file needs to be present
!     

program qmdffgen
use general ! general parameters
use qmdff  ! module with qmdff global variables
use evb_mod ! evb coupling parameters

implicit none
integer::next,freeunit,i,qmdffnumber,asciinum
integer::length,readstat
integer::iz1_two,iz1_three,iz2_two,iz2_three,ndummy
integer::parameters
real(kind=8)::dens_two,dens_three,c6_two,c6_three
real(kind=8),allocatable::xyz2(:,:)
real(kind=8)::e1_shifted,e2_shifted,e3_shifted,e1,e2,e3
real(kind=8)::e1_ref,e2_ref,e3_ref
real(kind=8)::vz(94)
real(kind=8)::r42,c6,dens
real(kind=8)::en_gauss1,en_gauss2 ! for read in of gaussian calculation result
integer::j,i1,i2,iz1,iz2,lin,n2,l1,l2,l3,k
character(len=80)::fname,names
real(kind=8)::scalehb(94),scalexb(94)
character(len=20)::keyword
character(len=1)::qmdffnum
character(len=120)::record
character(len=120)::string
character(len=80)::prefix,line
character(len=:),allocatable::prefix1,prefix2,prefix3
character(len=80)::pre1,pre2,pre3
character(len=80)::fffilen1,fffilen2,fffilen3,a80,header
character(len=85)::textout,textout2,buffer(4)
character(len=100)::sys_line ! for system call in case of gaussian reference
character(len=40)::commarg ! string for command line argument
real(kind=8)::en_new1
real(kind=8)::en_new2
real(kind=8),allocatable::internal1(:),wilson(:,:),dwilson(:,:,:) ! for Wilson matrix calculation
logical::prefix_key,exist
integer::fchkstat ! for gaussian output: status of fchk output conversion!    
integer::rank
logical::extra_ens  ! read in the energies from another level of theory!
!
!     Set MPI rank to zero for this program
!
rank=0

!
!     no MPI is used
!
use_mpi=.false.

!
!     The explore program is not used
!
use_explore = .false.
 
!
!     Show the nice promo-banner at the beginning!
!
call promo
!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("qmdffgen")
      stop
   else 
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'qmdffgen.x -help' or 'qmdffgen.x -h'."
   end if
   details=.false.
   if (trim(commarg) .eq. "-details") then
      details=.true.
   end if

end if

!
!     set the integer to tell the subroutine that qmdffgen is executed
!
!     Check if command line arguments were submitted
!
i = 0
do
  call get_command_argument(i, commarg)
  if (len_trim(commarg) == 0) exit
  i = i+1
end do

if (i .lt. 2) then
   whichprog = 1
else 
   whichprog = 2
end if

prefix_key = .false.
read_software = .false.
!
!     if VASP has been used with selective dynamics
!
vasp_hessian_sel = .false.
!
!     Read keywords from input key-file 
!
call getkey(rank)

!
!     Read the reference program/software from input file
!

do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'SOFTWARE ') then
      qmdffnumber=1
      read(record,*) prefix,software
!     possible options: O (orca), G (gaussian), T (turbomole), C (CP2K), V(VASP)
      read_software=.true.
   end if
end do
!
!     If the QMDFF names are stored already in a file
!
check_coord=.false.
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. '1QMDFF ') then
      qmdffnumber=1
      read(record,*) prefix,pre1
      allocate(character(len=LEN(TRIM(pre1))) :: prefix1)
      prefix1=trim(pre1)
      l1=LEN(TRIM(pre1))
      prefix_key=.true.
   else if (keyword(1:11) .eq. '2QMDFF ') then
      qmdffnumber=2
      read(record,*) prefix,pre1,pre2
      allocate(character(len=LEN(TRIM(pre1))) :: prefix1)
      allocate(character(len=LEN(TRIM(pre2))) :: prefix2)
      prefix1=trim(pre1)
      if (software .eq. "G") then
         textout = prefix1 // ".log"
      else
         textout = prefix1 // ".out"
      end if
      prefix2=trim(pre2)
      if (software .eq. "G") then
         textout = prefix2 // ".log"
      else
         textout = prefix2 // ".out"
      end if
      l1=LEN(TRIM(pre1))
      l2=LEN(TRIM(pre2))
      prefix_key=.true.
   else if (keyword(1:11) .eq. '3QMDFF ') then
      qmdffnumber=3
      read(record,*) prefix,pre1,pre2,pre3
      allocate(character(len=LEN(TRIM(pre1))) :: prefix1)
      allocate(character(len=LEN(TRIM(pre2))) :: prefix2)
      allocate(character(len=LEN(TRIM(pre3))) :: prefix3)
      prefix1=trim(pre1)
      prefix2=trim(pre2)
      prefix3=trim(pre3)
      l1=LEN(TRIM(pre1))
      l2=LEN(TRIM(pre2))
      l3=LEN(TRIM(pre3))
      prefix_key=.true.
!
!     No QMDFF will be generated, instead, the internal coordinates 
!      of the given structure will be built and the Wilson matrix 
!      and its derivative will be generated
!
   else if (keyword(1:11) .eq. 'CHECK_COORD ') then
      check_coord=.true.
      write(*,*) 
      write(*,*) "The CHECK_COORD option was activated! No Hessian will be "
      write(*,*) "  read in and no QMDFF will be generated, instead, the internal"
      write(*,*) "  coordinates of the structure will be obtained from the Wiberg-"
      write(*,*) "  Mayer bond orders and the Wilson matrix as well as its derivative"
      write(*,*) "  will be written."
   end if
end do
if (.not. read_software) then
   exist=.false.
   software="X"
   do while (.not. exist)
      write(iout,'(/," Used Software for reference calculations:",/, &
        &  " (O = orca, G = gaussian, C = CP2K, V = VASP): ",$)')
      read (*,'(A1)')  software
      if ((software .eq. "O") .or. (software .eq. "G") .or. &
         &  (software .eq. "T") .or. (software .eq. "C") .or. &
         &  (software .eq. "V")) then
         exist = .true.
         select case (software)
            case("O")
               write(*,*) "Reference is read in from ORCA output"
            case("G")
               write(*,*) "Reference is read in from GAUSSIAN output"
            case("T")
               write(*,*) "Reference is read in from TURBOMOLE output"
               write(*,*) "This is currently not supported!"
               call fatal
            case("C")
               write(*,*) "Reference is read in from CP2K output"
            case("V")
               write(*,*) "Reference is read in from VASP output"
            case default
               write(*,*) "None of the supported software was used! Try again!"
               exist=.false.
         end select
      end if
   end do
end if
!
!     OPTIONAL: if the energies were calculated with a higher level 
!      of theory than the other reference data, they can be read in if a file 
!      min_energies.dat with these values is present in the folder
!
write(*,*)
inquire(file="min_energies.dat", exist=extra_ens)
if (extra_ens) then
   write(*,*) "The file min_energies.dat was found. The reference energies will"
   write(*,*) " be read in from this file. If this behavior isn´t desired, remove"
   write(*,*) " that file!"
   open(unit=45,file="min_energies.dat",status="old")
   read(45,*) en_new1
   read(45,*) en_new2
   close(45)
else
   write(*,*) "Hint: if you want to use energies from a different level of theory"
   write(*,*) " generate the file min_energies.dat and store the energies into it."
end if

if (.not. prefix_key) then
   exist=.false.
   do while (.not. exist)
      write(iout,'(/," Number of QMDFF´s to generate: ",$)')
      read (*,*,iostat=readstat)  qmdffnumber
      if (readstat .ne. 0) cycle
!      asciinum = ICHAR(qmdffnum)
!
!     Convert vom ASCII to the real number (1=49,2=50,3=51)
!
!      select case (asciinum)
!      case (49)
!         qmdffnumber=1
!      case (50)
!         qmdffnumber=2
!      case (51)
!         qmdffnumber=3
!      end select
  !    write(*,*) "num",qmdffnumber
      if (qmdffnumber.eq.1 .or. qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
         exist=.true.
      end if 
   end do
!
!     Get the names of input files
!     --> distinguish between the different software packages 
!
   write(*,*) qmdffnum
   exist=.false.
   do while (.not. exist)
      if (.not. check_coord) then
         if (software .eq. "O") then
            write(iout,'(/, " Prefix of the first reference (name.hess &
                 and name.out must be in this folder!): ",$)')
         else if (software .eq. "G") then
            write(iout,'(/, " Prefix of the first reference (name.out &
                 and name.chk must be in this folder!): ",$)')
         else if (software .eq. "C") then
            write(iout,'(/, " Prefix of the first reference (name.out &
                 must be in this folder!): ",$)')
         else if (software .eq. "V") then
            write(iout,'(/, " Prefix of the first reference (name.OUTCAR &
                 and name.charges must be in this folder!): ",$)')
         end if
      else 
         if (software .eq. "O") then
            write(iout,'(/, " Prefix of the first reference (name.out &
                 must be in this folder!): ",$)')
         else if (software .eq. "G") then
            write(iout,'(/, " Prefix of the first reference (name.out &
                 and name.chk must be in this folder!): ",$)')
         else if (software .eq. "C") then
            write(iout,'(/, " Prefix of the first reference (name.out &
                 must be in this folder!): ",$)')
         else if (software .eq. "V") then
            write(iout,'(/, " Prefix of the first reference (name.OUTCAR &
                 and name.charges must be in this folder!): ",$)')
         end if
      end if
      read (*,'(A80)')  line
      allocate(character(len=LEN(TRIM(line))) :: prefix1)
      prefix1=trim(line)
      l1=LEN(TRIM(prefix1))
      if (software .eq. "G") then
         textout = prefix1 // ".log"
      else if (software .eq. "V") then
         textout = prefix1 // ".OUTCAR"
         textout2 = prefix1 // ".charges"
      else 
         textout = prefix1 // ".out"
      end if
      inquire(file=textout,exist=exist)
      if (exist) then
         inquire(file=textout2,exist=exist)
      end if
      if (.not. exist) then 
         deallocate(prefix1)
      end if
   end do
   if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,'(/," Prefix of the second reference: ",$)')
         read (*,'(A80)')  line
         allocate(character(len=LEN(TRIM(line))) :: prefix2)
         prefix2=trim(line)
         l2=LEN(TRIM(prefix2))
         if (software .eq. "G") then
            textout = prefix2 // ".log"
         else if (software .eq. "V") then
            textout = prefix1 // ".OUTCAR"
         else
            textout = prefix2 // ".out"
         end if
         inquire(file=textout,exist=exist)
         if (prefix1 .eq. prefix2) then
            exist = .false.
         end if
         if (.not. exist) then
            deallocate(prefix2)
         end if
      end do
   end if
   if (qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,'(/," Prefix of the third reference: ",$)')
         read (*,'(A80)')  line
         allocate(character(len=LEN(TRIM(line))) :: prefix3)
         prefix3=trim(line)
         l3=LEN(TRIM(prefix3))
         textout = prefix3 // ".out"
         inquire(file=textout,exist=exist)
         if (prefix1 .eq. prefix2 .or. prefix1.eq.prefix3) then
            exist=.false.
         end if
         if (.not. exist) then
            deallocate(prefix3)
         end if
      end do
   end if
end if
!
!     Read the reference energies of the single minima..
!
if (qmdffnumber .eq. 2) then
!
!     for orca output (eventually generated from gaussian output..)
!
   if (software .eq. "O") then
      open(unit=19,file=prefix1 // ".hess")
   
      read(19,'(a)')a80
      if(index(a80,'gauss2orca.pl').ne.0) then
         do i=1,6
            read(19,*) a80
         end do
         read(19,*) en_gauss1
      else 
         goto 52
      end if
  
      open(unit=20,file=prefix2 // ".hess") 
         do i=1,7
            read(20,*) a80
         end do
         read(20,*) en_gauss2
      close(20)
!
!   Write the energies in the last lines of both *.out files
!
      write(*,*)
      write(*,*) "I have noticed that your hessian matrix originates from"
      write(*,*) "a Gaussian calculation."
      write(*,*) "Therefore the gaussian energies were written at" 
      write(*,*) "the end of the *.out files." 
      write(*,*)
      write(sys_line,*) "echo 'FINAL SINGLE POINT ENERGY", en_gauss1,"' >> min1.out"
      call system(sys_line)
      write(sys_line,*) "echo 'FINAL SINGLE POINT ENERGY", en_gauss2,"' >> min2.out"
      call system(sys_line)
      52 continue
      close(19)
   end if
end if
! 
!     produce the QMDFFs 
!
!
!    Measure the needed time for production
!
call cpu_time(time1)
!
!     QMDFF1
!
write(*,*) "First QMDFF:  ",prefix1,".qmdff.."
write(*,*) "Logfile written to ",prefix1,"_qmdff.log"
write(*,*) "Calculating......."

length=len(trim(prefix1))
open(unit=10,file=prefix1//"_qmdff.log",status="unknown")

second_qmdff = .false.
call main_gen(prefix1,length)
fffilen1= prefix1 // ".qmdff"
call rdsolvff0(n_one,fffilen1)
natoms=n_one

allocate(at(n_one),xyz(3,n_one),at2(n_one),g_one(3,n_one),q(n_one),&
        c6xy(n_one,n_one),cn(n_one),imass(n_one),xyz2(3,n_one))
call copyc6 ! not molecule/FF specific
call setnonb(scalehb,scalexb,vz,sr42,zab,r0ab)
r094_mod=r094


call rdsolvff(n_one,xyz,at,q,imass,dens,scalehb,scalexb,fffilen1)
call ncoord_qmdff(n_one,rcov,at,xyz,cn,5000.0d0)
do i1=1,n_one
   iz1=at(i1)
   do i2=1,i1
      iz2=at(i2)
      call getc6(5,94,c6ab,maxci,iz1,iz2,cn(i1),cn(i2),c6)
      c6xy(i2,i1)=c6
      c6xy(i1,i2)=c6
   end do
end do

call ff_eg(n_one,at,xyz,e1,g_one)
call ff_nonb(n_one,at,xyz,q,r0ab,zab,r094_mod,sr42,c6xy,e1,g_one)
call ff_hb(n_one,at,xyz,e1,g_one)


if (software .eq. "O") then
   open(unit=42,file=prefix1 // ".out")
100   read(42,'(a)',end=99)a80
      if(index(a80,'FINAL SINGLE POINT ENERGY').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),buffer(4),e1_ref
      end if
      goto 100
99    close(42)
!    e-shifted = e_referenz - e_qmdff
   e1_shifted=e1_ref-e1
   close(10)
!
!     For CP2K: problem: in a vibrational analysis calculation no energy
!     for the undistorted structure is calculated! Therefore the energy of the 
!     minimum is not totally true in this case
!
else if (software .eq. "C") then
   open(unit=42,file=prefix1 // ".out")
101   read(42,'(a)',end=98)a80
      if(index(a80,'VIB|                                     Total Energy:').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),e1_ref
         goto 98
      end if
      goto 101
98    close(42)
   e1_shifted=e1_ref-e1
   close(10)
!
!     For gaussian output: First generate the name.fchk file!
!     The formchk-program needs to be present!
!
else if (software .eq. "G") then
!
!     To avoid useless Gaussian segmentation faults, remove the old 
!     name.fchk file if existent
!
   exist=.false.
   inquire(file=prefix1//".fchk",exist=exist)
   if (exist) then
      call system("rm "//prefix1//".fchk")
      call system("sleep 0.2")
   end if
   call system("formchk "//prefix1//".chk "//prefix1//".fchk", fchkstat)
   if (fchkstat .ne. 0) then
      write(*,*) "The formchk command couldn´t be executed! Check if all is set correctly!"
      call fatal
   end if
   open(unit=42,file=prefix1 // ".fchk")
102   read(42,'(a)',end=97)a80
      if (index(a80,'Total Energy                               R').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),e1_ref
         goto 97
      end if
      goto 102
97    close(42)
   e1_shifted=e1_ref-e1
   close(10)
else if (software .eq. "V") then
!
!     For VASP output: Read in the energy of the first calculated geometry 
!     (the undistorted geometry) and convert it from eV to Hartrees!
!
   open(unit=42,file=prefix1 // ".OUTCAR")
133   read(42,'(a)',end=96)a80
      if(index(a80,'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)').ne.0) then
         read(42,*)
         read(42,*)
         read(42,*)
         read(42,*) buffer(1),buffer(2),buffer(3),buffer(1),buffer(2),buffer(3),e1_ref
         e1_ref=e1_ref/evolt
         goto 96
      end if
      goto 133
96    close(42)
   e1_shifted=e1_ref-e1
   write(*,*) "energy1",e1_ref,e1
   close(10)
else 
   stop "No read in for energies implemented for that method so far!"
end if
!
!     if the energies are given from a seperate level of theory
!
if (extra_ens) then
   e1_shifted=en_new1-e1
end if

!
!     For one QMDFF, print out the Wilson matrix of the reference structure in the
!     QMDFF internal coordinates
!
!     First, fill the coord_def array
!
nat6=nbond+nangl+ntors
allocate(coord_def(nat6,5))
allocate(internal1(nat6))
allocate(wilson(nat6,3*n_one))
allocate(dwilson(nat6,3*n_one,3*n_one))
do i=1,nbond
   coord_def(i,1)=1
   coord_def(i,2:3)=bond(1:2,i)
   coord_def(i,4:5)=0
end do
do i=1,nangl
   coord_def(i+nbond,1)=2
   coord_def(i+nbond,2:4)=angl(1:3,i)
   coord_def(i+nbond,5)=0
end do
do i=1,ntors
   coord_def(i+nbond+nangl,1)=3
   coord_def(i+nbond+nangl,2:5)=tors(1:4,i)
end do

!
!    Second, calculate the internal coordinate array    
!
call xyz_2int(xyz,internal1,n_one)
!
!    Do the coordinate analysis steps only for small systems (less than 20 atoms)!
!
!    Third, calculcate the Wilson matrix and its derivative
!
if (check_coord) then
   call calc_wilson(xyz,internal1,wilson)
   call calc_dwilson(xyz,internal1,dwilson)
!
!
!
   open(unit=45,file="coord_analysis.dat",status="replace")
   write(45,*) "# Coordinate analysis for the molecule ",prefix1,":"
   write(45,*) "# (generated by qmdffgen.x)"
   write(45,*)
   write(45,*) "# Coordinate definitions and values:"
   write(45,*) "# bonds (coord. No.,   Atom1,   Atom2,   value (bohr))"
   do i=1,nbond
      write(45,*) i,coord_def(i,2),coord_def(i,3),internal1(i)
   end do
   if (nangl > 0) then
      write(45,*) "# angles (coord. No.,   Atom1,   Atom2,   Atom3,  value (radians))"
      do i=1,nangl
         write(45,*) i+nbond,coord_def(i+nbond,2),coord_def(i+nbond,3), &
                 & coord_def(i+nbond,4),internal1(i+nbond)
      end do
   end if
   if (ntors > 0) then
      write(45,*) "# dihedrals (coord. No.,   Atom1,   Atom2,   Atom3,   Atom4,   value (radians))"
      do i=1,ntors
         write(45,*) i+nbond,coord_def(i+nbond+nangl,2),coord_def(i+nbond+nangl,3), &
             & coord_def(i+nbond+nangl,4), &
             & coord_def(i+nbond+nangl,5),internal1(i+nbond+nangl)
      end do
   end if
   write(45,*)
   write(45,*) "# Wilson matrix for coordinate transformation:"
   write(45,*) "# No. int. coord.(i)   No. cart. coord(j)        B(i,j)"
   do i=1,nat6
      do j=1,3*n_one
         write(45,*) i,"      ",j,"      ",wilson(i,j)
      end do
   end do

   write(45,*)
   write(45,*) "# Wilson matrix derivative for Hessian matrix transformation:"
   write(45,*) "# No. int. coord.(i)   No. cart. coord(j)     No.cart.coord(k)       B'(i,j,k)"
   do i=1,nat6
      do j=1,3*n_one
         do k=1,3*n_one
            write(45,*) i,"      ",j,"      ",k,"      ",dwilson(i,j,k)
         end do
      end do
   end do
   write(*,*) "Internal coordinate details (+Wilson matrix) written to 'coord_analysis.dat'"
end if

close(45)
!
!     Write the internal coordinates to the coord_def.inp file for further usage
!
open(unit=46,file="coord_def.inp",status="replace")

do i=1,nbond
   write(46,*) coord_def(i,2),coord_def(i,3)
end do
do i=1,nangl
   write(46,*) coord_def(i+nbond,2),coord_def(i+nbond,3),coord_def(i+nbond,4)
end do
do i=1,ntors
   write(46,*) coord_def(i+nbond+nangl,2),coord_def(i+nbond+nangl,3), &
             & coord_def(i+nbond+nangl,4),coord_def(i+nbond+nangl,5),0
end do

close(46)
write(*,*) "Internal coordinates input file 'coord_def.dat' written."
!
!     QMDFF2
!

if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
   write(*,*) "Second QMDFF:  ",prefix2,".qmdff.."
   write(*,*) "Logfile written to ",prefix2,"_qmdff.log"   
   write(*,*) "Calculating......."

   length=len(trim(prefix2))
   open(unit=10,file=prefix2//"_qmdff.log",status="unknown")
   second_qmdff=.true.
   call main_gen(prefix2,length)
   fffilen2=prefix2 // ".qmdff"
   call rdsolvff0(n_two,fffilen2)
   if (n_two .ne. n_one) then
      write(*,*)'qmdff1 has ',n_one,' atoms'
      write(*,*)'qmdff2 has ',n_two,' atoms'
      write(*,*)'   This will not work!'
      stop 'atom number mismatch between FF1 and FF2'
   end if
   natoms=n_one
   allocate(at_two(n_one),xyz_two(3,n_one),g_two(3,n_one),q_two(n_one),&
           c6xy_two(n_one,n_one),cn_two(n_one),imass_two(n_one))

   call rdsolvff_two(n_one,xyz_two,at_two,q_two,imass_two,dens_two,scalehb,scalexb,fffilen2)
   call ncoord_qmdff(n_one,rcov,at_two,xyz_two,cn_two,5000.0d0)
   do i1=1,n_one
      iz1_two=at_two(i1)
      do i2=1,i1
         iz2_two=at_two(i2)
         call getc6(5,94,c6ab,maxci,iz1_two,iz2_two,cn_two(i1),cn_two(i2),c6_two)
         c6xy_two(i2,i1)=c6_two
         c6xy_two(i1,i2)=c6_two
      end do
   end do
   call ff_eg_two(n_one,at_two,xyz_two,e2,g_two)
   call ff_nonb_two(n_one,at_two,xyz_two,q_two,r0ab,zab,r094_mod,sr42,c6xy_two,e2,g_two)
   call ff_hb_two(n_one,at_two,xyz_two,e2,g_two)
   if (software .eq. "O") then
      open(unit=43,file=prefix2 // ".out")
110   read(43,'(a)',end=89)a80
      if(index(a80,'FINAL SINGLE POINT ENERGY').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),buffer(4),e2_ref
      end if
      goto 110
89    close(43)
   else if (software .eq. "C") then
      open(unit=42,file=prefix1 // ".out")
111   read(42,'(a)',end=88)a80
      if(index(a80,'VIB|                                     Total Energy:').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),e2_ref
         goto 88
      end if
      goto 111
88    close(42)
     !
!     For gaussian output: First generate the name.fchk file!
!     The formchk-program needs to be present!
!
   else if (software .eq. "G") then
!
!     To avoid useless Gaussian segmentation faults, remove the old 
!     name.fchk file if existent
!
      exist=.false.
      inquire(file=prefix2//".fchk",exist=exist)
      if (exist) then
         call system("rm "//prefix2//".fchk")
         call system("sleep 0.2")
      end if
      call system("formchk "//prefix2//".chk "//prefix2//".fchk", fchkstat)
      if (fchkstat .ne. 0) then
         write(*,*) "The formchk command couldn´t be executed! Check if all is set correctly!"
         call fatal
      end if
      open(unit=42,file=prefix2 // ".fchk")
103   read(42,'(a)',end=199)a80
      if (index(a80,'Total Energy                               R').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),e2_ref
         goto 199
      end if
      goto 103
199    close(42)
   else if (software .eq. "V") then
!
!     For VASP output: Read in the energy of the first calculated geometry 
!     (the undistorted geometry) and convert it from eV to Hartrees!
!
      open(unit=42,file=prefix2 // ".OUTCAR")
104   read(42,'(a)',end=96)a80
      if(index(a80,'FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)').ne.0) then
         read(42,*)
         read(42,*)
         read(42,*)
         read(42,*) buffer(1),buffer(2),buffer(3),buffer(1),buffer(2),buffer(3),e2_ref
         e2_ref=e2_ref/evolt
         goto 204
      end if
      goto 104
204   close(42)
      e2_shifted=e2_ref-e2
      write(*,*) "energy2",e2_ref,e2
   end if
   close(10)
!
!     if the energies are given from a seperate level of theory
!
   if (extra_ens) then
      e2_shifted=en_new2-e2
   else 
      e2_shifted=e2_ref-e2
   end if
end if

!
!     QMDFF3
!

if (qmdffnumber.eq.3) then
   write(*,*) "Third QMDFF:  ",prefix3,".qmdff.."
   write(*,*) "Logfile written to ",prefix3,"_qmdff.log"
   write(*,*) "Calculating......."

   length=len(trim(prefix3))
   open(unit=10,file=prefix3//"_qmdff.log",status="unknown")

   call main_gen(prefix3,length)
   fffilen3=prefix3 // ".qmdff"
   call rdsolvff0(n_three,fffilen3)
   if (n_three .ne. n_one) then
      write(*,*)'qmdff1 has ',n_one,' atoms'
      write(*,*)'qmdff3 has ',n_three,' atoms'
      write(*,*)'   This will not work!'
      stop 'atom number mismatch between FF1 and FF3'
   end if
   natoms=n_one
   allocate(at_three(n_one),xyz_three(3,n_one),g_three(3,n_one),q_three(n_one),&
           c6xy_three(n_one,n_one),cn_three(n_one),imass_three(n_one))

   call rdsolvff_three(n_one,xyz_three,at_three,q_three,imass_three,dens_three,scalehb,&
         scalexb,fffilen3)
   call ncoord_qmdff(n_one,rcov,at_three,xyz_three,cn_three,5000.0d0)

   do i1=1,n_one
      iz1_three=at_three(i1)
      do i2=1,i1
         iz2_three=at_three(i2)
         call getc6(5,94,c6ab,maxci,iz1_three,iz2_three,cn_three(i1),&
               cn_three(i2),c6_three)
         c6xy_three(i2,i1)=c6_three
         c6xy_three(i1,i2)=c6_three
      end do
   end do

   call ff_eg_three(n_one,at,xyz_three,e3,g_three)
   call ff_nonb_three(n_one,at,xyz_three,q_three,r0ab,zab,r094_mod,sr42,c6xy_three,e3,g_three)
   call ff_hb_three(n_one,at,xyz_three,e3,g_three)
   open(unit=44,file=prefix2 // ".out")
120      read(44,'(a)',end=79)a80
         if(index(a80,'FINAL SINGLE POINT ENERGY').ne.0) then
         read(a80,*) buffer(1),buffer(2),buffer(3),buffer(4),e3_ref
         end if
         goto 120
79 close(44)
   close(10)
!    e-shifted = e_referenz - e_qmdff
   e3_shifted=e3_ref-e3
end if
!
!    Write the energies in a generic caracal.key file for EVB-optimization
!
open(unit=81,file="caracal.key",status='unknown')
write(81,'(a)')"# This is an automatically generated keyword file for CARACAL "
write(81,'(a)')"# EVB-QMDFF optimizations with the program evbopt.x"
write(81,*)
!
!    Which kind of EVB coupling is used: for 2 QMDFFs always de_evb as default
!
if (qmdffnumber .eq. 1) then
   write(81,*) "pes qmdff"
else if (qmdffnumber .eq. 2) then
   write(81,*) "pes de_evb"
end if
!
!    QMDFF information
!
write(81,*)
write(81,*) "qmdff {"
if (qmdffnumber.eq.1) then
   write(81,*)"   ffnames  ",prefix1 // ".qmdff"
   write(81,*)"   eshift",e1_shifted
else if (qmdffnumber.eq.2) then
   write(81,*)"   ffnames  ",prefix1 // ".qmdff  ",prefix2 // ".qmdff"
   write(81,*)"   eshift ",e1_shifted,e2_shifted
else if (qmdffnumber.eq.3) then
      write(81,*)"   ffnames  ",prefix1 // ".qmdff  ",prefix2 // ".qmdff  ",prefix3 // ".qmdff"
   write(81,*)"   eshift",e1_shifted,e2_shifted,e3_shifted
end if 
write(81,*) "}"
!
!    calculate the needed time for optimization
!
call cpu_time(time2)

duration=time2-time1
write(*,*)
write(*,'(A, F12.3, A)') " The calculation needed a time of",duration," seconds."
write(*,*)

write(*,*) "Calculations finished! See in both logfiles for details."
close(81)
close(10)
end program qmdffgen
