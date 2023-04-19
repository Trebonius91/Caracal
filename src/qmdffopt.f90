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
!     ##  program qmdffopt  --  optimizing QMDFF parameters          ##
!     ##                                                             ##
!     #################################################################
!
!     "qmdffopt" optimizes parameters of already generated QMDFFs in order 
!     to reproduce reference QM data even better (e.g., for simulation
!     of liquids with QMDFF).
!     The optimization consists of two parts:     
!     1.: Read in the QMDFF and generate a number of dimer structures,
!       where relative orientations are sampled by chance. These structures
!       are written to file 
!     2.: After having calculated QM energies of the structures externally,
!       the second part of the program does MSLS-LM in order to optimize 
!       the chosen FF parameters to reproduce the reference as good as 
!       possible.
!
program qmdffopt 
use general
use evb_mod
use qmdff
implicit none
integer::i,j,k,l
integer::rank
character(len=70)::fffile1,fffile2,fffile3,line
character(len=40)::commarg
integer::nstruct ! number of reference structures to be written
real(kind=8)::r_min,r_max,r_act  ! the R coordinate (spherical coordinates)
real(kind=8)::theta_min,theta_max,theta_act  ! the theta angle (polar angle)
real(kind=8)::phi_min,phi_max,phi_act  ! the phi angle (azimuth angle)
real(kind=8)::alpha_min,alpha_max,alpha_act  ! the Eulerian Nutation angle alpha
real(kind=8)::beta_min,beta_max,beta_act   ! the Eulerian Precession angle beta
real(kind=8)::gamma_min,gamma_max,gamma_act ! the Eulerian rotation angle gamma
integer::readstat
logical::exist,printhelp
real(kind=8)::com_1(3),com_2(3)  ! center of masses for the molecules
real(kind=8),allocatable::xyz_m2(:,:)  ! coordinates of second, moved molecule
real(kind=8)::euler_mat(3,3)  ! the euler rotation matrix
real(kind=8)::mass_sum 
real(kind=8)::dist,dist_vec(3)  ! check for atomic distances after random placement
integer::prog_mode !  do either preparation or optimization
!  For the optimization part
character(len=:),allocatable::prefix1
integer,allocatable::units(:)  ! array with file input units of QMDFFs
character(len=20),allocatable::numbers(:)  ! written numeration of QMDFFs
character(len=50),allocatable::ff_names(:)  ! names of the QMDFF files
integer,allocatable::ff_atoms(:)  ! number of atoms for each reference QMDFF
real(kind=8),allocatable::line1(:,:)  ! first line of QMDFF file
real(kind=8),allocatable::atnum_ref(:,:)  ! atomic numbers of atoms in reference QMDFFs 
real(kind=8),allocatable::xyz_ref(:,:,:) ! xyz coordinates of reference QMDFFs
real(kind=8),allocatable::size_ref(:,:) ! x, y and z extents of the reference molecules
real(kind=8),allocatable::charge_ref(:,:) ! charges of atoms in reference QMDFFFs
integer,allocatable::nbonds2(:),nangs(:),ntors2(:),nhbnds(:),nncov(:),n12(:) ! number of terms
integer,allocatable::bond_atms(:,:,:),ang_atms(:,:,:),tors_atms(:,:,:) ! atoms of bonded interactions
integer,allocatable::hbnd_atms(:,:,:),ncov_atms(:,:,:) ! atoms of nonbonded interactions
real(kind=8),allocatable::bond_pars(:,:,:),ang_pars(:,:,:),tors_pars(:,:,:) ! parameters of bonded..
integer::natoms_max,nbonds_max,nangs_max,ntors_max,nhbnds_max,nncov_max,n12_max  ! number of FF terms 
integer,allocatable::hbnd_lines(:),ncov_lines(:)  ! number of full lines with h-bonds/noncov-interact.
character(len=2)::adum
!   for readin of keyfile
character(len=20) keyword
character(len=60)::names,a80
character(len=120) record
character(len=120) string
integer::next

integer,allocatable::linenum(:)
integer::numnci,remain,l1,nqmdffs
character(len=60)::fileenergy,filegeo,textout
real(kind=8)::box_shift  !size of shift between two molecules in the box
real(kind=8),allocatable::xyz_all(:,:)   ! Coordinates of all atoms in the box/final QMDFF 
integer,allocatable::atnum_all(:)  ! atomic numbers of box atoms
real(kind=8),allocatable::charge_all(:)  ! charges of box atoms
integer,allocatable::bond_atms_all(:,:),ang_atms_all(:,:),tors_atms_all(:,:) ! atoms of bonded...
integer,allocatable::hbnd_atms_all(:,:),ncov_atms_all(:,:) ! atoms of nonbonded interactions 
real(kind=8),allocatable::bond_pars_all(:,:),ang_pars_all(:,:),tors_pars_all(:,:) ! parameters...
integer,allocatable::nmol_all(:)  ! molecule number for each atom in the box
integer,allocatable::mol_inds(:)  ! Which molecule occupies which site in the box
real(kind=8),allocatable::abundance(:)  ! Relative abundancy of different molecules 
integer,allocatable::num_species(:),num_act(:)  !   Absolute number of the different molecules
real(kind=8),allocatable::distribute(:)  ! Auxiliary array for random sampling of components
integer::natoms_box,nbonds_box,nangs_box,ndiheds_box  ! number of terms in the full box
integer::nhbnd_box,nncov_box,n12_box  ! number of nonbonded terms in the full box
integer::atom_count,bond_count,ang_count,tors_count,hbnd_count,ncov_count  ! counter for box builtup
integer::act_index,atoms_first,num  ! index of the current molecule to be chosen
real(kind=8)::add_x
integer::points ! number of structures in the trainset
!
!     Start MPI parallel computation:
!     Since the whole program is executed by all processes, all "serial"
!     parts will be executed only by processor #0
!
!call mpi_init(ierr)
!call mpi_comm_size(mpi_comm_world,psize,ierr) ! total number of procs
!call mpi_comm_rank(mpi_comm_world,rank,ierr) ! index of current proc
rank=0
!
!     set up the structure and mechanics calculation
!
call initial(rank)
!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
prog_mode=0
printhelp=.true.
if (rank .eq. 0) then
   do i=1,10
      call get_command_argument(i, commarg,readstat)
      if (readstat .ne. 0) then
         if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
            call help("qmdffopt")
            printhelp=.false.
            stop
         end if
         if (trim(commarg) .eq. "-prepare") then
            write(*,*) "The prepare mode was started, a file 'trainset.xyz' with the"
            write(*,*) " structures of a trainset for following optimization will& 
                 & be written."
            write(*,*)
            prog_mode=1
         else if (trim(commarg) .eq. "-optimize") then
            write(*,*) "The optimization mode was started, the QMDFF parameters will"
            write(*,*) " be optimized with MSLS-LM."
            write(*,*)
            prog_mode=2
         end if
      end if
   end do
end if
if (printhelp) then
   write(*,*) "To show some basic infos about the program and a list of all"
   write(*,*) "used keywords in it, type 'qmdffopt.x -help' or 'qmdffopt.x -h'."
   write(*,*)
end if

if (prog_mode .eq. 0) then
   write(*,*)
   write(*,*) "Please choose either the preparation mode '-prepare' or the "
   write(*,*) " optimization mode '-optimize'!"
   stop
end if
!
!     Open the global output-file for EVB-term optimizations
!
if (rank .eq. 0) then
   open(unit=15,file="qmdffopt.log",status="unknown")
   call promo_log(15)
end if

!
!     Read keywords from the key-file: only for the optimizing mode!
!
if (prog_mode .eq. 2) then
   call getkey(rank)
end if

!
!     Read in the Force field Parameters
!
E_zero1=0d0
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'EQMDFF ') then
      read(record,*) names,E_zero1
   end if
end do
nqmdff=0
fffile1="" 
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'FFNAME ') then
      read(record,*) names,fffile1
      nqmdff=1
   end if
end do
if (nqmdff .eq. 0) then
   exist=.false.
   do while (.not. exist)
      write(iout,'(a,a,a)',advance="no") " Name of the QMDFF &
          & whose parameters shall be optimized:  "
      read (*,'(A80)')  fffile1
      inquire(file=fffile1,exist=exist)
   end do
end if

!-------------------------------------------------------------------------------
!       PART 1: Preparation of QM input 
!-------------------------------------------------------------------------------
if (prog_mode .eq. 1) then
!
!    Define the filenames of trainset
!
   filegeo="trainset.xyz"
   fileenergy="trainset_ens.dat"
   exist=.false.
   do while (.not. exist)
      write(iout,'(a)',advance="no") " How many structures shall be generated for &
                 &reference?: "
      read(*,*,iostat=readstat) nstruct
      if (readstat .eq. 0) then
         if (nstruct .gt. 0 .and. nstruct .lt. 1E7) then
            exist=.true.
         end if 
      end if
   end do

   exist=.false.
   do while (.not. exist)
      write(iout,'(a)',advance="no") " Minimum distance of the molecules?: "
      read(*,*,iostat=readstat) r_min
      if (readstat .eq. 0) then
         if (r_min .gt. -0.00001d0 .and. r_min .lt. 1E4) then
            exist=.true.
         end if
      end if
   end do
   exist=.false.
   do while (.not. exist)
      write(iout,'(a)',advance="no") " Maximum distance of the molecules?: "
      read(*,*,iostat=readstat) r_max
      if (readstat .eq. 0) then
         if (r_max .gt. 0.1d0 .and. r_max .lt. 1E4) then
            exist=.true.
         end if
      end if
   end do



!
!     Initialize the QMDFF
!
   call prepare (fffile1,fffile2,fffile3,1)
   natoms=n_one
!
!     Initialize parameters for structure sampling (might be read in later)
!
!     The Spherical coordinates for the COM movement of the second molecule
!
   r_min=r_min/bohr
   r_max=r_max/bohr
   theta_min=0.d0
   theta_max=pi
   phi_min=-pi
   phi_max=pi
!
!     The Eulerian angles for the rotation of the second molecule arount itself
!
   alpha_min=0.d0
   alpha_max=pi
   beta_min=0.d0
   beta_max=pi
   gamma_min=0.d0
   gamma_max=2*pi

!
!     Determine element symbols for trajectory printout and the atomic masses
!
   do i=1,natoms
      call atomname(at(i),name(i))
      call upcase(name(i))
      call atommass(i)
   end do

!
!     Calculate the center of mass of the central molecule 
!
   com_1=0.d0
   mass_sum=sum(mass(1:natoms))
   do i=1,natoms
      do j=1,3
         com_1(j)=com_1(j)+mass(i)*xyz(j,i)/mass_sum
      end do
   end do
!
!     Move molecule to the origin
!
   do i=1,natoms
      xyz(:,i)=xyz(:,i)-com_1(:)
   end do

!
!     Determine element symbols for trajectory printout 
!
   do i=1,natoms
      call atomname(at(i),name(i))
      call upcase(name(i))
   end do
!
!     Now enter the loop for generating structures for subsequent reference 
!     calculations    
!
!     Initialize the random number generator
!
   call random_init_local(rank)
   allocate(xyz_m2(3,natoms))

!
!     Open trajectory file for rotated molecules
!
   open(unit=14,file="trainset.xyz",status="replace")
!
!     Open orca xyz input file for rotated molecules 
!
   open(unit=15,file="trainset_orca.xyz",status="replace")

   do i=1,nstruct
      17 continue   ! if the molecules were too near, go there
!
!     Obtain random values for all six internal variables 
!
!     1. The radius: distance of both COMs     
      call random_number(r_act)
      r_act=r_min+r_act*(r_max-r_min)
!
!     2. The polar angle theta: COM rotation in z direction
!
      call random_number(theta_act)
      theta_act=theta_min+theta_act*(theta_max-theta_min)
!
!     3. The azimuth angle phi: COM rotation in x-y direction
!
      call random_number(phi_act)
      phi_act=phi_min+phi_act*(phi_max-phi_min)
!
!     4. The Eulerian Nutation angle alpha: rotate molecule 2 in z-acis
!
      call random_number(alpha_act)
      alpha_act=alpha_min+alpha_act*(alpha_max-alpha_min)
!
!     5. The Eulerian Precession angle beta: rotate molecule 2 in x-y plane
!
      call random_number(beta_act)
      beta_act=beta_min+beta_act*(beta_max-beta_min)
!
!     6. The Eulerian Rotation angle gamma: rotate molecule 2 in x' plane
!
      call random_number(gamma_act)
      gamma_act=gamma_min+gamma_act*(gamma_max-gamma_min)

!
!     Calculate the Euler rotation matrix
!
      euler_mat(1,1)=cos(gamma_act)*cos(beta_act)-sin(gamma_act)*&
                      & cos(alpha_act)*sin(beta_act)
      euler_mat(2,1)=-sin(gamma_act)*cos(beta_act)-cos(gamma_act)*&
                      & cos(alpha_act)*sin(beta_act)
      euler_mat(3,1)=sin(alpha_act)*sin(beta_act)
      euler_mat(1,2)=cos(gamma_act)*sin(beta_act)+sin(gamma_act)*&
                      & cos(alpha_act)*cos(beta_act)
      euler_mat(2,2)=-sin(gamma_act)*sin(beta_act)+cos(gamma_act)*&
                      & cos(alpha_act)*cos(beta_act)
      euler_mat(3,2)=-sin(alpha_act)*cos(beta_act)
      euler_mat(1,3)=sin(gamma_act)*sin(alpha_act)
      euler_mat(2,3)=cos(gamma_act)*sin(alpha_act)
      euler_mat(3,3)=cos(alpha_act)
!
!     Apply the Euler rotation to the coordinates of the centralized 
!        second molecule 
!
      xyz_m2=xyz
      do j=1,natoms
         xyz_m2(:,j)=matmul(xyz_m2(:,j),euler_mat)
      end do
!
!     Move COM of second molecule to the desired place (spherical coordinates)
!
      com_2(1)=r_act*sin(theta_act)*cos(phi_act)
      com_2(2)=r_act*sin(theta_act)*sin(phi_act)
      com_2(3)=r_act*cos(theta_act)


!
!     Move coordinates of second molecule to new COM
!
      do j=1,natoms
         xyz_m2(:,j)=xyz_m2(:,j)+com_2
      end do
!
!     Check if any collisions or too narrow atoms occur in the dimer
!     If yes, repeat the whole process for this cycle
!
      do j=1,natoms 
         do k=1,natoms
            dist_vec=xyz(:,j)-xyz_m2(:,k)
         
            dist=sqrt(dot_product(dist_vec,dist_vec))
            if (dist .lt. (rad(at(j))+rad(at(k)))*1.4) then
               goto 17
            end if           

         end do
      end do
!
!     Write actual finished structure to file 
!
      r_act=r_act*bohr
      theta_act=theta_act*180.d0/pi
      phi_act=phi_act*180.0/pi
      write(14,*) 2*natoms 
      write(14,'(a,i5,a,f11.6,a,f11.6,a,f11.6,a)') " Mol. ",i," r=",r_act," A, theta=", & 
                       & theta_act,"°, phi=",phi_act,"°" 
      write(15,*) 2*natoms 
      write(15,*)
      do j=1,natoms
         write(14,*) name(j),xyz(:,j)*bohr
         write(15,*) name(j),xyz(:,j)*bohr
      end do
      do j=1,natoms
         write(14,*) name(j),xyz_m2(:,j)*bohr
         write(15,*) name(j),xyz_m2(:,j)*bohr
      end do
      if (i .ne. nstruct) then
         write(15,'(a)') ">"
      else 
         write(15,*) " "
      end if
   end do
   close(14)
   close(15)
   write(*,*) 
   write(*,*) "Trainset generation was successful!"
   write(*,*) "Files 'trainset.xyz' and 'trainset_orca.xyz' written!"
   write(*,*) "Please calculate QM reference energies for these structures and "
   write(*,*) " submit them in the file trainset_ens.dat for the second part."
   write(*,*) 
end if
!-------------------------------------------------------------------------------
!       PART 2: Optimization of QMDFF parameters
!-------------------------------------------------------------------------------
if (prog_mode .eq. 2) then
!
!     First, generate an intermediary QMDFF for the molecule dimer
!     INFO: this code is taken from the mult_qmdff program in order to enable 
!     later generalizations for more than one QMDFF!
!
   allocate(numbers(10))
   numbers=(/ "first  ","second ","third  ","fourth ","fifth  ","sixth  ","seventh", &
          & "eigth  ","ninth  ","tenth  " /)

   nqmdffs=1
   allocate(ff_atoms(nqmdffs))
   allocate(line1(nqmdffs,2))
   allocate(units(nqmdffs))
   allocate(nbonds2(nqmdffs),nangs(nqmdffs),ntors2(nqmdffs),nhbnds(nqmdffs),&
               & nncov(nqmdffs),n12(nqmdffs))
   allocate(hbnd_lines(nqmdffs),ncov_lines(nqmdffs))
   allocate(size_ref(3,nqmdffs))

!
!     Read in the settings for the Multi-Start-Local-Search Levenberg-Marquardt 
!     algorithm!
!
   optsteps=100
   maxstart=100
   lower_bond=0.5d0
   upper_bond=1.5d0
   lm_par=0.01d0
   lm_threshold=1E-8
   lm_par_change=2d0
   diff_step=0.0001d0
   do i = 1, nkey
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
      end if
   end do

   allocate(ff_names(nqmdffs))
   do i=1,nqmdffs
      goto 127
      do while (.not. exist)
         write(iout,'(a,a,a)',advance="no") " Prefix of the ",trim(numbers(i)), &
              & " QMDFF (name.qmdff must be in this folder!): "
         read (*,'(A80)')  line
         allocate(character(len=LEN(TRIM(line))) :: prefix1)
         prefix1=trim(line)
         l1=LEN(TRIM(prefix1))
         textout = prefix1 // ".qmdff"
         inquire(file=textout,exist=exist)
         if (.not. exist) then
            deallocate(prefix1)
         end if
      end do
      127 continue
      units(i)=17+i
      textout=fffile1
      open(unit=units(i),file=textout,status="old")
      ff_names(i)=textout
      exist=.false.
 !     deallocate(prefix1)
!        Determine the number of atoms in the QMDFF
      read(units(i),*) ff_atoms(i),line1(i,1),line1(i,2)
      read(units(i),*)
   end do
!
!     Take the largest atom number as reference for array initializations
!
   natoms_max=maxval(ff_atoms(:))
!
!     Allocate arrays with atom numbers, coordinates and charges 
!
   allocate(atnum_ref(natoms_max,nqmdffs))
   allocate(xyz_ref(3,natoms_max,nqmdffs))
   allocate(charge_ref(natoms_max,nqmdffs))

   do i=1,nqmdffs
!
!     Read in the coordinate section of the qmdffs
!
      do j=1,ff_atoms(i)
         read(units(i),*) atnum_ref(j,i),xyz_ref(:,j,i),charge_ref(j,i)
      end do
!
!     Shift molecule towards coordinate origin for uniformity
!
      xyz_ref(1,:,i)=xyz_ref(1,:,i)-minval(xyz_ref(1,:,i))
      xyz_ref(2,:,i)=xyz_ref(2,:,i)-minval(xyz_ref(2,:,i))
      xyz_ref(3,:,i)=xyz_ref(3,:,i)-minval(xyz_ref(3,:,i))


!
!     Determine molecule sizes in x,y and z directions
!
      size_ref(1,i)=maxval(xyz_ref(1,:,i))
      size_ref(2,i)=maxval(xyz_ref(2,:,i))
      size_ref(3,i)=maxval(xyz_ref(3,:,i))
!
!     Read in the number of ffield terms
!
      read(units(i),*) nbonds2(i),nangs(i),ntors2(i),nhbnds(i),nncov(i),n12(i)
   end do

!
!     Bonds for parameter arrays of an arbitrary QMDFF number 
!
   nbonds_max=maxval(nbonds2(:))
   nangs_max=maxval(nangs(:))
   ntors_max=maxval(ntors2(:))
   nhbnds_max=maxval(nhbnds(:))
   nncov_max=maxval(nncov(:))
   n12_max=maxval(n12(:))

!
!      Allocate the bonded and nonbonded parameter arrays
!
   allocate(bond_atms(2,nbonds_max,nqmdffs))
   allocate(ang_atms(3,nangs_max,nqmdffs))
   allocate(tors_atms(6,ntors_max,nqmdffs))
   allocate(hbnd_atms(3,nhbnds_max,nqmdffs))
   allocate(ncov_atms(3,nncov_max,nqmdffs))
   allocate(bond_pars(3,nbonds_max,nqmdffs))
   allocate(ang_pars(2,nangs_max,nqmdffs))
   allocate(tors_pars(5,ntors_max,nqmdffs))
! 
!      read parameters of covalent bonded section
!
   do i=1,nqmdffs
      hbnd_lines(i)=nhbnds(i)/6
      ncov_lines(i)=nncov(i)/8
!
!     Read in the bond section
!
      do j=1,nbonds2(i)
         read(units(i),*) bond_atms(:,j,i),bond_pars(:,j,i)
      end do
!
!     Read in the angle section
! 
      do j=1,nangs(i)
         read(units(i),*) ang_atms(:,j,i),ang_pars(:,j,i)
      end do
!
!    Read in the torsion section
!
      do j=1,ntors2(i)
         read(units(i),*) tors_atms(:,j,i),tors_pars(:,j,i)
      end do
!
!    Read in the hydrogen bond section
!
      if (nhbnds(i) > 0) then
         numnci=0
         allocate(linenum(18))
         do k=1,int(nhbnds(i)/6)
            read(units(i),*) linenum
            do j=1,6
               numnci=numnci+1
               hbnd_atms(1:3,numnci,i)=linenum((j-1)*3+1:(j-1)*3+3)
            end do
         end do
         deallocate(linenum)
         remain=(nhbnds(i)-numnci)*3
         allocate(linenum(remain))
         read(units(i),*) linenum
         do j=1,remain
            numnci=numnci+1
            hbnd_atms(1:3,numnci,i)=linenum((j-1)*3+1:(j-1)*3+3)
         end do
         deallocate(linenum)
      end if
!
!    Read in the van der Waals and coulomb interactions
!
      if (nncov(i) > 0) then
         numnci=0
         allocate(linenum(24))
         do k=1,int(nncov(i)/8)
            read(units(i),*) linenum
            do j=1,8
               numnci=numnci+1
               ncov_atms(1:3,numnci,i)=linenum((j-1)*3+1:(j-1)*3+3)
            end do
         end do
         deallocate(linenum)
         remain=(nncov(i)-numnci)*3
         allocate(linenum(remain))
         read(units(i),*) linenum
         do j=1,remain
            numnci=numnci+1
            ncov_atms(1:3,numnci,i)=linenum((j-1)*3+1:(j-1)*3+3)
         end do
         deallocate(linenum)

      end if
   end do

   do i=1,nqmdffs
      close(units(i))
   end do

!
!    Write the duplicate QMDFF to file 'dimer.qmdff'
!
   nmols=2
!
!    Convert bohr 
!
   box_shift=2.d0
   box_shift=box_shift/bohr
!
!    Size that each molecule occupies: largest extend of any of them in any direction + 4 bohr!
!
   box_shift=maxval(size_ref(:,:))+box_shift

   allocate(mol_inds(nmols))

   allocate(abundance(nqmdffs))
   allocate(num_species(nqmdffs))
   allocate(num_act(nqmdffs))

   filegeo="trainset.xyz"
   fileenergy="trainset_ens.dat"
   num=0   ! the number/index of the current molecule 
   natoms_box=0
   nbonds_box=0
   nangs_box=0
   ndiheds_box=0
   nhbnd_box=0
   nncov_box=0
   n12_box=0

   call random_init_local(0)
   num_act(:)=0
   do i=1,2

      num=num+1
      act_index=1  ! only one QMDFF is used
!
!    Determine the number of atoms of the first molecule for subsequent shifts etc
!
      if (num .eq. 1) then
         atoms_first=ff_atoms(act_index)
      end if
      mol_inds(num)=act_index
      natoms_box=natoms_box+ff_atoms(act_index)
      nbonds_box=nbonds_box+nbonds2(act_index)
      nangs_box=nangs_box+nangs(act_index)
      ndiheds_box=ndiheds_box+ntors2(act_index)
      nhbnd_box=nhbnd_box+nhbnds(act_index)
      nncov_box=nncov_box+nncov(act_index)
      n12_box=n12_box+n12(act_index)
   end do
!
!     Allocate global arrays for final QMDFF output 
!

   allocate(xyz_all(3,natoms_box))
   allocate(atnum_all(natoms_box))
   allocate(charge_all(natoms_box))
   allocate(nmol_all(natoms_box))
   allocate(bond_atms_all(2,nbonds_box))
   allocate(ang_atms_all(3,nangs_box))
   allocate(tors_atms_all(6,ndiheds_box))
   allocate(hbnd_atms_all(3,nhbnd_box))
   allocate(ncov_atms_all(3,nncov_box))
   allocate(bond_pars_all(3,nbonds_box))
   allocate(ang_pars_all(2,nangs_box))
   allocate(tors_pars_all(5,ndiheds_box))
!
!     Fill all global arrays of the box!
!
   num=0
   atom_count=0
   bond_count=0
   ang_count=0
   tors_count=0
   hbnd_count=0
   ncov_count=0
   add_x=-box_shift
   do i=1,2
      num=num+1
      add_x=add_x+box_shift
!
!     First fill the coordinate section
!         
      do l=1,ff_atoms(act_index)
         atom_count=atom_count+1
         atnum_all(atom_count)=atnum_ref(l,act_index)
         charge_all(atom_count)=charge_ref(l,act_index)
         nmol_all(atom_count)=num
         xyz_all(1,atom_count)=xyz_ref(1,l,act_index)+add_x
         xyz_all(2,atom_count)=xyz_ref(2,l,act_index)
         xyz_all(3,atom_count)=xyz_ref(3,l,act_index)
      end do
!
!     Second, fill the bonds parameters section
!
      do l=1,nbonds2(act_index)
         bond_count=bond_count+1
         bond_atms_all(:,bond_count)=bond_atms(:,l,act_index)+atom_count-ff_atoms(act_index)
         bond_pars_all(:,bond_count)=bond_pars(:,l,act_index)
      end do
!
!     Third, fill the angle parameters section
!
      do l=1,nangs(act_index)
         ang_count=ang_count+1
         ang_atms_all(:,ang_count)=ang_atms(:,l,act_index)+atom_count-ff_atoms(act_index)
         ang_pars_all(:,ang_count)=ang_pars(:,l,act_index)
      end do
!
!     Fourth, fill the torsion parameters section
!
      do l=1,ntors2(act_index)
         tors_count=tors_count+1
         tors_atms_all(1:4,tors_count)=tors_atms(1:4,l,act_index)+atom_count-ff_atoms(act_index)
         tors_atms_all(5:6,tors_count) = tors_atms(5:6,l,act_index)
         tors_pars_all(:,tors_count)=tors_pars(:,l,act_index)
      end do
!
!     Fifth, fill the hydrogen bond parameters section
!
      do l=1,nhbnds(act_index)
         hbnd_count=hbnd_count+1
         hbnd_atms_all(:,hbnd_count)=hbnd_atms(:,l,act_index)+atom_count-ff_atoms(act_index)
      end do
!
!     Sixth, fill the noncovalent bond parameters section
!         
      do l=1,nncov(act_index)
         ncov_count=ncov_count+1
         ncov_atms_all(1:2,ncov_count)=ncov_atms(1:2,l,act_index)+atom_count-ff_atoms(act_index)
         ncov_atms_all(3,ncov_count)=ncov_atms(3,l,act_index)
      end do
   end do

!
!     Write out the final force field!
!
   open(unit=279,file="dimer.qmdff",status="replace")
   n=natoms_box
   write(279,'(i4,2F10.4)') n,real(nmols),0d0
   write(279,'(a)') trim("CM5*1.15")
   do i=1,n
      write(279,'(i5,4F18.12,i7)')atnum_all(i),xyz_all(1:3,i),charge_all(i),nmol_all(i)
   end do
   write(279,'(6i10)') nbonds_box,nangs_box,ndiheds_box,nhbnd_box,nncov_box,n12_box
   do i=1,nbonds_box
      write(279,'(2i7,3F18.12)')bond_atms_all(1,i),bond_atms_all(2,i),bond_pars_all(1:3,i)
   end do
   do i=1,nangs_box
      write(279,'(3i7,5F18.12)')ang_atms_all(1,i),ang_atms_all(2,i),ang_atms_all(3,i), &
             & ang_pars_all(1:2,i)
   end do
   do i=1,ndiheds_box
      write(279,'(4i7,2i4,2F16.12,10(2F4.1,E16.8))') &
        &   tors_atms_all(1:6,i),tors_pars_all(:,i)
   end do
   if (nhbnd_box.gt.0) write(279,'(6(3i6,2x))')hbnd_atms_all(1:3,1:nhbnd_box)
   write(279,'(8(3i6,2x))')ncov_atms_all(1:3,1:nncov_box)

   close(279)
!
!    Now call the Levenberg-Marquardt routine!
!    First, determine the number of items in the trainset
!
   open(unit=13,file=fileenergy,iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",trim(fileenergy)," could not be found!"
      call fatal
   end if        
   points=0
   do
      read(13,*,iostat=readstat) 
      if (readstat .ne. 0) exit
      points=points+1
   end do
   close(13)
   natoms=n
   call opt_qmdff_ser(filegeo,fileenergy,points,fffile1)

end if


end program qmdffopt
