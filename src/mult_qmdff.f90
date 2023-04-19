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
!     ##  program mult_qmdff  --  duplicate QMDFFs of molecules      ##
!     ##                                                             ##
!     #################################################################
!
!
!     "mult_qmdff" takes one or more QMDFFs of single molecules and 
!     builds up a larger ensemble of copies, e.g. in order to simulate 
!     liquids and solvent boxes.
!
program mult_qmdff
use general

implicit none
integer::i,j,k,l  ! loop indices
integer::nqmdffs  ! number of reference qmdfffs
character(len=80)::prefix,line  ! for read in of QMDFF names
character(len=:),allocatable::prefix1
integer::l1,readstat
logical::exist
integer::rank
character(len=40)::commarg ! string for command line argument
integer::istat
character(len=85)::textout
character(len=20),allocatable::numbers(:)  ! written numeration of QMDFFs
integer,allocatable::units(:)  ! array with file input units of QMDFFs
character(len=50),allocatable::ff_names(:)  ! names of the QMDFF files
integer,allocatable::ff_atoms(:)  ! number of atoms for each reference QMDFF
real(kind=8),allocatable::line1(:,:)  ! first line of QMDFF file
real(kind=8),allocatable::atnum_ref(:,:)  ! atomic numbers of atoms in reference QMDFFs 
real(kind=8),allocatable::xyz_ref(:,:,:) ! xyz coordinates of reference QMDFFs
real(kind=8),allocatable::size_ref(:,:) ! x, y and z extents of the reference molecules
real(kind=8),allocatable::charge_ref(:,:) ! charges of atoms in reference QMDFFFs
integer,allocatable::nbonds(:),nangs(:),ntors(:),nhbnds(:),nncov(:),n12(:) ! number of terms
integer,allocatable::bond_atms(:,:,:),ang_atms(:,:,:),tors_atms(:,:,:) ! atoms of bonded interactions
integer,allocatable::hbnd_atms(:,:,:),ncov_atms(:,:,:) ! atoms of nonbonded interactions
real(kind=8),allocatable::bond_pars(:,:,:),ang_pars(:,:,:),tors_pars(:,:,:) ! parameters of bonded..
integer::natoms_max,nbonds_max,nangs_max,ntors_max,nhbnds_max,nncov_max,n12_max  ! number of FF terms 
integer,allocatable::hbnd_lines(:),ncov_lines(:)  ! number of full lines with h-bonds/noncov-interact.
integer::box_size   ! size of the box (number of molecules in x,y,z)
real(kind=8)::box_shift  !size of shift between two molecules in the box
real(kind=8),allocatable::xyz_all(:,:)   ! Coordinates of all atoms in the box/final QMDFF 
real(kind=8),allocatable::mn_par_all(:)  ! noncovalent correction parameters for all molecules
real(kind=8),allocatable::mn_par(:,:)
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
real(kind=8)::rand  ! random number between 0 and 1
integer::atoms_first  ! shift for atom numbers in bond terms etc.
integer::act_index  ! index of the current molecule to be chosen
integer::natoms_box,nbonds_box,nangs_box,ndiheds_box  ! number of terms in the full box
integer::nhbnd_box,nncov_box,n12_box  ! number of nonbonded terms in the full box
integer::num,nmols   ! number of QMDFF copies/molecules in the box
integer::atom_count,bond_count,ang_count,tors_count,hbnd_count,ncov_count  ! counter for box builtup
real(kind=8)::add_x,add_y,add_z  ! local shifts for xyz of box
character(len=3)::atname  ! atomic symbols for xyz test printout
real(kind=8)::duration,time1,time2  ! time measurement for the builtup
character(len=2)::adum
logical::ff_mod_noncov
character(len=80)::corr_name
character(len=1)::decision
integer,allocatable::linenum(:)
integer::numnci,remain

nqmdffs=1

call promo


!
!     print out help informations if -help or -h is given as command line argument
!     call the respective subroutine, else, show infos about the possibility
!
rank=0
if (rank .eq. 0) then
   call get_command_argument(1, commarg)
   if (trim(commarg) .eq. "-help" .or. trim(commarg) .eq. "-h") then
      call help("mult_qmdff")
      stop
   else
      write(*,*) "To show some basic infos about the program and a list of all"
      write(*,*) "used keywords in it, type 'mult_qmdff.x -help' or 'mult_qmdff.x -h'."
      write(*,*)
   end if
end if


allocate(numbers(10))
numbers=(/ "first  ","second ","third  ","fourth ","fifth  ","sixth  ","seventh", &
          & "eigth  ","ninth  ","tenth  " /)
write(*,*) "Welcome to the program MULT_QMDFF.X!"
write(*,*) "Here you can build up QMDFF of multi-molecule systems"
write(*,*) "by copying QMDFFs of different single molecules!"
write(*,*) 
exist=.false.
do while (.not. exist) 
   write(*,'(a)',advance="no") " How many different molecules (and their QMDFFs) are present?: "
   read(*,*,iostat=readstat) nqmdffs
   if (readstat .eq. 0) then
      if (nqmdffs .gt. 10) then
         write(*,*) "So far, only up to ten different molecules can be chosen!"
         write(*,*) "Please try again!"
      else 
         exist=.true.
      end if
   end if
end do

allocate(ff_atoms(nqmdffs))
allocate(line1(nqmdffs,2))
allocate(units(nqmdffs))
allocate(nbonds(nqmdffs),nangs(nqmdffs),ntors(nqmdffs),nhbnds(nqmdffs),nncov(nqmdffs),n12(nqmdffs))
allocate(hbnd_lines(nqmdffs),ncov_lines(nqmdffs))
allocate(size_ref(3,nqmdffs))
!
!     Read in all reference QMDFF (arbitrary number)
!
exist=.false.
allocate(ff_names(nqmdffs))
do i=1,nqmdffs
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
   units(i)=17+i
  
   open(unit=units(i),file=textout,status="old")
   ff_names(i)=textout
   exist=.false.
   deallocate(prefix1)
!     Determine the number of atoms in the QMDFF
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
   read(units(i),*) nbonds(i),nangs(i),ntors(i),nhbnds(i),nncov(i),n12(i)   
end do


!
!     Bonds for parameter arrays of an arbitrary QMDFF number 
!
nbonds_max=maxval(nbonds(:))
nangs_max=maxval(nangs(:))
ntors_max=maxval(ntors(:))
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
   do j=1,nbonds(i)
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
   do j=1,ntors(i)
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

   !   do =1,hbnd_lines(i)
   !      do k=1,5
   !         read(units(i),'(3i4)',advance="no") hbnd_atms(:,(j-1)*6+k,i)
   !      end do
   !      read(units(i),'(3i4)') hbnd_atms(:,(j-1)*6+6,i)
   !   end do
   !   do j=1,nhbnds(i)-hbnd_lines(i)*6-1
   !      read(units(i),'(3i4)',advance="no") hbnd_atms(:,hbnd_lines(i)*6+j,i)
   !   end do
   !   read(units(i),'(3i4)') hbnd_atms(:,nhbnds(i),i)
   end if
!
!    Read in the van der Waals and coulomb interactions
!
   if (nncov(i) > 0) then
   !   read(units(i),'(8(3i4,2x))')ncov_atms(1:3,1:nncov(i),i)
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


 !     do j=1,ncov_lines(i)
 !        do k=1,7
 !           read(units(i),'(3i4,a2)',advance="no") ncov_atms(:,(j-1)*8+k,i),adum            
 !        end do
 !        read(units(i),'(3i4)') ncov_atms(:,(j-1)*8+8,i)
 !     end do
 !     do j=1,nncov(i)-ncov_lines(i)*8-1
 !        read(units(i),'(3i4)',advance="no") ncov_atms(:,ncov_lines(i)*8+j,i)
 !     end do
 !     read(units(i),'(3i4)') ncov_atms(:,nncov(i),i)
   end if
end do

do i=1,nqmdffs
   close(units(i))
end do
!
!    Now build the large QMDFF for all molecules!
!

exist=.false.
do while (.not. exist)
   write(*,'(a)',advance="no") " Give the size of the box! (number of mols per dimension): "
   read(*,*,iostat=readstat) box_size
   if (readstat .eq. 0) then
      if (box_size .gt. 20) then
         write(*,*) "So far, only up to twenty molecules can be taken per box dimension!"
         write(*,*) "Please try again!"
      else
         exist=.true.
      end if
   end if
end do

!
!    Determine the distance between two molecules in the box (default: 3 Angstrom)
!
box_shift=2.d0
exist=.false.
do while (.not. exist)
   write(*,'(a)',advance="no") " Give the distance between two mols. in the grid (def: 2 A): "
   read(*,*,iostat=readstat) box_shift
   if (readstat .eq. 0) then
      if (box_shift .lt. 0.d0) then
         write(*,*) "Only positive distances can be chosen!"
      else
         exist=.true.
      end if
   else 
      box_shift=2.d0
      exist=.true.
   end if
end do
!
!    Convert bohr 
!
box_shift=box_shift/bohr
!
!    Size that each molecule occupies: largest extend of any of them in any direction + 4 bohr!
!
box_shift=maxval(size_ref(:,:))+box_shift

!
!    Number of molecules in the box...
!

nmols=box_size**3
allocate(mol_inds(nmols))

!
!    Determine the exact composition of the solvent in the box!
!    First, read in relative abundancies of the different species
!
allocate(abundance(nqmdffs))
allocate(num_species(nqmdffs))
allocate(num_act(nqmdffs))

if (nqmdffs .gt. 1) then
   do i=1,nqmdffs
      exist=.false.
      do while (.not. exist)
         write(iout,'(a,a,a)',advance="no") " Relative abundancy of the ",trim(numbers(i)), &
              & " QMDFF (arbitrary number): "
         read (*,*,iostat=readstat) abundance(i)
         if (readstat .eq. 0) then
            exist=.true.
         end if
      end do
   end do
else 
   abundance(1)=1.d0
end if
!
!    Ask if noncovalent correction terms shall be included into the multiplication
!
exist=.false.
do while (.not. exist)
   write(iout,'(a,a,a)',advance="no") " Shall correction parameters be read in from &
                     &'[name]_mod.dat'? (y/n): "
   read (*,*,iostat=readstat) decision
   if (decision .eq. "y") then
      ff_mod_noncov=.true. 
      exist=.true.
   end if
   if (decision .eq. "n") then
      ff_mod_noncov=.false.
      exist=.true.
   end if
end do


!
!    Measure the needed time for production
!
call cpu_time(time1)

abundance(:)=abundance(:)/sum(abundance(:))

!
!    Calculate the numbers of different species based on the abundancies and correct 
!    them if they do not fit the total number of molecules
!
do i=1,nqmdffs
   num_species(i)=int(nmols*abundance(i))
end do

if (sum(num_species) .lt. nmols) then
   num_species(nqmdffs)=num_species(nqmdffs)+(nmols-sum(num_species))
else if (sum(num_species) .gt. nmols) then
   num_species(nqmdffs)=num_species(nqmdffs)-(nmols-sum(num_species))
end if

write(*,*)
write(*,*) "The system will contain the following randomly placed molecules:"
do i=1,nqmdffs
   write(*,'(i4,a,a)') num_species(i)," units from ",ff_names(i)
end do
!
!    Create array for random number to be placed sampled within the abundancies 
!    of the different species 
!
allocate(distribute(nqmdffs))
do i=1,nqmdffs
   distribute(i)=sum(abundance(1:i))
end do
!
!    Loop over all three dimensions, first for the determination of array sizes!
!    Determine the molecular composition of the system!
!
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
do i=1,box_size
   do j=1,box_size
      do k=1,box_size
         num=num+1
!
!    Sample the added molecules randomistically, but asure that their predefined 
!    numbers are not exceeded
!
22       call random_number(rand)
         do l=1,nqmdffs
            if (rand .lt. distribute(l)) then
               if (num_act(l) .eq. num_species(l)) then
                  goto 22
               end if
               act_index=l
               num_act(act_index)=num_act(act_index)+1
               exit
            end if
         end do
         !write(*,*)rand, act_index
         !act_index=1
!
!    Determine the number of atoms of the first molecule for subsequent shifts etc
!
         if (num .eq. 1) then
            atoms_first=ff_atoms(act_index)
         end if
         mol_inds(num)=act_index 
         natoms_box=natoms_box+ff_atoms(act_index)
         nbonds_box=nbonds_box+nbonds(act_index)
         nangs_box=nangs_box+nangs(act_index)
         ndiheds_box=ndiheds_box+ntors(act_index) 
         nhbnd_box=nhbnd_box+nhbnds(act_index)
         nncov_box=nncov_box+nncov(act_index)
         n12_box=n12_box+n12(act_index)
      end do
   end do
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
!     If the noncovalent corrections shall be applied
!
if (ff_mod_noncov) then
   allocate(mn_par_all(1+7*natoms_box))
   allocate(mn_par(nqmdffs,1+7*maxval(ff_atoms)))
   do i=1,nqmdffs
      corr_name=ff_names(i)(1:len(trim(ff_names(i)))-6) // "_mod.dat"
      open(unit=67,file=corr_name,status="old")
      read(67,*)
      read(67,*) mn_par(i,1)   ! the global charge scaling
      do j=1,ff_atoms(i)
         read(67,*) mn_par(i,2+(j-1)*7:1+j*7)
      end do
      close(67)
   end do
!
!     Currently: Take the global charge scaling of the first QMDFF
!      also as global charge scaling for the whole box 
!
   mn_par_all(1)=mn_par(1,1)
end if
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
do i=1,box_size
   add_x=add_x+box_shift
   add_y=-box_shift
   add_z=-box_shift
   do j=1,box_size
      add_y=add_y+box_shift
      add_z=-box_shift
      do k=1,box_size
         num=num+1
         act_index=mol_inds(num)
         add_z=add_z+box_shift
!
!     First fill the coordinate section
!         
         do l=1,ff_atoms(act_index)
            atom_count=atom_count+1
            atnum_all(atom_count)=atnum_ref(l,act_index)
            charge_all(atom_count)=charge_ref(l,act_index)
            nmol_all(atom_count)=num
            xyz_all(1,atom_count)=xyz_ref(1,l,act_index)+add_x
            xyz_all(2,atom_count)=xyz_ref(2,l,act_index)+add_y
            xyz_all(3,atom_count)=xyz_ref(3,l,act_index)+add_z
!
!     If ordered, fill the global array for noncovalent interaction corrections
!
            if (ff_mod_noncov) then
               mn_par_all(2+(atom_count-1)*7:1+atom_count*7)= &
                               & mn_par(act_index,2+(l-1)*7:1+l*7)
            end if

         end do
!
!     Second, fill the bonds parameters section
!
         do l=1,nbonds(act_index)
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
         do l=1,ntors(act_index)
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
   end do
end do

!
!     Write out the xyz coordinates of the generated box for visualization
!
open(unit=278,file="box.xyz",status="replace")
write(278,*) natoms_box
write(278,*)
do i=1,natoms_box
   call atomname(atnum_all(i),atname)
   call upcase(atname)
   write(278,*) atname,xyz_all(:,i)*bohr
end do
close(278)
write(*,*)
write(*,*) "The molecular configuration of the system was written to 'box.xyz'..."

!
!     Write out the final force field!
!
open(unit=279,file="box.qmdff",status="replace")
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

write(*,*) 
write(*,'(a,f12.4,a)') " The size of the box (for periodic calculations) is ", &
   &maxval(xyz_all)*bohr+2.d0, " Angstroms."
write(*,*)
write(*,*) "The new QMDFF was generated successfully! It is named 'box.qmdff'"

!
!    If noncovalent corrections were available, write out the resulting file for them
!
if (ff_mod_noncov) then
   open(unit=238,file="box_mod.dat",status="replace")
   write(238,*) "** This file contains correction parameters for QMDFF nonbonded interactions **"
   write(238,'(f14.8)') mn_par_all(1)
   do i=1,natoms_box
      write(238,'(7f14.8)') mn_par_all(2+(i-1)*7:1+i*7)
   end do
   close(238)
end if
!
!    calculate the needed time for optimization
!
call cpu_time(time2)

duration=time2-time1
write(*,*)
write(*,'(A, F12.3, A)') " The calculation needed a time of",duration," seconds."
write(*,*)



end program mult_qmdff

