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
!     subroutine init_int: This routine is called for every 
!     dQ-EVB or DG-EVB calculation
!     to initialize and specify the used set of internal coordinates
!     These can be read in from file coord_def.inp or my be defined 
!     on the fly by analyzing the fluctuations on the reaction path
!     (then there are num_coord coordinates)
!
!     part of EVB
!
subroutine init_int(filegeo,points,rank,mode)
use evb_mod
use qmdff
use general

implicit none
integer::i,j,k,l,o,p,points
integer::nat,nat3
integer::atind(natoms)  ! atomic indices
real(kind=8)::dist,ang,oop,dihed
real(kind=8),dimension(3,natoms)::xyz2,coord 
real(kind=8),dimension(:),allocatable::rmsd_path_int 
real(kind=8)::minimum,temp  
real(kind=8),allocatable::temp3(:)  
integer,allocatable::temp2(:)  
character(len=20) keyword
character(len=120) record
character(len=120) string 
integer::location,start,ihalf,next
integer,dimension(:,:),allocatable::coord_tmp,coord_reshape,bond_add
integer,dimension(:),allocatable::coord_types
real(kind=8),dimension(:,:,:),allocatable::xyz_path  ! xyz structures of the reactionpath
real(kind=8),dimension(:,:),allocatable::path_int
real(kind=8),dimension(:,:),allocatable::path_xyz,coord_test
real(kind=8),dimension(:),allocatable::int_path,int_coord  ! for internal coordinates
integer::mode  !  which procedure shall be done for defining coordinates: 
               !  1: read them from file coord_def.inp
               !  2: take the N coordinates with largest fluctuations on the path
               !  3: for rpmd.x with RP-EVB: analyze bond orders and fill angles 
integer::readstat  ! catch IO exceptions
integer::coordnums(3)  ! number of bonds, angles, dihedrals in one array
real(kind=8),dimension(:),allocatable::coord_vals,coord_vals_tmp  ! values of all internal coordinates
real(kind=8)::bondlentol  ! tolerance factor for covalent bondlengths in order to establish a bond
real(kind=8),dimension(:),allocatable::int_test(:),grad_int(:)  ! test for sanity of wilson matrix
real(kind=8)::grad_xyz(3*natoms)  ! test for sanity of wilson matrix
real(kind=8),dimension(:),allocatable::hess_int(:,:)  ! test for sanity of wilson matrix derivative
real(kind=8)::hess_xyz(3*natoms,3*natoms)  ! test for sanity of wilson matrix derivative

integer::nbonds_old,nbonds_add  ! stored result of bond number for the previous tolerance
integer::remkind,rempos(100)  ! type and number of coordinate to be removed if too many are there
integer::nat6_check  ! first number of checked internal coordinates (angles)
integer::double_bond  ! if a bond is already stored in bondorder matrix 1
integer::ncombi  ! number of combinations for reduced coordinate sets
real(kind=8)::rcombi ! number of combinations for reduced coordinate sets
integer::ndead  ! number of "dead" internal coordinates that shall be included
integer::nadd  ! number of added coordinates for the later stage of fillup
integer::nat6_base  ! number of coordinates in the basis set to be completed
integer::ignore ! switch if the current combination of coordinates shall be ignored
real(kind=8)::fact ! function fact to call (factorial)
integer::i1,i2,i3,i4,i5  ! for manual fill of array
integer,dimension(:,:),allocatable::outtakes,fillup  ! array with combinations to exclude/add
integer,dimension(:,:),allocatable::coord_base  ! found coordinate set that is still incomplete
logical::fillup_done ! if the coordinate set was filled up or not
integer::atm1,atm2,atm3,atm4  ! atom numbers for angle determination
character(len=*)::filegeo
character(len=40)::coord_line
logical::has_next
!  only for MPI calculations!
integer::rank

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
!
!     Convert radii to bohr
!
rad=rad/bohr
!rad(1)=rad(1)*1.5  ! shift H-atom value a bit to allow H+H2
!
!     Set default values for numbers of coordinate types
!
nbonds=0
nangles=0
noops=0
ndiheds=0
!
!     set tolerance factor for bondlengths (if a bond shall still be defined..)
!
bondlentol=1.2d0

! 
!     decide if Wilson matrix for conversion from/in internal coordinates shall
!     be calculated numerically or analytically
!
num_wilson=.true.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:15) .eq. 'NUM_WILSON ') then
      num_wilson=.true.
      if (rank .eq. 0) then
!         write(*,*) "You have activated NUM_WILSON, therefore the Wilson coordinate"
!         write(*,*) "transformation matrix will be calculated numerically (test mode)."
!         write(*,*)
      end if
   end if
end do


!
!     A: read the coordinates from file coord_def.inp
!       If the Wilson matrix shall be calculated in explore.x, the coordinates will
!       be read in from 'coord_wilson.inp'.
!
if (mode.eq.1) then
   if (filegeo .eq. "wilson") then
      write(*,*) "The set of internal coordinates for the calculation of a Wilson"
      write(*,*) " matrix will be read in from 'coord_wilson.inp'."
      allocate(coord_tmp(1000000,5))
      coord_tmp=0
      open (unit=165,file="coord_wilson.inp",status="old",iostat=readstat)
      if (readstat .ne. 0) then
         if (rank .eq. 0) then
            write(*,*) "The file coord_wilson.inp doesn´t exist!"
            call fatal
         end if
      end if
      i=1
   else 
      if (rank .eq. 0) then
         write(*,*) "The set of internal coordinates needed for DG-EVB/EVB-dQ coupling"
         write(*,*) "is read from the file 'coord_def.inp'"
      end if
      allocate(coord_tmp(1000000,5))
      coord_tmp=0
      open (unit=165,file="coord_def.inp",status="old",iostat=readstat)
      if (readstat .ne. 0) then
         if (rank .eq. 0) then
            write(*,*) "The file coord_def.inp doesn´t exist!"
            call fatal
         end if
      end if
      i=1
   end if
!
!   Use this temporary array for memory of coordinate types
!   Up to five integers: 2 for bonds, 3 for angles, 5 for torsions and out-of-plane 
!    (4 atoms plus one specifier of coordinate type)
!
   allocate(coord_types(50000))
   do
      read(165,"(A)",iostat=readstat) coord_line
      if (readstat .ne. 0)  exit

      read(coord_line,*,iostat=readstat) coord_tmp(i,5),coord_tmp(i,1),coord_tmp(i,2),coord_tmp(i,3),coord_tmp(i,4)
      if (readstat .ne. 0) then
         read(coord_line,*,iostat=readstat) coord_tmp(i,1),coord_tmp(i,2),coord_tmp(i,3)
         if (readstat .ne. 0) then
            read(coord_line,*,iostat=readstat) coord_tmp(i,1),coord_tmp(i,2)
            if (readstat .ne. 0)  then
               if (rank .eq. 0) then
                  write(*,*) "The format of a line in the file coord_def.inp seems to be corrupted!"
                  call fatal
               end if
            end if
            nbonds=nbonds+1
            coord_types(i)=1
         else
            nangles=nangles+1
            coord_types(i)=2
         end if
      else
!
!     dihedral angle: if the specifier is set to zero!
!
         if (coord_tmp(i,5) .eq. 0) then
            ndiheds=ndiheds+1
            coord_types(i)=3
!
!     out of plane angle: if the specifier is set to one!
!
         else if (coord_tmp(i,5) .eq. 1) then
            noops=noops+1
            coord_types(i)=4
         end if
      end if
      i=i+1
   end do
   if (i.eq.1) then
      write(*,*) "The file 'coord_def.inp' with internal coordinates is corrupted!"
      call fatal
   end if
!
!    only if the coordinates are not used for the rpmd.x program..
!
   if ((i.eq.2) .and. (.not. use_rpmd)) then
      write(*,*) "At least two internal coordinates are needed in 'coord_def.inp'!"
      call fatal
   end if
   nat6=i-1
!
!   Check if all array indices for internal coordinates are valid
!
   do j=1,500
      do k=1,4
         if (coord_tmp(j,k) .gt. natoms) then
            write(*,*) "There are higher indices than the number of atoms!"
            call fatal
         end if
      end do
   end do
!
!    Now define the final array of coordinate definitions
!
!     elements of the coord_def array:
!     1: bondlenghts
!     2: bend angles
!     3: dihedral angles
!     4: out-of-plane angles     
!  
   allocate(coord_def(nat6,5))
   coord_def=0
   do j=1,nat6
      if (coord_types(j).eq.1) then
         coord_def(j,1)=1
         coord_def(j,2:5)=coord_tmp(j,1:4)
      else if (coord_types(j).eq.2) then
         coord_def(j,1)=2
         coord_def(j,2:5)=coord_tmp(j,1:4)
      else if (coord_types(j).eq.3) then
         coord_def(j,1)=3
         coord_def(j,2:5)=coord_tmp(j,1:4)
      else if (coord_types(j).eq.4) then
         coord_def(j,1)=4
         coord_def(j,2:5)=coord_tmp(j,1:4)
      end if
   end do
   deallocate(coord_types)
   close(165)
   deallocate(coord_tmp)
!
!     B: define the coordinate set with num coord elements on the fly 
!
else if (mode .eq. 2) then
   if (rank .eq. 0) then
      write(*,*) "The set of internal coordinates needed for DG-EVB/EVB-dQ coupling"
      write(*,*) "is defined from the structure automatically!"
      write(*,'(a,i4,a)') " The",num_coord," coordinates with largest fluctuations along the"
      write(*,*)  "reaction path are chosen."
   end if
   o=0
!
!     First, read in a geometry in the middle of the reaction path
!     this is the reference for the building of possible internal 
!     coordinates
!

   nat=natoms
   nat3=3*nat
   allocate(xyz_path(3,nat,points))
   allocate(path_xyz(nat*3,points))
   allocate(coord_test(3,nat))
   open (unit=127,file=filegeo,status='old')
   do i=1,points
      call next_geo(coord,natoms,127,has_next)
      xyz_path(:,:,i)=coord
      if (.not.has_next) exit
   end do
   ihalf=i/2
   coord_test=xyz_path(:,:,ihalf)
!
!     calculate all internal coordinates
!
!     calculate bondlengths (max. 6 bohr)
!
   allocate(coord_tmp(50000,4))
   do i=1,natoms
      do j=1,i-1
         if (dist(i,j,coord_test) .lt. 6.0) then
            o=o+1
            coord_tmp(o,1)=i
            coord_tmp(o,2)=j
         end if
      end do
   end do
   nbonds=o  ! number of bonds
!
!     calculate bend angles (bonds max. 6 bohr)
!
   do i=1,natoms
      do j=1,i-1
         do k=1,j-1
            if ((dist(i,j,coord_test).lt.4.0) .and. (dist(j,k,coord_test).lt.4.0)) then
            o=o+1
            coord_tmp(o,1)=i
            coord_tmp(o,2)=j
            coord_tmp(o,3)=k
            end if
         end do
      end do
   end do
   nangles=o-nbonds ! number of angles
!
!    calculate torsional angles (bonds max. 6 bohr)
!
   do i=1,natoms
      do j=1,i-1
         do k=1,j-1
            do l=1,k-1
               if ((dist(i,j,coord_test).lt.4.0) .and. &
                  & (dist(j,k,coord_test).lt.4.0) .and. &
                  & (dist(k,l,coord_test).lt.4.0)) then
                  o=o+1
                  coord_tmp(o,1)=i
                  coord_tmp(o,2)=j
                  coord_tmp(o,3)=k
                  coord_tmp(o,4)=l
               end if
            end do
         end do
      end do
   end do
   ndiheds=o-nbonds-nangles
!
!    initialize the array with all internal coordinates (index for type)
!
   allocate(coord_def(o,5))

!
!     elements of the coord_def array:
!     1: bondlenghts
!     2: bend angles
!     3: dihedral angles
!          
      do j=1,nbonds
         coord_def(j,1)=1
         coord_def(j,2:5)=coord_tmp(j,:)
      end do
      do j=1,nangles
         coord_def(j+nbonds,1)=2
         coord_def(j+nbonds,2:5)=coord_tmp(j+nbonds,:)
      end do
      do j=1,ndiheds
         coord_def(j+nbonds+nangles,1)=3
         coord_def(j+nbonds+nangles,2:5)=&
                   & coord_tmp(j+nbonds+nangles,:)
      end do
      nat6=nbonds+nangles+ndiheds
   close(165)
   deallocate(coord_tmp)
!
!     now calculate internal variables for all structures on the reference 
!         reaction path
!  

!     short forms for variables
   allocate(path_int(nat6,points))
   allocate(int_path(nat6))
   allocate(rmsd_path_int(nat6))
   open (unit=127,file=filegeo,status='old')
   do i=1,points
      coord=xyz_path(:,:,i)
      call xyz_2int(coord,int_path,nat)
      if (.not.has_next) exit
      do j=1,nat6
         path_int(j,i)=int_path(j)
      end do
   end do
   close(127)
!
!     calculate RMS change for all internal variables and all steps
!
   rmsd_path_int=0
   do i=2,points
      do j=1,nat6
         rmsd_path_int(j)=rmsd_path_int(j)+((path_int(j,i)-path_int(j,i-1)))**2
      end do
   end do
!
!     take the absolute value for total array
!
   rmsd_path_int=abs(rmsd_path_int)
!
!     sort the entries in ascending number
!  
   allocate(temp2(points))
   allocate(temp3(points))
   do i=1,nat6-1
      start=i
      minimum=rmsd_path_int(start)
      location=start
      do j=start+1,nat6
         if (rmsd_path_int(j) .ge. minimum) then
            minimum=rmsd_path_int(j)
            location=j
         end if
      end do
      temp=rmsd_path_int(i)
      temp2=coord_def(i,:)
      temp3=path_int(i,:)

      rmsd_path_int(i)=rmsd_path_int(location)
      coord_def(i,:)=coord_def(location,:)
      path_int(i,:)=path_int(location,:)

      rmsd_path_int(location)=temp
      coord_def(location,:)=temp2
      path_int(location,:)=temp3
   end do
   deallocate(temp2,temp3)
!
!     use only the number of coordinates defined via num_coord
!     reshape the size of the coord_def array!
!
   nat6=num_coord
   allocate(coord_reshape(nat6,5))
   do i=1,nat6
      coord_reshape(i,:)=coord_def(i,:)
   end do
   deallocate(coord_def)
   allocate(coord_def(nat6,5))
   coord_def=coord_reshape
   deallocate(coord_reshape)
   deallocate(path_int,int_path)
!
!    C: define as many internal coordinates as possible (maximal number: 3N-6)
!

else if (mode .eq. 3) then
   if (rank .eq. 0) then
      write(*,*) "Coordinate initialization - mode 3: The bond orders for both"
      write(*,*) "QMDFF minima will be compared and filled with angle coords." 
      write(*,*) "The maximum number will be 3N-6 coordinates! (3N-7 for planar)"
      write(15,*) "Coordinate initialization - mode 3: The bond orders for both"
      write(15,*) "QMDFF minima will be compared and filled with angle coords."
      write(15,*) "The maximum number will be 3N-6 coordinates! (3N-7 for planar)"
   end if

   o=0
!
!     allocate test array with large size for fill in of coordinates
!
   allocate(coord_tmp(50000,4))
   allocate(bond_add(10000,4))
   coord_tmp=0
!
!     Determine element numbers of the systems atoms
!
   do i=1,natoms
      atind(i)=at(i)
   end do
! 
!     Read in the structure of the TS guess for coordinate definitions
!
   open (unit=127,file="ts.xyz",status='old')
   do i=1,points
      call next_geo(coord,natoms,127,has_next)
   end do
!
!     Convert coordinates to bohr
!
   coord=coord/bohr
   nbonds_old=0
!
!     set tolerance factor for bondlengths (if a bond shall still be defined..)
!
   bondlentol=1.2d0

   if (rank .eq. 0) then
      write(15,*) "Begin generation of natural bonding coordinate set."
      write(15,*) "Based on the QMDFF reference bond orders, covalent bonds"
      write(15,*) "will be picked as bondlengths, then all possible angles"
      write(15,*) "Between neighbored bonds are generated."
      write(15,*) "If more than 6N-6 coordinates are generated, remove all"
      write(15,*) "additional bonds."
   end if
!
!     Compare found bonds for QMDFFs of both minima! Each bond that appears 
!     at least once will be taken as bondlength coordinate
!     Take only bonds that are not too long! (i.e., covalent bondlengths)
!     Further, store all possible longer bonds in array bond_add for later fill up of 
!     coordinate definition array (if needed)
!
      
   p=0
   do i=1,nbond
      if (vbond(1,i) .lt. (rad(atind(bond(1,i)))+rad(atind(bond(2,i))))*bondlentol) then
         o=o+1
         coord_tmp(o,1:2)=bond(:,i)
      else 
         p=p+1
         bond_add(p,1:2)=bond(:,i) 
      end if
   end do
!
!     check for double appearences of bond in bond order matrices
!     Fill larger bond_add array with all bonds that are not doubled
!
   do i=1,nbond_two
      double_bond=0
      do j=1,nbond
         if ((bond_two(1,i) .eq. bond(1,j)) .and. (bond_two(2,i) .eq. bond(2,j))) then
            double_bond=1
         else
            if ((bond_two(1,i) .eq. bond(2,j)) .and. (bond_two(2,i) .eq. bond(1,j))) then 
               double_bond=1
            else
            end if
         end if
      end do
      if (double_bond .eq. 0) then
         if (vbond_two(1,i) .lt. (rad(atind(bond_two(1,i)))+rad(atind(bond_two(2,i))))*&
                    & bondlentol) then
            o=o+1
            coord_tmp(o,1:2)=bond_two(:,i)
         else 
            p=p+1
            bond_add(p,1:2)=bond_two(:,i)
         end if
      end if
   end do
!
!     Determine total number of covalent bonds and of longer (hypothetical) bonds 
!     to be added
!
   nbonds=o
   nbonds_add=p
!
!     Now, take possible angle coordinates (connecting neighbored bonds) until
!     3N-6 or 3N-7 coordinates are defined 
!
!     Three atoms: only two coordinates! (both bondlengths)
!
   o=0
   if (natoms .eq. 3) then
      nangles=0   
   else 
      nangles=0
      do i=1,nbonds
         do j=i,nbonds
            if (i .ne. j) then
               atm1=coord_tmp(i,1)
               atm2=coord_tmp(i,2)
               atm3=coord_tmp(j,1)
               atm4=coord_tmp(j,2)
               if ((atm1 .ne. atm2) .and. (atm1 .ne. atm3) .and. (atm1 .ne. atm4) .and. &
                   & (atm2 .ne. atm3) .and. (atm2 .ne. atm4) .and. (atm3 .ne. atm4)) then

               else
                  o=o+1
                  if (atm1 .eq. atm2) then
                     coord_tmp(nbonds+o,1)=atm3
                     coord_tmp(nbonds+o,2)=atm1
                     coord_tmp(nbonds+o,3)=atm4
                  else if (atm1 .eq. atm3) then
                     coord_tmp(nbonds+o,1)=atm2
                     coord_tmp(nbonds+o,2)=atm1
                     coord_tmp(nbonds+o,3)=atm4
                  else if (atm1 .eq. atm4) then
                     coord_tmp(nbonds+o,1)=atm2
                     coord_tmp(nbonds+o,2)=atm1
                     coord_tmp(nbonds+o,3)=atm3
                  else if (atm2 .eq. atm3) then
                     coord_tmp(nbonds+o,1)=atm1
                     coord_tmp(nbonds+o,2)=atm2
                     coord_tmp(nbonds+o,3)=atm4
                  else if (atm2 .eq. atm4) then
                     coord_tmp(nbonds+o,1)=atm1
                     coord_tmp(nbonds+o,2)=atm2
                     coord_tmp(nbonds+o,3)=atm3
                  else if (atm3 .eq. atm4) then
                     coord_tmp(nbonds+o,1)=atm1
                     coord_tmp(nbonds+o,2)=atm3
                     coord_tmp(nbonds+o,3)=atm2
                  end if
               end if
            end if
         end do
      end do
      nangles=o
   end if
!
!    Check if too many coordinates were determined and remove redundant angles
!   
   if (nangles + nbonds .gt. 3*natoms-6) then
      if (rank .eq. 0) then
         write(15,*) "Too many coordinates were defined for this system!"
         write(15,*) "3N-6 coordinates:",3*natoms-6,", actual coordinates:",nangles+nbonds
         write(15,*) "Redundant angles will be removed.."
      end if
      nangles=nangles-((nangles+nbonds)-(3*natoms-6))
   end if
!
!    initialize the array with all internal coordinates (index for type)
!
   allocate(coord_def(nbonds+nangles,5))
   coord_def=0
   
!
!     elements of the coord_def array:
!     1: bondlenghts
!     2: bend angles
!     3: dihedral angles
!      
   if (rank .eq. 0) then   
      write(15,*) 
      write(15,*) "Summary of pre-generated coordinate set:"
      write(15,*) "3*N-6 coordinates :",3*natoms-6
      write(15,*) "Number of bonds (nbonds):",nbonds
      write(15,*) "Number of angles (nangles):",nangles
      write(15,*) "Number of internals (nint):",nbonds+nangles 
      write(15,*)
   end if
!
!     The chosen set: check if the Wilson matrix has any singular values 
!     (nearby zero), then, check, what coordinate(s) are resposible for that
!     behavior
! 
   nat6=nbonds+nangles
   write(15,*) "The current full coordinate set:" 
   write(15,*) "----- Bonds -----"
   do j=1,nbonds
      coord_def(j,1)=1
      coord_def(j,2:5)=coord_tmp(j,:)
      write(15,'(a,i5,a,2i5)') "  * ",j,": ",coord_def(j,2:3)
   end do
   write(15,*) "----- Angles -----"
   do j=1,nangles
      coord_def(j+nbonds,1)=2
      coord_def(j+nbonds,2:5)=coord_tmp(j+nbonds,:)
      write(15,'(a,i5,a,3i5)') "  * ",j,": ",coord_def(j+nbonds,2:4)
   end do
   write(15,*) 
   if (rank .eq. 0) write(*,*)
   if (rank .eq. 0) write(*,*) "Check sanity of chosen coordinate set:"
!
!     Calculate the internal gradient for a given uniform cartesian gradient 
!     (should be numerically zero). If not, the Wilson matrix is broken
!
   allocate(int_test(nat6),grad_int(nat6))
   call xyz_2int(coord,int_test,natoms)
   grad_xyz=0.001d0
   call grad2int(coord,int_test,grad_int,grad_xyz)
!   allocate(hess_int(nat6,nat6))
!   hess_xyz=0.001d0
!   call hess2int(coord,int_test,hess_int,hess_xyz,grad_int)

!   write(*,*) "hess_sum",sum(hess_int**2)
!
!     If the actual coordinate set is bad, remove coordinates in order 
!     to check which are responsible for that
!
   if (sum(grad_int**2) .gt. 1D-10) then
      if (rank .eq. 0) then
         write(15,*) " The Wilson matrix is singular for this set!"
         write(15,*) " ---> Which coordinate is responsible for that?"
      end if
!
!     Check only angle coordinates in order to save calculation efforts!
!     Set full number of internal coordinates to bonds+angles
!
      nat6_check=nangles
      nat6=nangles+nbonds
!
!     Outer loop: decrease the number of used coordinates (or: increase the 
!     number of "dead" coordinates)
!
      change_num: do k=1,nat6_check
         if (rank .eq. 0) then
            if (k .eq. 1) then
               write(15,*) "************"
               write(15,*) "One coordinate will be dismissed per round."
            else 
               write(15,*) "************"
               write(15,'(a,i1,a)') " ",k," coordinates will be dismissed per round."
            end if    
         end if
         nat6=(nangles+nbonds)-k
!
!     Number of combinations in total: n!/((n-k)!*k!)
!       (n: total number of coords, k: dead coords)
!
         if (rank .eq. 0) write(*,*) "Remaining number of internal coordinates:",nat6
         ndead=nat6_check-(nat6-nbonds)
         ncombi=int(fact(nat6_check)/(fact(nat6_check-ndead)*fact(ndead)))
         if (rank .eq. 0) write(*,*) "Number of combinations to be checked:",ncombi
!
!     (Re-) allocate array with actual dead coordinates 
!
         if (allocated(outtakes)) deallocate(outtakes)
         if (.not. allocated(outtakes)) allocate(outtakes(ndead,ncombi))
!
!     For 1 to 5 dead coordinates: manually fill the array
!
         o=0
!
!     One coordinate is bad
!
         if (ndead .eq. 1) then
            do i1=1,nat6_check
               o=o+1
               outtakes(1,o)=i1
            end do
!
!     At least two coordinates are bad
!
         else if (ndead .eq. 2) then
            do i1=1,nat6_check
               do i2=i1,nat6_check
                  if (i1 .ne. i2) then
                     o=o+1
                     outtakes(1,o)=i1
                     outtakes(2,o)=i2
                   end if                           
               end do
            end do
!
!     At least three coordinates are bad 
!
         else if (ndead .eq. 3) then
            do i1=1,nat6_check
               do i2=i1,nat6_check
                  do i3=i2,nat6_check
                     if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3)) then
                        o=o+1
                        outtakes(1,o)=i1
                        outtakes(2,o)=i2
                        outtakes(3,o)=i3
                     end if
                  end do
               end do
            end do
!
!     At least four coordinates are bad
!
         else if (ndead .eq. 4) then
            do i1=1,nat6_check
               do i2=i1,nat6_check
                  do i3=i2,nat6_check
                     do i4=i3,nat6_check
                        if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3) &
                             & .and. (i1 .ne. i4) .and. (i2 .ne. i4) .and. & 
                             & (i3 .ne. i4)) then
                           o=o+1
                           outtakes(1,o)=i1
                           outtakes(2,o)=i2
                           outtakes(3,o)=i3
                           outtakes(4,o)=i4
                        end if
                     end do
                  end do
               end do
            end do
!
!     At least five coordinates are bad
!
         else if (ndead .eq. 5) then
            do i1=1,nat6_check
               do i2=i1,nat6_check
                  do i3=i2,nat6_check
                     do i4=i3,nat6_check
                        do i5=i4,nat6_check
                           if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3) &
                                & .and. (i1 .ne. i4) .and. (i2 .ne. i4) .and. &
                                & (i3 .ne. i4) .and. (i1 .ne. i5) .and. (i2 .ne. i5) .and. &
                                & (i3 .ne. i5) .and. (i4 .ne. i5) ) then
                              o=o+1
                              outtakes(1,o)=i1
                              outtakes(2,o)=i2
                              outtakes(3,o)=i3
                              outtakes(4,o)=i4
                              outtakes(5,o)=i5
                           end if
                        end do
                     end do
                  end do
               end do
            end do
!
!     If more than five are bad, concel the calculations (too many combinations...)
!
         else  
            write(*,*) "Too many coordinates lead to singular values in the Wilson matrix!"
            write(*,*) "Automatic generation of internals failed! Please do it manually"
            write(*,*) "and inform Julien Steffen that this error occured!"
            write(15,*) "Too many coordinates lead to singular values in the Wilson matrix!"
            write(15,*) "Automatic generation of internals failed! Please do it manually"
            write(15,*) "and inform Julien Steffen that this error occured!"
            call fatal
         end if
!
!     Go into inner loop: check all precalculated combinations of dismissed 
!     coordinates 
!     Reallocate local arrays for calculations         
!
         change_set: do i=1,ncombi
            if (allocated(coord_def)) deallocate(coord_def)
            if (allocated(int_test)) deallocate(int_test)
            if (allocated(grad_int)) deallocate(grad_int)
            if (.not. allocated(coord_def)) allocate(coord_def(nat6,5))
            if (.not. allocated(int_test)) allocate(int_test(nat6))
            if (.not. allocated(grad_int)) allocate(grad_int(nat6))                  
            coord_def=0
            o=0
!
!     Fill all initial bondlength coordinates into the coordinate array
!
            do j=1,nbonds
               o=o+1
               coord_def(o,1)=1
               coord_def(o,2:5)=coord_tmp(j,:)
            end do
!
!     Loop over all angle coordinates and remove all combinations of them in order 
!     test for validiy
!                  
            do j=1,nangles
               ignore=0
               do l=1,ndead
                  if (j .eq. outtakes(l,i)) then
                     ignore=1
                  end if
               end do
               if (ignore .eq. 0) then
                  o=o+1
                  coord_def(o,1)=2
                  coord_def(o,2:5)=coord_tmp(j+nbonds,:)
               end if
            end do
            call xyz_2int(coord,int_test,natoms)
            grad_xyz=0.001d0
            call grad2int(coord,int_test,grad_int,grad_xyz) 
!
!     If finally a good set with good behavior of the Wilson matrix was found,
!     exit both loops!
!                       
            if (sum(grad_int**2) .le. 1D-10) then
               write(15,*)
               write(15,*) "The bad coordinate(s) were found! They are:"
               write(15,*) "----- Angles -----"
               write(*,*) "the bad coordinates were:"
               do j=1,ndead
                  write(*,*) j,":",outtakes(j,i)+nbonds
                  write(15,'(a,i5,a,3i5)') "  * ",j,": ",coord_tmp(outtakes(j,i)+nbonds,1:3)
               end do
               write(15,*) 
               exit change_num
            end if
         end do change_set
         write(15,*)  
         write(15,*) "No coordinate(s) were found to remove the error!"
         write(15,*) "One more will be removed in the next round."
      end do change_num
   else 
      if (rank .eq. 0) write(*,*) " ... done! The coordinate set works well!"
   end if

!
!     The bad coordinates were extracted. If the remaining set is smaller than 3N-6 coordinates 
!     (and contains more than three atoms), additional bond coordinates will be added 
!     until a full good set was found. For that manner, sanity checks will be done for 
!     each added coordinate and the whole set at this moment!
!
   nat6_base=nat6
   if (.not. allocated(coord_base)) allocate(coord_base(nat6,5))
   fillup_done=.false.
   if ((3*natoms-6 .gt. nat6) .and. (natoms .gt. 3)) then
      if (rank .eq. 0) then
         write(15,*) "Too few coordinates are in the set!"
         write(15,'(a,i4,a,i4,a)') " There are ",3*natoms-6," coordinates needed and ", &
                     & nat6," currently present." 
         write(*,*) "Too few coordinates are in the set!"
         write(*,'(a,i4,a,i4,a)') " There are ",3*natoms-6," coordinates needed and ", &
                     & nat6," currently present."
      end if
!
!     The number of added coordinates 
!    
      nadd=(3*natoms-6)-nat6

      if (rank .eq. 0) then
         write(*,*) "Therefore,",nadd," internal coordinates will be added from the remaining"
         write(*,*) "pool of longer bonds in the QMDFF bond order matrices."
         write(15,*) "Therefore,",nadd," internal coordinates will be added from the remaining"
         write(15,*) "pool of longer bonds in the QMDFF bond order matrices."

         if (nbonds_add .eq. 0) then
            write(15,*) 
            write(15,*) "!!!Unfortunately, no further bond lengths remained!!!!"
         else 
            write(15,*) "This pool contains to following bonds:"
            do i=1,nbonds_add
               write(15,'(a,i5,a,2i5)') "  * ",i,": ",bond_add(i,1:2)
            end do
            write(15,*)  
         end if
      end if
!
!     Store the found set from above into a temporary array
!
      if (allocated(coord_base)) deallocate(coord_base)
      if (.not. allocated(coord_base)) allocate(coord_base(nat6,5))
      do i=1,nat6
         coord_base(i,:)=coord_def(i,:)
      end do
!
!     Similar to above, set up all possible combinations of coordinates to be added for 
!     upfollowing testing of these added combinations 
!
!     Number of possible addition combinations: n!/((n-k)!*k!) 
!     (n: total number of addible coords, k: added coords)
!
      ncombi=int(fact(nbonds_add)/(fact(nbonds_add-nadd)*fact(nadd)))  
      if (rank .eq. 0) then
         if (ncombi .gt. 0) then
            write(15,'(a,i8,a,i1,a)') " Now, ",ncombi," combinations of ",nadd, &
                    & " coordinates from this set will"
            write(15,*) "be added to the incomplete set until the Wilson matrix remains wellbehaved."
         end if 
      end if
!
!
!     For 1 to 5 dead coordinates: manually fill the array
!
      o=0
      if (.not. allocated(fillup)) allocate(fillup(nadd,ncombi))
!
!     One coordinate need to be added
!
      if (nadd .eq. 1) then
         do i1=1,nbonds_add
            o=o+1
            fillup(1,o)=i1
         end do
!
!     Two coordinates need to be added
!
      else if (nadd .eq. 2) then
         do i1=1,nbonds_add
            do i2=i1,nbonds_add
               if (i1 .ne. i2) then
                  o=o+1
                  fillup(1,o)=i1
                  fillup(2,o)=i2
                end if
            end do
         end do
!
!     Three coordinates need to be added
!
      else if (nadd .eq. 3) then
         do i1=1,nbonds_add
            do i2=i1,nbonds_add
               do i3=i2,nbonds_add
                  if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3)) then
                     o=o+1
                     fillup(1,o)=i1
                     fillup(2,o)=i2
                     fillup(3,o)=i3
                  end if
               end do
            end do
         end do
!
!     Four coordinates need to be added
!
      else if (nadd .eq. 4) then
         do i1=1,nbonds_add
            do i2=i1,nbonds_add
               do i3=i2,nbonds_add
                  do i4=i3,nbonds_add
                     if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3) &
                          & .and. (i1 .ne. i4) .and. (i2 .ne. i4) .and. &
                          & (i3 .ne. i4)) then
                        o=o+1
                        fillup(1,o)=i1
                        fillup(2,o)=i2
                        fillup(3,o)=i3
                        fillup(4,o)=i4
                     end if
                  end do
               end do
            end do
         end do
!
!     Five coordinates need to be added
!
      else if (nadd .eq. 5) then
         do i1=1,nbonds_add
            do i2=i1,nbonds_add
               do i3=i2,nbonds_add
                  do i4=i3,nbonds_add
                     do i5=i4,nbonds_add
                        if ((i1 .ne. i2) .and. (i1 .ne. i3) .and. (i2 .ne. i3) &
                             & .and. (i1 .ne. i4) .and. (i2 .ne. i4) .and. &
                             & (i3 .ne. i4) .and. (i1 .ne. i5) .and. (i2 .ne. i5) .and. &
                             & (i3 .ne. i5) .and. (i4 .ne. i5) ) then
                           o=o+1
                           fillup(1,o)=i1
                           fillup(2,o)=i2
                           fillup(3,o)=i3
                           fillup(4,o)=i4
                           fillup(5,o)=i5
                        end if
                     end do
                  end do
               end do
            end do
         end do
!
!     If too many coordinates need to be added (and the scaling is too large)
!
      else 
         write(*,*) "There are too few internal coordinates in the set and too many that"
         write(*,*) "shall be added in order to reach the number of 3N-6 coordinates!"
         write(*,*) "Set up the internal coordinate set manually (coord_def.inp) and contact"
         write(*,*) "Julien Steffen to report this strange case to him!"
         write(15,*) "There are too few internal coordinates in the set and too many that"
         write(15,*) "shall be added in order to reach the number of 3N-6 coordinates!"
         write(15,*) "Set up the internal coordinate set manually (coord_def.inp) and contact"
         write(15,*) "Julien Steffen to report this strange case to him!"
         call fatal
      end if
!
!     Now, begin with the loop over all combinations and add these to the incomplete but 
!     wellbehaved coordinate set until a full coordinate set behaves also well!
!
      nat6_base=nat6
      nat6=3*natoms-6

      if (allocated(coord_def)) deallocate(coord_def)
      if (allocated(int_test)) deallocate(int_test)
      if (allocated(grad_int)) deallocate(grad_int)
      if (.not. allocated(coord_def)) allocate(coord_def(nat6,5))
      if (.not. allocated(int_test)) allocate(int_test(nat6))
      if (.not. allocated(grad_int)) allocate(grad_int(nat6))
      do i=1,ncombi
!
!     Refill the coordinate definition array with both the working base set and the 
!     filled up additions
!
         coord_def=0
         do j=1,nat6_base
            coord_def(j,:)=coord_base(j,:)
         end do
         do j=1,nadd
            coord_def(j+nat6_base,1)=1
            coord_def(j+nat6_base,2:3)=bond_add(fillup(j,i),1:2)
         end do
         call xyz_2int(coord,int_test,natoms)
         grad_xyz=0.001d0
         call grad2int(coord,int_test,grad_int,grad_xyz)
         if (sum(grad_int**2) .le. 1D-10) then
            if (rank .eq. 0) then
               write(15,*)
               write(15,*) "Good fill up bond coordinates were found! These are"
               do k=1,nadd
                  write(15,'(a,i5,a,2i5)') "  * ",k,": ",bond_add(fillup(k,i),1:2)
               end do
               write(*,*)
               write(*,*) "Good fill up bond coordinates were found! These are"
               do k=1,nadd
                  write(*,'(a,i5,a,2i5)') "  * ",k,": ",bond_add(fillup(k,i),1:2)
               end do
            end if
            fillup_done=.true.
            exit
         end if
         if (i .eq. ncombi) then
            if (rank .eq. 0) then
               write(15,*) 
               write(15,*) "Coordinate fillup wasn't successful, the Wilson matrix became"
               write(15,*) "singular for each combination!"
               write(15,*) "This does not have to be bad, the uncompleted set might still"
               write(15,*) "be useful!"
               write(*,*)
               write(*,*) "Coordinate fillup wasn't successful, the Wilson matrix became"
               write(*,*) "singular for each combination!"
               write(*,*) "This does not have to be bad, the uncompleted set might still"
               write(*,*) "be useful!"
            end if
            fillup_done=.false.
         end if
      end do
   else 
!
!     If there were already enough internal coordinates, store them also in temporary array!
!
      if (allocated(coord_base)) deallocate(coord_base)
      if (.not. allocated(coord_base)) allocate(coord_base(nat6,5))
      do i=1,nat6
         coord_base(i,:)=coord_def(i,:)
      end do
   end if
!
!     Do final evaluation and information about the used coordinate set
!     Further, write the found set of coordinates into file coord_def.inp!
!
   open(unit=99,file="coord_def.inp",status="unknown")
   write(15,*) 
   write(15,*) "---RESULT: The calculated internal coordinate set is:----"
   if (fillup_done) then
      write(15,*) "(The coordinate set has been filled up with backup coordinates.)"
      nat6=3*natoms-6
      do i=1,nat6  
         if (coord_def(i,1) .eq. 1) then
            write(15,'(a,i5,a,2i5)') "  * ",i,": ",coord_def(i,2:3)   
            write(99,*) coord_def(i,2:3) 
         else if (coord_def(i,1) .eq. 2) then
            write(15,'(a,i5,a,3i5)') "  * ",i,": ",coord_def(i,2:4)
            write(99,*) coord_def(i,2:4)
         end if
      end do
   else 
      write(15,*) "(No coordinate fillup has been done.)"
      nat6=nat6_base
      if (allocated(coord_def)) deallocate(coord_def)
      if (.not. allocated(coord_def)) allocate(coord_def(nat6,5))
      do i=1,nat6
         if (coord_base(i,1) .eq. 1) then
            write(15,'(a,i5,a,2i5)') "  * ",i,": ",coord_base(i,2:3)     
            write(99,*) coord_base(i,2:3)     
         else if (coord_base(i,1) .eq. 2) then
            write(15,'(a,i5,a,3i5)') "  * ",i,": ",coord_base(i,2:4)
            write(99,*) coord_base(i,2:4)
         end if
         coord_def(i,:)=coord_base(i,:)
      end do

   end if
!   call xyz_2int(coord,int_test,natoms)
!   grad_xyz=0.01d0
!   call grad2int(coord,int_test,grad_int,grad_xyz)
!   allocate(hess_int(nat6,nat6))
!   hess_xyz=0.01d0
!   call hess2int(coord,int_test,hess_int,hess_xyz,grad_int)

!   write(*,*) "hess_sum",sum(hess_int**2),sum(grad_int**2)
!   stop "GPgp"
   if (rank .eq. 0) then
      write(15,*) "-------------------------------------------"
      write(15,*) "The coordinate set was written to file coord_def.inp"
      write(15,*) 
      close(99)

      write(*,*)
      write(*,*) "The set of internal coordinates has been determined automatically!"
      write(*,*) "The set was written to file coord_def.inp for other calculations."
      write(*,*) "Look into rpmd.log for further details!"
      write(*,*)
   end if
      
end if

!
!     D: If the coordinate set shall be a simple distance matrix of all atoms!
!       e.g.:    1-2 1-3 1-4 
!                    2-3 2-4
!                        3-4
!        This might be done if all other options fail!
!
if (mode.eq.4) then
   nat6=natoms*(natoms-1)/2
   if (rank .eq. 0) then
      write(*,*) "The set of internal coordinates needed for DG-EVB/EVB-dQ coupling"
      write(*,*) "will be a simple distance matrix of all atoms (only bondlengths)"
      write(*,'(a,i6,a)') " These are in total ",nat6," bondlenghts."
      write(15,*) "The set of internal coordinates needed for DG-EVB/EVB-dQ coupling"
      write(15,*) "will be a simple distance matrix of all atoms (only bondlengths)"
      write(15,'(a,i6,a)') " These are in total ",nat6," bondlenghts."
   end if
!
!     Define the set of coordinates in the global array!
!

   allocate(coord_def(nat6,5))
   coord_def=0
   k=1
   do i=1,natoms  
      do j=i+1,natoms
         coord_def(k,1)=1
         coord_def(k,2)=i
         coord_def(k,3)=j
         k=k+1
      end do
   end do
!
!     Write the determined set of coordinates to file
!
   if (rank .eq. 0) then
      open(unit=99,file="coord_def.inp",status="unknown")
      write(15,*)
      write(15,*) "---RESULT: The calculated internal coordinate set is:----"
      write(15,*) "(The coordinate set has been filled up with backup coordinates.)"
      do i=1,nat6
         write(15,'(a,i5,a,2i5)') "  * ",i,": ",coord_def(i,2:3)
         write(99,*) coord_def(i,2:3)
      end do
      write(15,*) "-------------------------------------------"
      write(15,*) "The coordinate set was written to file distance_matrix.inp"
      write(*,*) "The coordinate set was written to file distance_matrix.inp"
      write(15,*)
      close(99)
   end if
end if

end subroutine init_int

