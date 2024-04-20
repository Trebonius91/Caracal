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
!     subroutine getxyz: asks for a Cartesian coordinate file name,
!     then reads in the coordinates file
!
!     part of EVB
!
subroutine getxyz
use general
use evb_mod
use pbc_mod
implicit none
integer::ixyz,i,j,k,m,next,ind
integer::freeunit
integer::readstat
integer,allocatable::list(:) ! local array for tff-connections
logical::exist
character(len=20)::keyword
character(len=120)::xyzfile
character(len=120)::record
character(len=120)::string
character(len=240)::xyzline
character(len=1)::check_sel ! if T or F is noted in VASP POCARs
integer::idummy ! dummy variable for line numbers..
real(kind=8)::x_tmp,y_tmp,z_tmp  ! temporary coordinates of current atom
integer::nmax
integer::tag2(1000) ! local array for atom numbering
logical::reorder

!
!     ask for the user specified input structure filename
!
exist= .false.

do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'XYZSTART ') then
      call getword (record,xyzfile,next)
      call basefile (xyzfile)
      call suffix (xyzfile,'xyz','old')
      inquire(file=xyzfile,exist=exist)
   end if
end do
!
!     if it isn´t given in the command line
!
do while (.not. exist)
   write (iout,'(/," Enter Cartesian Coordinate File Name of Start Structure :  ",$)')
   read (input,'(a120)')  xyzfile
   call basefile (xyzfile)
   call suffix (xyzfile,'xyz','old')
   inquire (file=xyzfile,exist=exist)
end do
!
!     first open and then read the Cartesian coordinates file
!
!
!     If the name of the file is POSCAR, treat it as VASP input file 
!     and read it accordingly!
!
ixyz = freeunit ()
coord_vasp = .false.
if (trim(xyzfile) .eq. "POSCAR") then 
   write(*,*)
   write(*,*) "The geometry is given in file POSCAR. VASP format will be assumed!"
   coord_vasp = .true.
   open (unit=ixyz,file=xyzfile,status='old')
   rewind (unit=ixyz)
   read(ixyz,*,err=60,end=60) 
   read(ixyz,*,err=60,end=60) vasp_scale
   read(ixyz,*,err=60,end=60) vasp_a_vec
   read(ixyz,*,err=60,end=60) vasp_b_vec
   read(ixyz,*,err=60,end=60) vasp_c_vec
   vasp_names="XX"
!
!     Read the elements for all atoms
!
   read(ixyz,'(a)',err=60,end=60) string
   read(string,*,iostat=readstat) vasp_names
   vasp_numbers = 0
   read(ixyz,'(a)',err=60,end=60) string
   read(string,*,iostat=readstat) vasp_numbers
   ind=0
   nelems_vasp=0
   do i=1,20
      if (vasp_names(i) .eq. "XX") exit
      nelems_vasp=nelems_vasp+1
      do j=1,vasp_numbers(i)
         ind=ind+1
         name(ind)=vasp_names(i)
         call atommass(ind)
         call upcase(name(ind))
      end do    
   end do
   n=sum(vasp_numbers)
   natoms=n
   allocate(elem_index(natoms))
!
!    Store element indices 
!
   do i=1,natoms
      call elem(name(i),elem_index(i))
   end do
!
!    Check if selective dynamics or not
!
!    Check if POSCAR has direct or cartesian coordinates
!
   vasp_selective = .false.
   vasp_direct = .false.
   fix_atoms = .false.
   
   read(ixyz,'(a)',err=60,end=60) string 
   if (trim(string) .eq. "Selective dynamics" .or. & 
          & trim(string) .eq. "selective dynamics" .or. &
          & trim(string) .eq. "Selective" .or. &
          & trim(string) .eq. "selective" .or. &
          & trim(string) .eq. "Selective Dynamics" .or. &
          & trim(string) .eq. "Selective Dynamics") then
      vasp_selective = .true.
      fix_atoms = .true.
      allocate(fix_list(natoms))
   else
      if (trim(string) .eq. "direct" .or. trim(string) .eq. "Direct") then
         vasp_direct = .true.
      end if
   end if
   if (vasp_selective) then
      read(ixyz,'(a)',err=60,end=60) string
      if (trim(string) .eq. "direct" .or. trim(string) .eq. "Direct") then
         vasp_direct = .true.
      end if
   end if
!
!     Now read in the coordinates
!     If selective is activated, read in the first flag as well and 
!     decide if the atom is active or not (no coordinate distinction!)
!     Add the indices of the fixed atoms to the fix_list array
!
   ind=0   
   do i=1,natoms
      if (vasp_selective) then
         read(ixyz,*,err=60,end=60) x(i),y(i),z(i),check_sel
         if (check_sel .eq. "F") then
            ind=ind+1
            fix_list(ind)=i
         end if 
      else
         read(ixyz,*,err=60,end=60) x(i),y(i),z(i)
      end if
!
!     If direct coordinates are given, transform the coordinates to Angstrom!
!     
      if (vasp_direct) then
         x_tmp=x(i)
         y_tmp=y(i)
         z_tmp=z(i)
         x(i)=(x_tmp*vasp_a_vec(1)+y_tmp*vasp_b_vec(1)+z_tmp*vasp_c_vec(1))*vasp_scale
         y(i)=(x_tmp*vasp_a_vec(2)+y_tmp*vasp_b_vec(2)+z_tmp*vasp_c_vec(2))*vasp_scale
         z(i)=(x_tmp*vasp_a_vec(3)+y_tmp*vasp_b_vec(3)+z_tmp*vasp_c_vec(3))*vasp_scale         
      end if 
   end do
   if (vasp_selective) then
      fix_num=ind 
   end if

   close (unit=ixyz)
   return
end if

!
!    Read in a usual xyz file
!
open (unit=ixyz,file=xyzfile,status='old')
rewind (unit=ixyz)

read(ixyz,*,err=60,end=60) n
natoms=n
read(ixyz,*,err=60,end=60) 
allocate(elem_index(natoms))
do i=1,n
   read(ixyz,*,err=60,end=60) name(i),x(i),y(i),z(i)
   call atommass(i)
   call upcase(name(i))
!
!    Store element indices 
!
   call elem(name(i),elem_index(i))

end do

close (unit=ixyz)

return
60 continue
write(*,*) "Error in input file for cartesian coordinates!"
close (unit=ixyz)
call fatal

return
end subroutine getxyz
