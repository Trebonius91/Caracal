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
implicit none
integer::ixyz,i,j,k,m,next
integer::freeunit
integer,allocatable::list(:) ! local array for tff-connections
logical::exist
character(len=20)::keyword
character(len=120)::xyzfile
character(len=120)::record
character(len=120)::string
character(len=240)::xyzline
integer::idummy ! dummy variable for line numbers..
integer::nmax
integer::tag2(1000) ! local array for atom numbering
logical::reorder

!
!     ask for the user specified input structure filename
!
exist= .false.
do i = 1, nkey
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

ixyz = freeunit ()
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
