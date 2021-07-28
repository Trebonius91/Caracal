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
   if (keyword(1:11) .eq. 'XYZFILE ') then
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
!     If the Müller-Brown suface is used, read in only one line 
!
if (mueller_brown) then
   open (unit=ixyz,file=xyzfile,status='old')
   read(ixyz,*,err=60,end=60) x(1),y(1)
   z(1)=0.d0
   return
end if
!
!     first open and then read the Cartesian coordinates file
!

ixyz = freeunit ()
open (unit=ixyz,file=xyzfile,status='old')
rewind (unit=ixyz)

read(ixyz,*,err=60,end=60) n
natoms=n
read(ixyz,*,err=60,end=60) 
do i=1,n
   read(ixyz,*,err=60,end=60) name(i),x(i),y(i),z(i)
   call atommass(i)
end do

close (unit=ixyz)

return
60 continue
write(*,*) "Error in input file for cartesian coordinates!"
close (unit=ixyz)
call fatal

return
end subroutine getxyz
