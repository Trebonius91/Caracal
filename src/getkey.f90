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
!     subroutine getkey: Finds the keyfile and stores its content 
!        into subsequent line images 
!
!
!     Based on:
!     TINKER molecular modeling package
!     COPYRIGHT (C)  1990  by  Jay William Ponder   
!     All Rights Reserved 
!
!     part of EVB
!
subroutine getkey(rank)
use general
use qmdff

implicit none

include 'mpif.h'
integer::i,ikey
integer::next,length
integer::freeunit
integer::trimtext
logical::exist,header
character(len=20)::keyword
character(len=120)::comment
character(len=120)::record
character(len=120)::string
character(len=20)::blank
integer::ierr
integer::rank  ! the current MPI rank
!
!     obtain all command line arguments
!
!     initialize command line arguments as blank strings
!
narg = 0
blank = '                    '
do i = 0, maxarg
   arg(i) = blank//blank//blank
end do
!
!     get the number of arguments and store each in a string
!
narg = iargc ()
if (narg .gt. maxarg)  narg = maxarg
do i = 0, narg
   call getarg (i,arg(i))
end do

!if (details) then
!   narg=0
!end if

!
!     check for keyfile specified on command line
!     loop over all command line arguments
!
narg=1
exist = .false.

do i = 1, narg
   string = arg(i)
   keyfile = string 
   if (whichprog .eq. 1) exit
   if (keyfile .eq. "-details") exit
   inquire (file=keyfile,exist=exist)
   if (.not. exist) then
      if (rank .eq. 0) then
         write(*,*)
         write (iout,*) "The keyfile that was specified on command", &
              &  " line was not found!"
         call mpi_barrier(mpi_comm_world,ierr)
         call fatal
      end if
   end if
end do
if (.not. exist) then
   if (rank .eq. 0) then
      if (whichprog .ne. 1) then
         write(*,*)
         write(iout,*) "No keyfile was specified on command line!"
         call fatal
      end if
   end if
end if

!
!     read the keyfile and store it for latter use
!     into the keyline(:) array
!
if (exist) then
   ikey = freeunit ()
   open (unit=ikey,file=keyfile,status='old')
   rewind (unit=ikey)
   do while (.true.)
      read (ikey,'(A120)',err=40,end=40)  record
      nkey = nkey + 1
      keyline(nkey) = record
   end do
   40    continue
   close (unit=ikey)
end if
!
!     TEST: divide the nkey word through two!  
!
!nkey=nkey/2
return
end subroutine getkey

