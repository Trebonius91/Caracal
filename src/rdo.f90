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
!     subroutine rdo: Read in coordinates and atom types from orca output
!
!     part of QMDFF
!
subroutine rdo(echo,fname,outname,n,xyz,iat,check_coord)
implicit real(kind=8) (a-h,o-z)
dimension::xyz(3,n),iat(n),xx(10)
integer::idum,readstat
real(kind=8)::rdum
character(len=128)::line,adum
character(len=2)::a2
character(len=*)::fname,outname
logical::echo,check_coord,struc_exist

ich=142
if (check_coord) then
   if (echo) then
      write(10,*) '========================='
      write(10,*) 'reading ... ',trim(outname)
      write(10,*) '========================='
   end if

   ich=142
   struc_exist=.false.
   open(unit=ich,file=outname)
   do
      read(ich,'(a)',iostat=readstat) line
      if (readstat .ne. 0) then
         exit
      end if
      if(index(line,'CARTESIAN COORDINATES (A.U.)').ne.0) then
         read(ich,'(a)')line
         read(ich,'(a)')line
         do i=1,n
            read(ich,*) idum,a2,rdum,idum,rdum,xyz(1:3,i)
            call elem(a2,iat(i))
         end do
         struc_exist=.true.
        ! exit
      end if
   end do
   close(ich)
   if (.not. struc_exist) then
      write(*,*) "The file ",trim(outname)," seems to be corrupted!"
      call fatal
   end if
else 
   if (echo) then
      write(10,*) '========================='
      write(10,*) 'reading ... ',trim(fname)
      write(10,*) '========================='
   end if

   ich=142
   open(unit=ich,file=fname)
   do
      read(ich,'(a)') line
      if(index(line,'$atoms').ne.0) then
         read(ich,'(a)')line
         do i=1,n
            read(ich,*) a2,dum,xyz(1:3,i)
            call elem(a2,iat(i))
         end do
         exit
      end if
   end do
   close(ich)
end if

return
end subroutine rdo
