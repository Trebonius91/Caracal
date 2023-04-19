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
!     subroutine rdc: Read in coordinates and atom types from CP2K output
!
!     part of QMDFF
!
subroutine rdc(echo,fname,n,xyz,iat)
implicit real(kind=8) (a-h,o-z)
dimension::xyz(3,n),iat(n),xx(10)
real(kind=8)::bohr
character(len=128)::line
character(len=2)::a2
character(len=*)::fname
logical::echo
parameter (bohr=0.52917721092d0)

if (echo) then
   write(10,*) '========================='
   write(10,*) 'reading ... ',trim(fname)
   write(10,*) '========================='
end if

ich=142
open(unit=ich,file=fname)
do
   read(ich,'(a)') line
   if(index(line,' MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom') &
       &  .ne.0) then
      read(ich,'(a)')line
      read(ich,'(a)')line
      read(ich,'(a)')line
      do i=1,n
!     element numbers are already included in this case!
         read(ich,*) dum,dum,a2,dum,xyz(1:3,i)
         call elem(a2,iat(i))
      end do
      exit
   end if
end do
!
!     convert coordinates into bohr!
!
xyz=xyz/bohr

close(ich)
return
end subroutine rdc
