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
!     rdo0: Read in name of file and number of atoms from orca output
!
!     part of QMDFF
!
subroutine rdo0(fname,outname,n,check_coord)
implicit real(kind=8) (a-h,o-z)
character(len=128)::line,adum
character(len=*)::fname,outname
logical::ex,check_coord
!
!     If the check_coord mode is activated, no Hessian is required, and 
!      therefore, no name.hess file must be there, name.out will be opened 
!      instead
!
n=0
ich=142

if (check_coord) then
   inquire(file=outname,exist=ex)
   if(.not.ex) return
   open(unit=ich,file=outname)
   do
      read(ich,'(a)') line
      if (index(line,'Number of atoms').ne.0) then
         read(line,*) adum,adum,adum,adum,n
         close(ich)
         return
      end if
   end do
   close(ich)
!
!     Normal mode: read in the name.hess file
!
else 
   inquire(file=fname,exist=ex)
   if(.not.ex) return
   open(unit=ich,file=fname)
   do
      read(ich,'(a)') line
      if (index(line,'$atoms').ne.0) then
         read(ich,*) n
         close(ich)
         return
      end if
   end do
   close(ich)
end if

return

end subroutine rdo0
