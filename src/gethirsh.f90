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
!     subroutine gethirsh: Read in Hirshfeld charges from orca ouput file
!
!     part of QMDFF
!
subroutine gethirsh(n,chir,ok,fname2,fname3)
use qmdff
implicit none
integer::n,intsum,numi,idum
real(kind=8)::fdum
real(kind=8)::chir(n)
real(kind=8)::csum
logical::ok

character(len=80)::a,a80,header,attype,fname2,fname3
real(kind=8)::xx(20)
integer::nn,i
integer::readstat

!chir=0
ok=.false.
!
!     Read the Hirshfeld charges directly from the reference output
!     (no charges file needed)!
!     If the Hirshfeld-section appears more than once in the file
!     the latest appearing is printed out.
!
!
!     for orca output
!
if (software .eq. "O") then
   open(unit=142,file=fname2)
   chir=0.d0
   do
      read(142,'(a)',end=99) a80
      if(index(a80,'HIRSHFELD ANALYSIS').ne.0) then
         do i=1,4
            read(142,*) header
         end do
         do i=1,n
            read(142,*) numi,attype,chir(i)
            csum=csum+chir(i)
         end do
      end if
   end do
!
!     for CP2K output
!
else if (software .eq. "C") then
   open(unit=142,file=fname2)
   chir=0.d0
   do
      read(142,'(a)',end=99) a80
      if(index(a80,'Hirshfeld Charges').ne.0) then
         read(142,*) header
         do i=1,n
!
!    Two different formats possible (?)
!
            read(142,*,iostat=readstat) numi,attype,idum,fdum,fdum,fdum,fdum,chir(i)
            if (readstat .ne. 0) then 
               read(142,*,iostat=readstat) numi,attype,idum,fdum,fdum,chir(i)
            end if
            if (readstat .ne. 0) then
               write(*,*) "ERROR! The CP2K charge section has a wrong formate!"
               write(*,*) " Presumably you have used the wrong CP2K version, this program"
               write(*,*) " assumes that you are using CP2K4.1. The formate should be:"
               write(*,*) " Int.,Char.,Int.,Real,Real,Real,Real,Charge"
               call fatal
            end if
            csum=csum+chir(i)
         end do
      end if
   end do
!
!     for Gaussian output
!
else if (software .eq. "G") then
   open(unit=142,file=fname3)
   chir=0
   do
      read(142,'(a)',end=99) a80
      if(index(a80,'Hirshfeld charges, spin densities').ne.0) then
         do i=1,1
            read(142,*) header
         end do
         do i=1,n
            read(142,*) numi,attype,chir(i)
            csum=csum+chir(i)
         end do
      end if
   end do
 

end if
!     (never executed...)
intsum=NINT(csum)
if (abs(intsum-csum)/n.gt.2.d-2) stop "error in charges!"
write(10,*) 'sum of used atomic charges:',csum
99  close(142)
if(sum(abs(chir)).gt.1.d-6)ok=.true.

end subroutine gethirsh

