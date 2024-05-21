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
!     subroutine rdv: Read in coordinates and atom types from VASP output
!
!     part of QMDFF
!
subroutine rdv(echo,fname,n,xyz,iat)
implicit real(kind=8) (a-h,o-z)
dimension::xyz(3,n),iat(n),xx(10)
real(kind=8)::bohr
character(len=128)::line
character(len=2)::a2
character(len=*)::fname
character(len=100)::dum
character(len=2)::potcar_syms(10)
integer::potcar_num
integer::ion_nums(10)
logical::echo
parameter (bohr=0.52917721092d0)

if (echo) then
   write(10,*) '========================='
   write(10,*) 'reading ... ',trim(fname)
   write(10,*) '========================='
end if

ich=142
potcar_num=0
open(unit=ich,file=fname)
do
   read(ich,'(a)') line
!
!    First read in the atomic symbols appearing in the system (and the POTCAR)
!
   if(index(line,' VRHFIN ') .ne.0) then 
      potcar_num=potcar_num+1
      read(ich,'(a)')line
      read(ich,'(a)')line
      read(ich,'(a)')line
      read(ich,*) dum,dum,dum,potcar_syms(potcar_num)
   end if
!
!    Then read the number of ion times and fill the element symbol array
!
   if(index(line,' ions per type ') .ne.0) then
      read(line,*) dum,dum,dum,dum,ion_nums(1:potcar_num)
   end if
!
!    Then read the positions of the atoms (cartesian coordinates in Angstrom)
!
   if(index(line,' position of ions in cartesian coordinates ') .ne.0) then
      do i=1,n
         read(ich,*) xyz(:,i)
      end do
      exit
   end if
end do
!
!     Fill the element names array from the POTCAR information
!
k=1
do i=1,potcar_num
   do j=1,ion_nums(i)      
      call elem(potcar_syms(i),iat(k))
      k=k+1
   end do
end do
!
!     convert coordinates into bohr!
!
xyz=xyz/bohr

close(ich)
return
end subroutine rdv
