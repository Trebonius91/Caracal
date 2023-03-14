!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   EVB-QMDFF - RPMD molecular dynamics and rate constant calculations on
!               black-box generated potential energy surfaces
!
!   Copyright (c) 2021 by Julien Steffen (steffen@pctc.uni-kiel.de)
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
!     subroutine read_geo: Read in cartesian geometry from grad_hess.dat
!     file and convert to internals. The erngies of the actual structure 
!     is also read.
!
!     part of EVB
!
subroutine read_geo(atom_num,geo_int,geo_xyz,geo_xyz1,energy)
use evb_mod
implicit none
character(len=80) a80
integer::atom_num,k,i,maxcycle,rest,nat,j
real(kind=8),dimension(nat6)::geo_int
real(kind=8),dimension(atom_num*3)::geo_xyz
real(kind=8),dimension(3,atom_num)::geo_xyz1
real(kind=8)::energy,bohr
character(len=40)::names
logical::redun
redun=.false.
nat=atom_num
bohr=1.889725989
do
   read(ref_input_unit,'(a)') a80
   if(index(a80,'ENERGY ').ne.0) then
      read(a80,*) names,energy
   else if(index(a80,'GEOMETRY ').ne.0) then 
      maxcycle=(atom_num*3)/5
      do i=1,maxcycle
         k=(i-1)*5
         read(ref_input_unit,*) geo_xyz(k+1),geo_xyz(k+2),geo_xyz(k+3),geo_xyz(k+4),geo_xyz(k+5)
      end do
      rest=(atom_num*3)-maxcycle*5
      k=maxcycle*5
      if (rest .eq. 4) then
         read(ref_input_unit,*) geo_xyz(k+1),geo_xyz(k+2),geo_xyz(k+3),geo_xyz(k+4)
      else if (rest .eq. 3) then
         read(ref_input_unit,*) geo_xyz(k+1),geo_xyz(k+2),geo_xyz(k+3)
      else if (rest .eq. 2) then
         read(ref_input_unit,*) geo_xyz(k+1),geo_xyz(k+2)
      else if (rest .eq. 1) then
         read(ref_input_unit,*) geo_xyz(k+1)
      end if
      goto 30
   end if
end do
30 do i=1,nat
   do j=1,3
      geo_xyz1(j,i)=geo_xyz((i-1)*3+j)!/bohr
   end do
end do
call xyz_2int(geo_xyz1,geo_int,nat)

!write(*,*) geo_int
!call grad2int(geo_xyz1,geo_int)
return
end subroutine read_geo
