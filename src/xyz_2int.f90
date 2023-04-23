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
!     subroutine xyz2int: converts an array of cartesian coordinates to internal 
!     coordinates
!
!     part of EVB
!

subroutine xyz_2int(xyz5,internal,n)
use evb_mod
implicit none


real(kind=8)::xyz5(3,n),internal(nat6)
integer::i,n,l
real(kind=8)::x,y,z
real(kind=8)::ax,ay,az,bx,by,bz,cx,cy,cz,degree
real(kind=8)::sin_val,cos_val
real(kind=8)::dist,ang,dihed,oop

i=0
!
!     separate arrays for bondlenghts, angles, ooplanes  and dihedrals
!     store results into array for all internal coordinates
!
do i=1,nat6
   if (coord_def(i,1) .eq. 1) then
!
!  calculate bondlengths between chosen atoms
!
      internal(i)=dist(coord_def(i,2),coord_def(i,3),xyz5)
   else if (coord_def(i,1) .eq. 2) then
!
!  calculate bond angles between chosen atoms
!     
      internal(i)=ang(coord_def(i,2),coord_def(i,3),coord_def(i,4),xyz5)
   else if (coord_def(i,1) .eq. 3) then
!
!  calculate dihedral angles between chosen atoms
!
      internal(i)=dihed(coord_def(i,2),coord_def(i,3),coord_def(i,4),coord_def(i,5),xyz5)
   else if (coord_def(i,1) .eq. 4) then
!
!  calculate out-pf-plane angles between chosen atoms
!
      internal(i)=oop(coord_def(i,2),coord_def(i,3),coord_def(i,4),coord_def(i,5),xyz5)
   end if
end do

return
end subroutine xyz_2int

