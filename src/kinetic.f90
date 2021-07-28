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
!     subroutine kinetic: computes the total kinetic energy and kinetic energy
!     contributions to the pressure tensor
!
!     part of EVB
!
subroutine kinetic (eksum,ekin)
use general
implicit none
integer::i,j,k
integer::start,stop
real(kind=8)::eksum,weigh
real(kind=8)::term,value
real(kind=8)::xr,yr,zr
real(kind=8)::x2,y2,z2
real(kind=8)::xcm,ycm,zcm
real(kind=8)::ekin(3,3)
real(kind=8)::inert(3,3)
!
!     zero out the total kinetic energy and its outer product
!
eksum = 0.0d0
do i = 1, 3
   do j = 1, 3
      ekin(j,i) = 0.0d0
   end do
end do
!
!     get the total kinetic energy and tensor for atomic sites
!
do i = 1, n
   term = 0.5d0 * mass(i) / convert
   do j = 1, 3
      do k = 1, 3
         value = term * v(j,i) * v(k,i)
         ekin(k,j) = ekin(k,j) + value
      end do
   end do
end do
eksum = ekin(1,1) + ekin(2,2) + ekin(3,3)
!
!     get the kinetic energy for Bussi-Parrinello barostat
!
if (isobaric .and. barostat.eq.'BUSSI') then
   term = dble(nfree) * gasconst * kelvin * taupres * taupres
   value = 0.5d0 * term * eta * eta
   do j = 1, 3
      ekin(j,j) = ekin(j,j) + value/3.0d0
   end do
   eksum = eksum + value
end if

return
end subroutine kinetic
