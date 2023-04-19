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
!     subroutine ff_hb_two: energies and gradients of hydrogen/halogen
!         bond terms of the second QMDFF
!
!     part of QMDFF
!
subroutine ff_hb_two(n,at,xyz,eh,g)
use qmdff
implicit none
integer::n,at(*)
real(kind=8)::xyz(3,n),g(3,n),eh
integer::i1,i2,i,k,j
real(kind=8)::c1,c2,r,er,el,step,edum,dum

if (nhb_two.lt.1) return
!
!     loop over all hydrogen/halogen bonds:
!     first determine distance between donor/acceptor
!     count only if distance is smaller than 15 bohr
!
do k=1,nhb_two
   i1 =hb_two(1,k)
   i2 =hb_two(2,k)
   r=sqrt((xyz(1,i1)-xyz(1,i2))**2 &
     &   +(xyz(2,i1)-xyz(2,i2))**2 &
     &   +(xyz(3,i1)-xyz(3,i2))**2)
   if (r.gt.15.0d0) cycle
   i=hb_two(3,k)
!
!     if an H-atom takes part, calculate hydrogen bond
!
   if (at(i).eq.1) then
      c1 =vhb_two(1,k)
      c2 =vhb_two(2,k)
      call eabhag(n,i1,i2,i,xyz,c1,c2,eh,g)
   else
      c1 =vhb_two(1,k)
      call eabxag(n,i1,i2,i,xyz,c1,eh,g)
   end if
end do

return
end subroutine ff_hb_two
