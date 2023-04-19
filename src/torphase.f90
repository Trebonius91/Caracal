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
!     subroutine "torphase" sets the n-fold amplitude and phase values
!     for each torsion via sorting of the input parameters
!
!     part of others
!
subroutine torphase (ft,vt,st)
implicit none
integer::i,k
integer::ft(*)
real(kind=8)::ampli(6)
real(kind=8)::phase(6)
real(kind=8)::vt(*),st(*)
!
!     copy the input fold, amplitude and phase angles
!
do i = 1, 6
   ampli(i) = vt(i)
   phase(i) = st(i)
   vt(i) = 0.0d0
   st(i) = 0.0d0
end do
!
!     shift the phase angles into the standard range
!
do i = 1, 6
   do while (phase(i) .lt. -180.0d0)
      phase(i) = phase(i) + 360.0d0
   end do
   do while (phase(i) .gt. 180.0d0)
      phase(i) = phase(i) - 360.0d0
   end do
end do
!
!     convert input torsional parameters to storage format
!
do i = 1, 6
   k = ft(i)
   if (k .eq. 0) then
      exit
   else if (k .le. 6) then
      vt(k) = ampli(i)
      st(k) = phase(i)
   end if
end do
return
end subroutine torphase
