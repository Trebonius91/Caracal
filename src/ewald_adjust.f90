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
!     subroutine adjust: modifies site bounds on the PME grid and returns
!     an offset into the B-spline coefficient arrays
!
subroutine ewald_adjust (offset,nfft,nchk,amin,amax,cmin,cmax)
implicit none
integer::offset
integer::nfft,nchk
integer::amin,amax
integer::cmin,cmax
!
!     modify grid offset and bounds for site at edge of chunk
!
offset = 0
if (nchk .ne. 1) then
   if (amin.lt.cmin .or. amax.gt.cmax) then
      if (amin.lt.1 .or. amax.gt.nfft) then
         if (cmin .eq. 1) then
            offset = 1 - amin
            amin = 1
         else if (cmax .eq. nfft) then
            amax = nfft
            amin = amin + nfft
         end if
      else
         if (cmin .gt. amin) then
            offset = cmin - amin
            amin = cmin
         else
            amax = cmax
         end if
      end if
   end if
end if
offset = offset + 1 - amin
return
 
end subroutine ewald_adjust

