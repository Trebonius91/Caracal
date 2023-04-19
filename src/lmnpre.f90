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
!     subroutine lmnpre: Precalculation of EHT overlap integral terms 
!         (depends only on lm)
!
!     part of QMDFF
!
subroutine lmnpre(l,m,n,lmnexp,lmnfak)
implicit real*8(a-h,o-z)
real*8  lmnfak
integer lmnexp
dimension dftr(7)
data dftr /1.d0,1.d0,3.d0,15.d0,105.d0,945.d0,10395.d0/
lmnfak=0
lmnexp=0
lh=l/2
!
!     Fortran 66: Arithmetic IF:
!     if(2*lh-l) 50,1,50
!     (result negative: goto 50, zero: goto 1, positive: goto 50)
!                                                                    
if (2*lh-l .eq. 0d0) then
   mh=m/2
else
   return
end if
if (2*mh-m .eq. 0d0) then
   nh=n/2
else
   return
end if
if(2*nh-n .eq. 0d0) then
   lmnexp=lh+mh+nh
else
   return
end if
lmnfak=dftr(lh+1)*dftr(mh+1)*dftr(nh+1)

return
end subroutine lmnpre

