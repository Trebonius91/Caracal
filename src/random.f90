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
!     function random: generates a random number on [0,1] via a long
!     period generator due to L'Ecuyer with Bays-Durham shuffle!
!     literature references:
!     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)!
!     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
!     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
!     University Press, 1992, Section 7.1
!
!     part of EVB
!
function random ()
use general
implicit none
integer::im1,ia1,iq1,ir1
integer::im2,ia2,iq2,ir2
integer::big,nshuffle
integer::imm1,ndiv
real(kind=8)::factor
parameter (im1=2147483563)
parameter (ia1=40014)
parameter (iq1=53668)
parameter (ir1=12211)
parameter (im2=2147483399)
parameter (ia2=40692)
parameter (iq2=52774)
parameter (ir2=3791)
parameter (big=141803398)
parameter (nshuffle=32)
parameter (imm1=im1-1)
parameter (ndiv=1+imm1/nshuffle)
parameter (factor=1.0d0/im1)
integer::i,k,iy,next
integer::seed,seed2
integer::year,month,day
integer::hour,minute,second
integer::ishuffle(nshuffle)
real(kind=8)::random
logical::first
character(len=20)::keyword
character(len=120)::record
character(len=120)::string
save first
save seed,seed2
save iy,ishuffle

first=.true.
!
!     random number seed is first set to a big number,
!     then incremented by the seconds elapsed this decade
!
if (first) then
   first = .false.
   seed = big
   call calendar (year,month,day,hour,minute,second)
   year = mod(year,10)
   seed = seed + 32140800*year + 2678400*(month-1)
   seed = seed + 86400*(day-1) + 3600*hour
   seed = seed + 60*minute + second
!
!     search the keywords for a random number seed
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      if (keyword(1:11) .eq. 'RANDOMSEED ') then
         string = record(next:120)
         read (string,*,err=10)  seed
         seed = max(1,seed)
      end if
      10  continue
   end do
!
!     warm up and then load the shuffling table
!
   seed2 = seed
   do i = nshuffle+8, 1, -1
      k = seed / iq1
      seed = ia1 * (seed-k*iq1) - k*ir1
      if (seed .lt. 0)  seed = seed + im1
      if (i .le. nshuffle)  ishuffle(i) = seed
   end do
   iy = ishuffle(1)
end if
!
!     get a new random number value each call
!
k = seed / iq1
seed = ia1*(seed-k*iq1) - k*ir1
if (seed .lt. 0)  seed = seed + im1
k = seed2 / iq2
seed2 = ia2*(seed2-k*iq2) - k*ir2
if (seed2 .lt. 0)  seed2 = seed2 + im2
i = 1 + iy/ndiv
iy = ishuffle(i) - seed2
ishuffle(i) = seed
if (iy .lt. 1)  iy = iy + imm1
random = factor * iy

return
end function random
