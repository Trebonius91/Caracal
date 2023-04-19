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
!     subroutine getrot: determine the rotational symmetry number for bond ij
!     --> use hybridizations and ringsizes for that manner
!
!     part of QMDFF
!
subroutine getrot(n,at,i,j,hyb,cring,ringsize,xn)
implicit none
integer::n,at(n),hyb(n),i,j
integer::cring(8,n),ringsize(n)
real(kind=8)::xn
logical::samering
integer::metal(94)
!
!     default barrier at 90    
!     --> if no other criteria can be fulfilled...
!   
xn = 2
!
!     hybridizations of atoms::
!     double bond
!
if (hyb(i).eq.2.and.hyb(j).eq.2) return
!
!     sp2-sp3   
!   
if (hyb(i).eq.2.and.hyb(j).eq.3) then
   xn = 2
   return
end if
if (hyb(j).eq.2.and.hyb(i).eq.3) then
   xn = 2
   return
end if
!
!     sp3-sp3, default barrier at 60,120,180 ...
!
if (hyb(i).eq.3.and.hyb(j).eq.3) then
    xn = 3
!
!     small 4,5-rings have various values 
!            
    if (samering(n,i,j,cring,ringsize)) then
       if (ringsize(i).eq.3)  xn = 2
       if (ringsize(i).eq.4)  xn = 0
!
!     this choice makes cyclopentane right
!
       if(ringsize(i).eq.5)   xn = 0
       if(ringsize(i).eq.6)   xn = 0
!
!     cycloheptane is a bit too floopy but better than for xn=3   
!      
       if(ringsize(i).gt.6)   xn = 0
   end if
   return
end if
!
!     metals have more often 4-fold axis!
!
if(metal(at(i)).eq.1.or.metal(at(j)).eq.1) xn = 4

return
end subroutine getrot
