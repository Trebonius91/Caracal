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
!     setvalel: Defines number of valence electrons for EHT calculations
!
!     part of QMDFF
!
subroutine setvalel
use qmdff
implicit none
integer::i,at

do i=1,54
   at=float(i)
!     H
   z(i)=at
   if(i.le.10.and.i.gt.2)then
!     Li-F
      z(i)=at-2
   endif
   if(i.le.13.and.i.gt.10)then
!     Na-Al
      z(i)=at-10
   endif
   if(i.le.18.and.i.gt.13)then
!     Si-Ar
      z(i)=at-10
   endif
   if ((i.le.36.and.i.ge.30).or.i.eq.19.or. i.eq.20) then
!     K,Ca,Zn-Kr 
!     4s4p shells (?)
      z(i)=at-18
      if(i.ge.30)z(i)=at-28
   endif
   if(i.le.29.and.i.ge.21) then
!     Sc-Cu
      z(i)=at-18
   endif
   if(i.le.47.and.i.ge.37) then
!     Y-Ag  
      z(i)=at-36
   endif
   if(i.gt.47.and.i.le.54) then
!     In-Xe 
      z(i)=at-46
   endif
enddo
!
!     not well defined for these metals
!     some side group elements etc.
!
do i=57,80
   z(i)=4.0
enddo

do i=81,86
   z(i)=float(i)-78.0
enddo

z(87)=1
z(88)=2

z(89:94)=4

end subroutine setvalel
