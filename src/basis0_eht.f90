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
!     basis0_eht: first step to setup EHT basis: determine limits
!     of single functions
!
!     part of QMDFF
!
subroutine basis0_eht(n,at,nel,nbf)
use qmdff
implicit none
integer::n,nbf,nao,nel
integer::at(n)
integer::i
real(kind=8)::z1(n) !is now in qmdff.f90! 
nbf=0
nel=0
z1=0

!     
!      same procedure as in subroutine valel (ff_tool.f90)
!      (at least for the z(n) array...)
!
do i=1,n
!     H: 1s basis (?)
   if (at(i).le.2)then
      nbf=nbf+1
      z1(i)=at(i)
   end if
   if (at(i).le.10.and.at(i).gt.2) then
!     Li-F: 2s2p
      nbf=nbf+4
      z1(i)=at(i)-2
   end if
   if (at(i).le.13.and.at(i).gt.10) then
!     Na-Al: 3s3p
      nbf=nbf+4
      z1(i)=at(i)-10
   end if
   if (at(i).le.18.and.at(i).gt.13) then
!     Si-Ar: 3s3p
      nbf=nbf+4
      z1(i)=at(i)-10
   end if
   if ((at(i).le.36.and.at(i).ge.30).or. &
     &  at(i).eq.19.or. at(i).eq.20) then
!     K,Ca,Zn-Kr: 4s4p
      nbf=nbf+4
      z1(i)=at(i)-18
      if (at(i).ge.30)z(i)=at(i)-28
   end if
   if (at(i).le.29.and.at(i).ge.21) then
!     Sc-Cu: 4s4p3d
      nbf=nbf+10
      z1(i)=at(i)-18
   end if
   if (at(i).le.47.and.at(i).ge.39) then
!     Y-Ag: 5s5p4d 
      nbf=nbf+10
      z1(i)=at(i)-36
   end if
   if (at(i).gt.47.and.at(i).le.54) then
!     In-Xe: 5s5p
      nbf=nbf+4
      z1(i)=at(i)-46
   end if

end do
!
!     (old command idint(..) converted into int(..)) 
!     MODDED 15.03.17: 
!     manually calculate the number of electrons in the system!!
!
nel=0
do i=1,n
   nel=nel+z1(i)
end do
!nel = int(sum(z))

if(nbf.gt.maxao) stop 'TB basis too large'
end subroutine basis0_eht
