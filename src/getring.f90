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
!     subroutine getring: pick an atom of QMDFF reference and check its 
!        memberance to a ring using bond orders
!
!     part of QMDFF
!

subroutine getring(n,nbin,a0,c,iring)
implicit none
integer::n,nbin(20,n),a0,iring,i,nb(20,n)
integer::i1,i2,i3,i4,i5,i6,i7,i8
integer::n0,n1,n2,n3,n4,n5,n6,n7,n8
integer::a1,a2,a3,a4,a5,a6,a7,a8
integer::c(8)
logical::chk

nb=nbin
do i=1,n
   if(nb(20,i).eq.1)nb(20,i)=0
enddo

iring =0
c(1:8)=0

n0=nb(20,a0)
!
!    great loop: check for ringsizes up to seven atoms!
!
do i1=1,n0
   a1=nb(i1,a0)
   if (a1.eq.a0) cycle
   n1=nb(20,a1)
   do i2=1,n1
      a2=nb(i2,a1)
      if (a2.eq.a1) cycle
      n2=nb(20,a2)
      do i3=1,n2
         a3=nb(i3,a2)
         n3=nb(20,a3)
         if (a3.eq.a2) cycle
         c(1)=a1
         c(2)=a2
         c(3)=a3
         if (a3.eq.a0.and.chk(n,3,c)) then
            iring=3
            goto 99
         end if
         do i4=1,n3
            a4=nb(i4,a3)
            n4=nb(20,a4)
            if (a4.eq.a3) cycle
            c(4)=a4
            if (a4.eq.a0.and.chk(n,4,c)) then
               iring=4
               goto 99
            end if
            do i5=1,n4
               a5=nb(i5,a4)
               n5=nb(20,a5)
               if (a5.eq.a4) cycle
               c(5)=a5
               if (a5.eq.a0.and.chk(n,5,c)) then
                  iring=5
                  goto 99
               end if
               do i6=1,n5
                  a6=nb(i6,a5)
                  n6=nb(20,a6)
                  if (a6.eq.a5) cycle
                  c(6)=a6
                  if (a6.eq.a0.and.chk(n,6,c)) then
                     iring=6
                     goto 99
                  end if
                  do i7=1,n6
                     a7=nb(i7,a6)
                     n7=nb(20,a7)
                     if (a7.eq.a6) cycle
                     c(7)=a7
                     if (a7.eq.a0.and.chk(n,7,c)) then
                        iring=7
                        goto 99
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do
end do

return

99 continue
call isort(iring,c)

return
end subroutine getring

