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
!     subroutine bonde: Calculate energies for all single bonds 
!           in the system
!
!     part of QMDFF
!
subroutine bonde(n,at,xyz,pair,nb,cring,ringsize)
use qmdff
implicit none
integer::n,at(n)
real(kind=8)::xyz(3,n)
integer::pair(n*(n+1)/2),nb(20,n)
integer::cring(8,n),ringsize(n)
character(len=2)::asym
character(len=16)::tag(n*(n+1)/2)
integer::i,j,k,m,ni,nj,i1,i2,n1,n2,kk
real(kind=8)::sumb,e,be(n*(n+1)/2)
logical::samering

be=0
e=0
k=0
kk=0
do i=1,n
   do j=1,i
      k=k+1
      ni=nb(20,i)
      nj=nb(20,j)
      sumb=0
      if (pair(k).eq.1) then

         kk=kk+1
!
!     1,2 terms  
!          
         do m=1,nbond
            if (bond(1,m).eq.i.and.bond(2,m).eq.j) then
               sumb=sumb+vbond(2,m)
            end if
            if(bond(2,m).eq.i.and.bond(1,m).eq.j) then
               sumb=sumb+vbond(2,m)
            end if
         end do
!
!     1,3 terms  
!          
         do n1=1,ni
            i1=nb(n1,i)
            if (i1.eq.i) cycle
            do m=1,nbond
               if (bond(1,m).eq.i1.and.bond(2,m).eq.j) then
                  sumb=sumb+vbond(2,m)*0.5
               end if
               if (bond(2,m).eq.i1.and.bond(1,m).eq.j) then
                  sumb=sumb+vbond(2,m)*0.5
               end if
            end do
         end do
         do n2=1,nj
            i2=nb(n2,j)
            if (i2.eq.j) cycle
            do m=1,nbond
               if (bond(1,m).eq.i2.and.bond(2,m).eq.i) then
               sumb=sumb+vbond(2,m)*0.5
               end if
               if (bond(2,m).eq.i2.and.bond(1,m).eq.i) then
               sumb=sumb+vbond(2,m)*0.5
               end if
            end do
         end do
         be(kk)=sumb*627.51
         write(tag(kk),'(a2,i3,2x,a2,i3)') &
             &  asym(at(i)),i,asym(at(j)),j
         e=e+sumb
      end if
   end do
end do
write(10,*)
write(10,*) '======================'
write(10,*) 'covalent bond energies'
write(10,*) '======================'
write(10,*) 
write(10,*) "--------------------------------"
write(10,*) "El1 #1 El2 #2      E(kcal)"
write(10,*) "--------------------------------"
do k=1,kk
    write(10,'(A,A,F9.2)') " ",tag(k),be(k)
!   write(10,'(4(A14,''='',F7.1,4x))')(trim(tag(k)),be(k),k=1,kk)
end do
write(10,*) "--------------------------------"
write(10,'('' sum BE (kcal):'',F12.2)')e*627.51
write(10,*) "--------------------------------"

return
end subroutine bonde
