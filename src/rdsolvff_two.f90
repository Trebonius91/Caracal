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
!     subroutine rdsolvff_two: read in informations for QMDFF2 
!
!     part of QMDFF
!
subroutine rdsolvff_two(n,xyz,at,q,imass,d,scalehb,scalexb,fname)
use qmdff
implicit none  
integer::n,at(n),imass(n)
real(kind=8)::xyz(3,n),q(n),d
real(kind=8)::scalehb(94),scalexb(94)
character(len=*)::fname

integer::i,j,nn,nterm,idum
real(kind=8)::qdum,r,xx(10),dum1,dum2
character(len=80)::atmp

open(unit=1,file=fname)

read (1,'(a)') atmp   
call readl(atmp,xx,nn)
n=idint(xx(1))
qdum=xx(2)
d=xx(3)

if (nn.eq.3.and.d.lt.0) then
   d=xx(3)
   write(*,*) 'taking solvent density from file ',d
end if
read (1,'(a)') atmp   

do i=1,n
   read (1,*)at(i),xyz(1:3,i),q(i),imass(i)
end do
read (1,*) nbond_two,nangl_two,ntors_two,nhb_two,nnci_two
if (nbond_two.gt.ndim.or.nangl_two.gt.ndim.or.ntors_two.gt.ndim) &
    &  stop 'too many FF terms'
do i=1,nbond_two
   read (1,*)bond_two(1:2,i),vbond_two(1:3,i)
end do
do i=1,nangl_two
   read (1,*)angl_two(1,i),angl_two(2,i),angl_two(3,i), &
     &     vangl_two(1:2,i)
enddo
do i=1,ntors_two
   read (1,*)tors_two(1:6,i),vtors_two(1:2,i), &
       &   (vtors_two(3*(j-1)+3,i), &
       &    vtors_two(3*(j-1)+4,i), &
       &    vtors_two(3*(j-1)+5,i),j=1,tors_two(5,i))
end do
do i=1,ntors_two
   do j=1,tors_two(5,i)
      vtors_two(3*(j-1)+4,i)=vtors_two(3*(j-1)+4,i)*pi
   end do
end do
if(nhb_two.gt.0) read(1,'(6(3i4,2x))')hb_two(1:3,1:nhb_two)
do i=1,nhb_two
   if (at(hb_two(3,i)).eq.1) then
!
!     subroutine hbpara is totally FF-unspecific!
!
      call hbpara(10.0d0,5.0d0,q(hb_two(1,i)),dum1)
      vhb_two(1,i)=dum1*scalehb(at(hb_two(1,i)))
      call hbpara(10.0d0,5.0d0,q(hb_two(2,i)),dum1)
      vhb_two(2,i)=dum1*scalehb(at(hb_two(2,i)))
   else
      call hbpara(-6.5d0,1.0d0,q(hb_two(3,i)),dum1)
      vhb_two(1,i)=scalexb(at(hb_two(3,i)))*dum1
   end if
end do
read(1,'(8(3i4,2x))')nci_two(1:3,1:nnci_two)
close(1)

q=q*qdum 
 
return
end subroutine rdsolvff_two
