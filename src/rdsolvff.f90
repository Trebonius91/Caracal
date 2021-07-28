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
!     subroutine rdsolvff: read in information for QMDFF1
!
!     part of QMDFF
!
subroutine rdsolvff(n,xyz,at,q,imass,d,scalehb,scalexb,fname)
use qmdff
implicit none  
integer::n,at(n),imass(n)
real(kind=8)::xyz(3,n),q(n),d
real(kind=8)::scalehb(94),scalexb(94)
character(len=*)::fname

integer::i,j,nn,nterm,idum
real(kind=8)::qdum,r,xx(10),dum1,dum2
character(len=160)::atmp

integer,allocatable::linenum(:)
integer::numnci,remain

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
!write(10,*)'                           intramolecular solvent FF:'
!write(10,*)'method used for charges          ',trim(atmp)
!
!     read in coordinates, charges and molecule numbers
!     Changed from original; there was imass(i) instead of molnum(i)
!
do i=1,n
   read (1,*)at(i),xyz(1:3,i),q(i),molnum(i)
end do
!
!     Determine the number of molecules in the QMDFF
!
nmols=maxval(molnum)
!
!     For box/solvent QMDFF: allocate and fill the global charge array!
!
if (nmols .gt. 1) then
   allocate(q_glob(n))
   q_glob(:)=q(1:n)
end if
!
!     read in number of FF terms, check if too many
!
read (1,*) nbond,nangl,ntors,nhb,nnci
if (nbond.gt.ndim.or.nangl.gt.ndim.or.ntors.gt.ndim) &
    &  stop 'too many FF terms'
!
!     read in bond terms
!
do i=1,nbond
   read (1,*)bond(1:2,i),vbond(1:3,i)
enddo
!
!     read in bend angle terms
!
do i=1,nangl
 read (1,*)angl(1,i),angl(2,i),angl(3,i), &
     &    vangl(1:2,i)
enddo
!
!     read in torsion terms
!
do i=1,ntors
 read (1,*)tors(1:6,i),vtors(1:2,i), &
      &  (vtors(3*(j-1)+3,i), &
      &   vtors(3*(j-1)+4,i), &
      &   vtors(3*(j-1)+5,i),j=1,tors(5,i))
enddo
do i=1,ntors
   do j=1,tors(5,i)
      vtors(3*(j-1)+4,i)=vtors(3*(j-1)+4,i)*pi
   enddo
enddo
!
!     read in hydrogen/halogen bond terms
!
!if (nhb.gt.0) read(1,'(6(3i4,2x))')hb(1:3,1:nhb)
if (nhb .gt. 0) then
   numnci=0
   allocate(linenum(18))
   do i=1,int(nhb/6)
      read(1,*) linenum
      do j=1,6
         numnci=numnci+1
         hb(1:3,numnci)=linenum((j-1)*3+1:(j-1)*3+3)
      end do
   end do
   deallocate(linenum)
   remain=(nhb-numnci)*3
   if (remain .gt. 0) then
      allocate(linenum(remain))
      read(1,*) linenum
      do j=1,remain
         numnci=numnci+1
         hb(1:3,numnci)=linenum((j-1)*3+1:(j-1)*3+3)
      end do
      deallocate(linenum)
   end if
end if

do i=1,nhb
 if (at(hb(3,i)).eq.1) then
   call hbpara(10.0d0,5.0d0,q(hb(1,i)),dum1)
   vhb(1,i)=dum1*scalehb(at(hb(1,i)))
   call hbpara(10.0d0,5.0d0,q(hb(2,i)),dum1)
   vhb(2,i)=dum1*scalehb(at(hb(2,i)))
 else
   call hbpara(-6.5d0,1.0d0,q(hb(3,i)),dum1)
   vhb(1,i)=scalexb(at(hb(3,i)))*dum1
 endif
enddo
!
!     read in noncovalent interactions
!
if (nnci .gt. 0) then
   numnci=0
   allocate(linenum(24))
   do i=1,int(nnci/8)
      read(1,*) linenum
      do j=1,8
         numnci=numnci+1
         nci(1:3,numnci)=linenum((j-1)*3+1:(j-1)*3+3)
      end do
   end do
   deallocate(linenum)
     remain=(nnci-numnci)*3
   if (remain .gt. 0) then
      allocate(linenum(remain))
      read(1,*) linenum
      do j=1,remain
         numnci=numnci+1
         nci(1:3,numnci)=linenum((j-1)*3+1:(j-1)*3+3)
      end do
      deallocate(linenum)
   end if
end if 


do i=1,nnci
   if (nci(1,i) .eq. 0 .or. nci(2,i) .eq. 0 .or. nci(3,i) .eq. 0) then
      write(*,*) "ERROR! The noncovalent interactions were not read in correctly!"
      write(*,*) "Please check the formate of the QMDFF!"
      call fatal
   end if
end do

close(1)

q=q*1.d0!qdum 
!
!     Write global arrays for box/solvent QMDFFs
!
scalexb_glob(:)=scalexb(:)
scalehb_glob(:)=scalehb(:)

!write(10,*)'nbond,nangl,ntors   :',nbond,nangl,ntors
!write(10,*)'NCI terms added     :',nnci                     
!write(10,*)'HB/XB  terms added  :',nhb

return
end subroutine rdsolvff
