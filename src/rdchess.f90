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
!     subroutine rdchess: read cp2k QMDFF reference hessian
!
!     part of QMDFF
!
subroutine rdchess(nat,nat3,h,fname)
use general
use qmdff
implicit none
integer::nat,nat3
real(kind=8)::h(nat3,nat3)
real(kind=8)::mass3(nat3)
! array with all gradients for numerical hessian calculation:
real(kind=8)::grads(3*nat3,nat3) 
real(kind=8)::diffstep ! step for numerical differentiation
!integer::iat(nat)
!real(kind=8)::bohr,hartree,joule
character(len=*)::fname
integer::iunit,i,j,mincol,maxcol,idum,m,k
character(len=5)::adum
character(len=80)::a80
character(len=2)::asym
logical::snfcp2k
!parameter (bohr=0.52917721092d0)
!parameter (hartree=627.5094743d0)
!parameter (joule=4.1840d0)


snfcp2k=.true.
!write(10,*)
!write(10,*) 'reading <',trim(fname),'>'
iunit=12
open(unit=iunit,file=fname)
!
!     define number of 6 column lines
!
m=nat3/5
if (mod(nat3,5).gt.0) m=m+1

!
!     read in the hessian from file
!

!100 continue  
!read(iunit,'(a)',end=300)a80
!if (index(a80,'$act_energy').ne.0) snforca=.false.
!if (index(a80,' VIB| Hessian in cartesian coordinates').ne.0) then
!   read(iunit,'(a)',end=300) a80
!   maxcol = 0
!   do k=1,m
!      read(iunit,'(a)',end=300)a80  ! is the end statement really needed?
!      read(iunit,'(a)')a80
!      mincol = maxcol + 1
!      maxcol = min(maxcol+5,nat3)
!      do i=1,nat3
!         read(iunit,*) idum,adum,(h(i,j),j=mincol,maxcol)
!      end do
!   end do
!end if
!goto 100
!
!300 continue 
!
!     TEST1: read in gradients and calculate hessian from them
!
100 continue
read(iunit,'(a)',end=300)a80
!if (index(a80,'$act_energy').ne.0) snforca=.false.
if (index(a80,' VIB| Vibrational Analysis Info').ne.0) then
!   read(iunit,'(a)',end=300) a80
   do i=1,2*nat3
!
!    read the header lines of each block
!
      read(iunit,'(a)') a80
      read(iunit,'(a)') a80
      read(iunit,'(a)') a80
!
!    fill the big gradient array with the actual gradient
!
      do k=1,nat
         read(iunit,*) adum,adum,grads(i,(k-1)*3+1),grads(i,(k-1)*3+2),grads(i,(k-1)*3+3)
      end do
   end do
   goto 300
end if

goto 100

300 continue
diffstep=1.0d-2
do i=1,nat3
   do j=1,nat3
      h(i,j)=-(grads(i,j)-grads(i+nat3,j))/(2d0*diffstep)
   end do
end do

!
!     TEST: convert hessian entries from Ht/Ang^2 to Ht/bohr^2 (12.05.2017)
!
!h=h/(hartree)
!h=h*bohr*bohr/(hartree*joule)
!
!    determine masses of all atoms in the system via the element symbols
!
!do i=1,nat
 !  call atommass(i)
!   write(*,*) i,iat(i),asym(iat(i))
!   name(i)=asym(iat(i))
!   call atommass(i)
!end do
!
!   Define massses for each of the 3N cartesian coordinates
!
!do i=1,nat
!   do j=1,3
!      mass3((i-1)*3+j)=mass(i)
!   enddo
!enddo

!
!     multiplie with the masses to extinct the mass weighted hessian
!
!do i=1,nat3   
!   do j=1,nat3
!      h(i,j)=h(i,j)*sqrt(mass3(i)*mass3(j))  
!      if (i .ne. j) then
!         h(i,j)=-h(i,j)
!      end if
!   end do
!end do
 
close(iunit,status='keep')
!
!     scale it by global defined value
!

h=h*scalh

return
end subroutine rdchess
