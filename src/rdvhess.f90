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
!     subroutine rdvhess: read VASP QMDFF reference hessian
!
!     part of QMDFF
!
subroutine rdvhess(nat3,h,fname,iat)
use qmdff
implicit none
integer::nat,nat3
real(kind=8)::h(nat3,nat3)
real(kind=8),allocatable::row_act(:)
real(kind=8)::h_act
real(kind=8)::atmass(nat3/3)
character(len=*)::fname
integer::iat(nat3/3)
integer::iunit,i,j,mincol,maxcol,idum,m,k,l
integer::readstat
integer::nat_active
integer::ind_trans(nat3/3)
character(len=5)::adum
character(len=80)::a80
character(len=5000)::along
character(len=5)::head_inds(nat3)
logical::snforca

character(len=1)::JOBZ,UPLO
integer::Nn,LDA,LDU,LDV,LWORK,INFO
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::A_mat


nat=nat3/3
snforca=.true.
!write(10,*)
!write(10,*) 'reading <',trim(fname),'>'
iunit=12
open(unit=iunit,file=fname)
!
!     Predefine hessian as diagonal matrix with low diagonal elements 
!     for the case that come atoms are not activated for dynamics
!
h=0.d0
do i=1,nat3
   h(i,i)=0.001d0
end do

!
!     read in the hessian from file
!

head_inds="XXXXX"
ind_trans=0
100 continue  
read(iunit,'(a)',end=300)a80
if (index(a80,'SECOND DERIVATIVES (NOT SYMMETRIZED) ').ne.0) then
   write(*,*) "Gguopgup"
!
!     First read the comment line with the symbols to extract the number of atoms 
!     activated for dynamics
!
   read(iunit,*)
   read(iunit,'(a)',iostat=readstat) along 
   read(along,*,iostat=readstat) head_inds 
   k=0
   do i=1,nat
      do j=1,nat3
         adum=head_inds(j)
         if (adum .eq. "XXXXX") exit
         do m=1,5
            if (adum(m:m) .eq. "X") adum(m:m)=" "
            if (adum(m:m) .eq. "Y") adum(m:m)=" "
            if (adum(m:m) .eq. "Z") adum(m:m)=" "
         end do
         read(adum(1:4),*) idum
         if (idum .eq. i) then
            k=k+1
            ind_trans(k) = idum
            exit
         end if
      end do
   end do
   nat_active=k
   allocate(row_act(3*nat_active))
!
!     Now read in the Hessian for the activated atoms
!
   do i=1,nat_active
      do j=1,3
         read(iunit,*) adum, row_act
         do k=1,nat_active
            do l=1,3
               h((ind_trans(i)-1)*3+j,(ind_trans(k)-1)*3+l)=row_act((k-1)*3+l)
            end do
         end do
      end do
   end do
end if
goto 100

300 continue
close(iunit,status='keep')
do i=1,nat3
   write(28,*)i,h(i,:)
end do
!
!     convert it to atomic units, with advice from:
!     https://www.vasp.at/forum/viewtopic.php?t=19024
!
!     Step 1: symmetrization
!
do i=1,nat3
   do j=i+1,nat3
      h_act=0.5d0*(h(i,j)+h(j,i))
      h(i,j)=h_act
      h(j,i)=h_act 
   end do
end do

!
!     Step 2: Calculate mass-normalized force constant (standard Hessian)
!
!     First, fill array with atomic masses 
!
do i=1,nat
   call atomname(iat(i),adum)
   call atommass_loc(adum,atmass(i))
end do
!
!     Then, divide the hessian through the respective mass root products
!
do i=1,nat
   do j=1,nat
      do k=1,3
         do l=1,3
            h((i-1)*3+k,(j-1)*3+l)=h((i-1)*3+k,(j-1)*3+l)/sqrt(atmass(i)*atmass(j))
         end do
      end do
   end do
end do
!
!     Step 3: Convert the Hessian to atomic units!
!
h=-h!*1000


!
!   TEST    Diagonalize the hessian
!
!   therafter: W=frequencies, A=normal coordinates
!
   JOBZ='V' !eigenvalues and eigenvectors(U)
   UPLO='U' !upper triangle of a
   Nn=3*nat
   LDA=nat3
   INFO=0
   LWORK=nat3*nat3-1
   allocate(A_mat(nat3,nat3))
   allocate(W(Nn))
   allocate(WORK(LWORK))
   A_mat=h
   call DSYEV(JOBZ,UPLO,nat3,A_mat,LDA,W,WORK,LWORK,INFO)
   if (info .ne. 0) then
      write(*,*) "The diagonalization of the hessian matrix with Lapack returned"
      write(*,*) " an error code!"
      stop
   end if
!
!     Obtain frequencies with correct unit etc.
!
   do i = 1,nat3
      W(i)=sign(sqrt(abs(W(i))),W(i))
   end do
!
!     Do conversion according to Mopac manual: 
!     http://openmopac.net/manual/Hessian_Matrix.html
!
!     1: from eV/Ang^2 to kcal/(mol*Ang^2)
!
   W=W*sqrt(23.0609d0)
!
!     2: from kcal/(mol*Ang^2) to millidynes/Ang (and to Newton/m/kg)
!
   W=W*sqrt(1E8*4184d0/(1E-10*6.023E23))
!
!     3: from millidynes/Ang to cm^-1
!
   W=W*1302.79d0

   do i=1,nat3
      write(*,*) i,W(i)
   end do

h=h*0.0367493/0.52917721092d0**2

return
end subroutine rdvhess
