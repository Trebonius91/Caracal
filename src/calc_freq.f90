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
!     subroutine calc_freq: performs a normal coordinate analysis on a given 
!     hessian matrix
!
!     part of EVB
!

subroutine calc_freq(hess,freqs,eigvecs,printfreq)
use evb_mod
use general
use pbc_mod

implicit none
integer::i,j,k,l,inc,block,first
integer::dg_evb_mode,mat_size,nat,nat3
real(kind=8)::hess(3*natoms,3*natoms)
real(kind=8)::hess_backup(3*natoms,3*natoms)
real(kind=8)::freqs(3*natoms),mass3(3*natoms)
real(kind=8)::eigvecs(3*natoms,3*natoms)
! for matrix diagonalization
character(len=1)::JOBZ,UPLO
character(len=:),allocatable::filename1
character(len=50)::syscall
integer::M,Nn,LDA,LDU,LDV,LWORK,INFO
real(kind=8),dimension(:),allocatable::WORK,W
real(kind=8),dimension(:,:),allocatable::A_mat,coord
real(kind=8),dimension(:,:,:),allocatable::n_vib
real(kind=8),dimension(:,:),allocatable::U
real(kind=8)::com_coord(3)
integer::NRHS,LDB
integer::ind
logical::printfreq ! if normal mode frequencies shall be printed to log
real(kind=8)::pi
real(kind=8)::factor
!   for normal mode intensities
real(kind=8),allocatable::n_vib1(:,:,:),n_vib2(:,:,:)
real(kind=8)::deltaxyz
real(kind=8),allocatable::dipgrad(:,:,:),dxyz(:,:),ddip(:,:,:)
real(kind=8),allocatable::intens(:,:),totintens(:)
real(kind=8),allocatable::coord_mass(:)
real(kind=8)::cnvint
real(kind=8),allocatable::dip_list2(:,:)
!
!     unit conversion factors
!
real(kind=8)::amu2au,au2cm
parameter (amu2au=1.66053886E-27/9.10938215E-31)
parameter (au2cm =219474.63067d0)
nats=natoms
nat=nats
nat3=3*nats
pi=3.1415926535897932384d0
!
!   Define massses for each of the 3N cartesian coordinates
!
do i=1,nats
   do j=1,3
      mass3((i-1)*3+j)=1.d0/sqrt(mass(indi(i)))
   enddo
enddo

!
!    Store hessian in backup to return it unchanged
!
hess_backup=hess
!
!    For pGFN-FF and QMDFF: multiply hessian twice by bohr factor
!
if (.not. gfn_xtb) then
   hess=hess/bohr**2
end if
!
!   Include mass weighted coordinates!
!
do i=1,nats*3
   do j=1,nats*3
       hess(i,j)=hess(i,j)*(mass3(i)*mass3(j))
   enddo
enddo
!
!   diagonalize the matrix
!   therafter: W=frequencies, A=normal coordinates
!
JOBZ='V' !eigenvalues and eigenvectors(U)
UPLO='U' !upper triangle of a
Nn=3*(nats-fix_num)
LDA=Nn
INFO=0
LWORK=Nn*Nn-1

allocate(A_mat(Nn,Nn))
allocate(W(Nn))
allocate(WORK(LWORK))
allocate(n_vib(Nn,nats,3))
allocate(coord(3,natoms))
A_mat=hess
call DSYEV(JOBZ,UPLO,Nn,A_mat,LDA,W,WORK,LWORK,INFO)
if (info .ne. 0) then
   write(*,*) "The diagonalization of the hessian matrix in calc_freq.f90 returned"
   write(*,*) " an error code!"
   call fatal
end if
!
!   The freqequencies can be expressed as: 1/(2*pi)*sqrt(f)
!   where f is the eigenvalue
!   convert frequencies into cm^-1
!
!   Presumably, an arbitrary factor of 1.0823 is needed to produce 
!   good frequencies! 
!   ---> corrected 07.11.2017!! (square root was applied too early)
!
!do i = 1, 3*nats
!   W(i)=au2cm*sign(sqrt(abs(W(i))),W(i))/sqrt(amu2au)
!end do
!if (printfreq) then
!   write(15,*)"*The normal mode frequencies (cm-1) of the structure are:"
!   write(15,*)
!   do i=1,3*nats
!      write(15,*) i, W(i)
!   end do
!   write(15,*)
!end if
!
!    for return of frequencies
!
!freqs=W
!
!   Now, if desired, calculate the normal mode cartesian displacements!!
!
!eigvecs=A_Mat
!
!     Obtain frequencies with correct unit etc.
!
do i = 1, 3*(natoms-fix_num)
   W(i)=sign(sqrt(abs(W(i))),W(i))
end do
!
!     Do conversion according to Mopac manual: 
!     http://openmopac.net/manual/Hessian_Matrix.html
!
!     1: from Hartree/Ang^2 to kcal/(mol*Ang^2)
!
W=W*sqrt(627.503)
!
!     2: from kcal/(mol*Ang^2) to millidynes/Ang (and to Newton/m/kg)
!
W=W*sqrt(1E8*4184d0/(1E-10*6.023E23))
!
!     3: from millidynes/Ang to cm^-1
!
W=W*1302.79d0
!
!     Calculate normal mode vibration vectors (for all atoms and 
!      only for active atoms)
!
allocate(coord_mass(3*(natoms-fix_num)))
allocate(n_vib1(3,natoms-fix_num,Nn))
allocate(n_vib2(3,natoms,Nn))
k=0
do i=1,natoms
   if (at_move(i)) then
      k=k+1
      coord_mass((k-1)*3+1)=mass(i)
      coord_mass((k-1)*3+2)=mass(i)
      coord_mass((k-1)*3+3)=mass(i)
   end if
end do
do i=1,3*(natoms-fix_num)
   do j=1,natoms-fix_num
      do k=1,3
         n_vib1(k,j,i)=A_mat((j-1)*3+k,3*(natoms-fix_num)-i+1)/sqrt(coord_mass((j-1)*3+k))
      end do
   end do
end do
!
!     Fill full (all n atoms) normal mode vectors
!
do i=1,3*(natoms-fix_num)
   ind=1
   do j=1,natoms
      if (at_move(j)) then
         n_vib2(:,j,i)=n_vib1(:,ind,i)
         ind=ind+1
      else
         n_vib2(:,j,i)=0.d0
      end if
   end do
end do

!
!    Calculate the normal mode intensities, if available   
!    adapted from split_freq from VASP4CLINT, 21.04.2024
!
if (calc_freq_int) then
   allocate(dxyz(3,natoms))
   allocate(ddip(3,3,natoms))
   allocate(dipgrad(3,3,natoms))
   allocate(dip_list2(3,3*natoms))
!
!    Calculate the dipol derivatives
!

   do i=1,(natoms-fix_num)*3
      do j=1,natoms
         do k=1,3
            deltaxyz=int_pos_vecs(k,j,(i-1)*2+1)-int_pos_vecs(k,j,(i-1)*2+2)
            if (deltaxyz .gt. 1E-6) then
               dxyz(k,j) = deltaxyz
               do l=1,3
                  ddip(l,k,j)=dip_list(l,(i-1)*2+1)-dip_list(l,(i-1)*2+2)
                  dipgrad(l,k,j)=ddip(l,k,j) / dxyz(k,j)
               end do
            end if
         end do
      end do
   end do
!
!    Calculate normal mode intensities
!
   cnvint = 974.88d0
   allocate(intens(3,3*(natoms-fix_num)))
   allocate(totintens(3*(natoms-fix_num))) 
   intens=0.d0
   totintens=0.d0  
   do i=1,3*(natoms-fix_num)
      do j=1,3
         do k=1,natoms
            do l=1,3
               intens(j,i)=intens(j,i)+dipgrad(j,l,k)*n_vib2(l,k,i)
            end do
         end do
         intens(j,i)=intens(j,i)*intens(j,i)*cnvint
         totintens(i)=totintens(i)+intens(j,i)
      end do
   end do
end if
!
!     Shift structure to origin for better visualization with molden
!
do i=1,natoms
   coord(1,i)=x(i)
   coord(2,i)=y(i)
   coord(3,i)=z(i)
end do
com_coord=0.d0
do i=1,natoms
   do j=1,3
      com_coord(j)=com_coord(j)+mass(i)*coord(j,i)
   end do
end do
com_coord=com_coord/sum(mass(1:natoms))
do i=1,natoms
   coord(:,i)=coord(:,i)-com_coord(:)
end do

!
!     Write molden format output file for full system 
!
!
if (calc_freq_int) then
   open(unit=49,file="explore.molden",status="replace")
   write(49,*) "[Molden Format]"
   write(49,*) "[Atoms] Angs"
   do i=1,natoms
!   write(*,*) " ",name(i),i,elem_index(i),coord(:,i)
      write(49,'(a,a,i7,i4,3f12.6)') " ",name(i),i,elem_index(i), &
                & coord(:,i)
   end do
   write(49,*) "[FREQ]"
   do i=1,3*(natoms-fix_num)
      write(49,'(f10.4)') W(3*(natoms-fix_num)-i+1)
   end do
   write(49,*) "[INT]"

   if (sum(totintens) .lt. 1d0) then
      do i=1,3*(natoms-fix_num)
         write(49,'(f10.4)') 1.d0
      end do   
   else
      do i=1,3*(natoms-fix_num)
         write(49,'(f10.4)') totintens(i)
      end do
   end if
   write(49,*) "[FR-COORD]"
   do i=1,natoms
      write(49,'(a,a,3f11.6)') " ",name(i),coord(:,i)/bohr
   end do
   write(49,*) "[FR-NORM-COORD]"
   do i=1,3*(natoms-fix_num)
      write(49,'(a,i7)') "vibration",i
      ind=1
      do j=1,natoms
         write(49,'(3f11.6)') n_vib2(:,j,i)
      end do
   end do
   close(49)
   write(*,*)
   write(*,*) "File 'explore.molden' for mode visualization (all atoms) written!"
end if
!
!    return unchanged hessian
!
hess=hess_backup 

!
!    for return of frequencies
!
freqs=W
!
!   Now, if desired, calculate the normal mode cartesian displacements!!
!
eigvecs=A_Mat

 
deallocate(A_mat,W,work,n_vib,coord)
return
end subroutine calc_freq
