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
!     subroutine calc_freq: performs a normal coordinate analysis on a given 
!     hessian matrix
!
!     part of EVB
!

subroutine calc_freq(hess,freqs,eigvecs,printfreq)
use evb_mod
use general

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
integer::NRHS,LDB
logical::printfreq ! if normal mode frequencies shall be printed to log
real(kind=8)::pi
real(kind=8)::factor
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
Nn=3*nats
LDA=Nn
INFO=0
LWORK=Nn*Nn-1

allocate(A_mat(Nn,Nn))
allocate(W(Nn))
allocate(WORK(LWORK))
allocate(n_vib(Nn,nats,3))
allocate(coord(nats,3))
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
do i = 1, 3*nats
   W(i)=au2cm*sign(sqrt(abs(W(i))),W(i))/sqrt(amu2au)
end do
if (printfreq) then
   write(15,*)"*The normal mode frequencies (cm-1) of the structure are:"
   write(15,*)
   do i=1,3*nats
      write(15,*) i, W(i)
   end do
   write(15,*)
end if
!
!    for return of frequencies
!
freqs=W
!
!   Now, if desired, calculate the normal mode cartesian displacements!!
!   maximal number: 99 normal modes
!
eigvecs=A_Mat
if (print_nm) then
   call system("rm -r normal_modes")
   call system("mkdir normal_modes")
   do i=1,nat3
!
!   define output file
!
      if (i .lt. 10) then
         allocate(character(len=5)::filename1)
         write(filename1,"(I1,A4)") i,".xyz"
      else 
         allocate(character(len=6)::filename1)
         write(filename1,"(I2,A4)") i,".xyz"
      end if
      open(unit=87,file=filename1,status="unknown")
!
!   load input structure
!
      do j=1,nats
         coord(j,1)=x(indi(j))
         coord(j,2)=y(indi(j))
         coord(j,3)=z(indi(j))
      end do
!
!   define normal mode displacement vector (factor 1/2 for better view)
!
      do j=1,nats
         do k=1,3
            n_vib(i,j,k)=A_mat((j-1)*3+k,i)/2
         end do
      end do
      coord=coord-n_vib(i,:,:)
!
!   write normal modes (20 frames) to output file
!
      do j=1,20
         write(87,*) nats
         write(87,*) 
         do k=1,nats
           write(87,*) name(indi(k)),coord(k,:)
         end do
         coord=coord+n_vib(i,:,:)/10
      end do
      do j=1,20
         write(87,*) nats
         write(87,*)
         do k=1,nats
           write(87,*) name(indi(k)),coord(k,:)
         end do
         coord=coord-n_vib(i,:,:)/10
      end do

      close(87)
      write(syscall,"(A3,A6,A13)") "mv ",filename1," normal_modes"
      call system(syscall)
      deallocate(filename1)
   end do
   write(*,*)"Normal mode vibrations are written to folder normal_modes"
   write(15,*)"*Normal mode vibrations are written to folder normal_modes"
   write(15,*)
end if
!
!    return unchanged hessian
!
hess=hess_backup 
 
deallocate(A_mat,W,work,n_vib,coord)
return
end subroutine calc_freq
