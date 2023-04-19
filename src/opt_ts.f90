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
!     Subroutine opt_ts: optimize a transition state structure using 
!     the TS optimization algorithm described by Baker (als also used 
!     in orca 4.0)
!
!     part of EVB
!
subroutine opt_ts
use general
use evb_mod
implicit none 

integer::i,j   ! loop indices
integer::niter  ! geoopt steps
integer::atind(natoms) ! array with element numbers 
integer::nlines,nrest  ! print out of normal mode frequencies
real(kind=8)::xyz2(3,natoms)   ! local structure array
real(kind=8)::xyz_1d(3*natoms)  ! 1D structure array for algorithm
real(kind=8)::energy  ! the current energy
real(kind=8)::grad(3,natoms)  ! the current gradient vector
real(kind=8)::grad_1d(3*natoms)  ! one dimensional gradient vector
real(kind=8)::deltq(3,natoms)   ! cartesian optimization step
real(kind=8)::eigval_cart(3*natoms) ! cartesian eigenvalues
real(kind=8)::eigvecs(3*natoms,3*natoms)  ! cartesian eigenvectors
real(kind=8)::hess(3*natoms,3*natoms)  ! hessian matrix
real(kind=8)::hess_inv(3*natoms,3*natoms)  ! inverse of the hessian (test) 
real(kind=8)::units(nat6,nat6)  ! test unit matrix


real(kind=8)::F_bar(nat6)  ! projected eigenvectors
real(kind=8)::mass3(3*natoms)  ! test for mass weighted hessian..
real(kind=8)::crit_abs ! the first convergence criterion: maximal abs value
real(kind=8)::steplen  ! actual length of geometry change vector
real(kind=8)::gradnorm  ! length of actual gradient vector
real(kind=8)::steplenmax,gradmax  ! largest components of step/gradient vectors
! internal coordinates
real(kind=8)::grad_int(nat6),hess_int(nat6,nat6)
real(kind=8)::deltq_int(nat6),coord_int(nat6),coord_int2(nat6)
real(kind=8)::hess_int_inv(nat6,nat6) ! inverse of the internal hessian
real(kind=8)::unit_mat(nat6,nat6) ! unit matrix for test reasons
! for matrix diagonalization
character(len=1)::JOBZ,UPLO
integer::M,Nn,LDA,LDU,LDV,LWORK,INFO
integer,dimension(:),allocatable::ipiv
integer::imags  ! number of imaginary eigenvalues
real(kind=8),dimension(:),allocatable::WORK,Eigval
real(kind=8),dimension(:,:),allocatable::V_mat
real(kind=8)::step(3*natoms)
logical::eig_use(3*natoms)   ! which eigenvalue is large enough
real(kind=8)::lam_n,lam_p,lam_n_old  ! shift parameter (P-RFO)
!
!     Read in the start structure from global array
!     and convert it into bohr
!
do i=1,natoms
   xyz2(1,i)=x(i)/bohr
   xyz2(2,i)=y(i)/bohr
   xyz2(3,i)=z(i)/bohr
end do
!allocate(indi(natoms))
!
!     allocate optional indices for fractionated optimization
! 
do i=1,natoms
   indi(i)=i
end do
!
!     allocate needed lapack arrays
!
JOBZ='V' !eigenvalues and eigenvectors(U)
UPLO='U' !upper triangle of a
Nn=nat6
LDA=Nn
INFO=0
LWORK=3*Nn-1

allocate(V_mat(Nn,Nn))
allocate(Eigval(Nn))
allocate(WORK(LWORK))
allocate(ipiv(Nn))

!
!     (1) Calculate gradient and hessian for the input structure
!
open(unit=23,file="ts_traj.xyz",status="unknown")
!
!     Write start structure to trajectory file
!
write(23,*) natoms 
write(23,*) "TS-optimization, initial structure "
do i=1,natoms
   write(23,*) name(i),xyz2(:,i)*bohr
end do
write(15,*) " --- TS-OPTIMIZATION: --- "
write(15,*) 
if (newton_raphson) then
   write(15,*) " The simple Newton-Raphson optimization method will be used."
else 
   write(15,*) " The P-RFO optimization method as implemented in orca will be used."
end if
write(15,*) 
do niter=1,maxiter 
   write(15,'(a,i6)') " * Geoopt step No. ",niter
!
!     Convert the actual structure to internal coordinates
!
   call xyz_2int(xyz2,coord_int,natoms)
!
!     calculate actual cartesian gradient
!
   call gradient(xyz2,energy,grad,1)
!   write(*,*) "grad_xyz",grad
!
!     convert gradient to internal coordinates
!  
   call grad2int(xyz2,coord_int,grad_int,grad)
!   write(*,*) "grad_int",grad_int
!
!     calculate actual cartesian hessian 
!
   call hessevb(xyz2,hess)
!
!     Calculate frequencies in cartesian coordinates for information
!
   call calc_freq(hess,eigval_cart,eigvecs,.false.)
!
!     convert hessian to internal coordinates
!
   call hess2int(xyz2,coord_int,hess_int,hess,grad_int)
!
!     (2) Diagonalize the hessian in order to determine local 
!         surface characteristics
!
   V_mat=hess_int
!   write(*,*) "internal hessian:"
!   do i=1,nat6
!      write(*,'(6e17.8)') hess_int(i,:)
!   end do
!   stop "Hhfff"
   deallocate(work)
   allocate(WORK(LWORK))
   call DSYEV(JOBZ,UPLO,Nn,V_mat,LDA,Eigval,WORK,LWORK,INFO)
   if (info .ne. 0) then
      write(*,*) "A diagonalization of a hessian matrix in opt_ts.f90 failed!"
      call fatal
   end if
!
!     A: The P-RFO method for explicit splitting of positive/negative
!        eigenvalues/vectors
!
!     Set small eigenvalues to zero
!     If they were set to zero, do not use their eigenvectors for further 
!     calculations
! 
   if (.not. newton_raphson) then
      imags=0
      eig_use=.true.
      do i=1,nat6
         if (abs(Eigval(i)) .lt. 1E-8) then
            Eigval(i)=0.d0
            eig_use(i)=.false.
         end if
!
!     Then look how many eigenvalues are greater than zero
!
         if (Eigval(i) .lt. 0.d0) then
            imags=imags+1
            eig_use(i)=.false.
         end if 
      end do
      write(15,'(a,i4)') "  Number of (internal) imaginary modes in the actual strucure: ",imags
!
!     Print the actual eigenvalues in mass weighted coordinates for test reasons
!
      nlines=(3*natoms)/6
      nrest=mod(3*natoms,6)
      write(15,*) " The cartesian normal modes in units of cm^-1:"
      write(15,*) "     ------------     "
      do i=1,nlines
         write(15,'(6f13.6)') eigval_cart((i-1)*6+1:i*6) 
      end do
      if (nrest .ne. 0) then
         write(15,'(5f13.6)') eigval_cart(nlines*6+1:3*natoms)
      end if
      write(15,*) "     ------------     "

!
!     Transform gradient vector from cartesian space into normal mode space
!
      F_bar=matmul(V_mat,grad_int)
!
!     (3) determine the step h, that shall be taken next in coordinate
!         space in order to reach the TS structure 
!
!     The Lambda_p for the single negative eigenvalue 
!
      lam_p=0.5d0*Eigval(1)+0.5*sqrt(Eigval(1)*Eigval(1)+4.d0*F_bar(1)*F_bar(1))
!
!     The Lambda_n: one global value for all positive eigenvalues
!
      lam_n_old=1.d0
      j=1
      do 
         lam_n=0.d0
         do i=2,nat6
            if (eig_use(i)) then
               lam_n=lam_n+F_bar(i)*F_bar(i)/(lam_n_old-Eigval(i))
            end if 
         end do
         if (abs(lam_n-lam_n_old) .lt.1D-4) exit
         lam_n_old=lam_n
         j=j+1
      end do
   
!
!     Set together the geometry step to be taken:
!     deltq=-F1*V1/(b1-lam_p)-sum_{i=2}^N(F_i*V_i/(b_i-lam_n))
!
      deltq_int=0.d0
      deltq_int=-F_bar(1)*V_mat(:,1)/(Eigval(1)-lam_p)

      do i=2,nat6
         if (eig_use(i)) then
            deltq_int=deltq_int-F_bar(i)*V_mat(:,i)/(Eigval(i)-lam_n)
         end if
      end do
!
!     B) The usual Newton-Raphson method in order to converge to the 
!        nearest extremum (might be TS or minimum!)
!        --> formula similar to above, but no Lambda-correction parameters
!   
   else 
!
!     Formula: h=-V_i^T*g*V_i/b_i
!
      deltq_int=0.d0

!
!     Print the actual eigenvalues in mass weighted coordinates for test reasons
!
      nlines=(3*natoms)/6
      nrest=mod(3*natoms,6)
      write(15,*) " The cartesian normal modes in units of cm^-1:"
      write(15,*) "     ------------     "
      do i=1,nlines
         write(15,'(6f13.6)') eigval_cart((i-1)*6+1:i*6)
      end do
      if (nrest .ne. 0) then
         write(15,'(5f13.6)') eigval_cart(nlines*6+1:3*natoms)
      end if
      write(15,*) "     ------------     "

!
!     If the eigenvalue is too small, neglect the fraction of the step belonging 
!     to this one
!
!      do i=1,nat6
!         if (abs(eigval(i)) .gt. 1E-8) then
!            write(*,*) "eigval:",i,eigval(i)
!            deltq_int=deltq_int-dot_product(V_mat(:,i),grad_int)*V_mat(i,:)/Eigval(i)
!         end if
!      end do
!      write(*,*)"eigenval", deltq_int
!
!     TEST: direct application of formula: h = -H^-1*g
!
      call pseudoinv(hess_int,nat6,nat6,hess_int_inv)

      !units=matmul(hess_int,hess_int_inv)
      deltq_int=-matmul(hess_int_inv,grad_int)
!      write(*,*) "direct",deltq_int
!      stop

   end if
! 
!     Convert the resulting step to bohr!
!
   deltq_int=deltq_int*bohr
!
!     Calculate the absolute value/length of the optimization step vector
!
   steplen=sqrt(dot_product(deltq_int,deltq_int))
   gradnorm=sqrt(dot_product(grad_int,grad_int))
   gradmax=maxval(abs(grad_int))
   steplenmax=maxval(abs(deltq_int))
!
!     Print actual status to logfile
!
   write(15,'(a,e18.11,a)') "  Actual energy: ",energy," Eh."
   if (steplen .lt. dthr) then
      write(15,'(a,e11.4,a,e11.4,a)') "  Length of the coordinate step:      ", &
                  & steplen," (tol:",dthr,") -->  OK! "
   else 
      write(15,'(a,e11.4,a,e11.4,a)') "  Length of the coordinate step:      ", &
                  & steplen," (tol:",dthr,") -->  X "
   end if
   if (gradnorm .lt. gthr)  then
      write(15,'(a,e11.4,a,e11.4,a)') "  Norm of the actual gradient vector: ", &
                  & gradnorm," (tol:",gthr,") -->  OK! "
   else 
      write(15,'(a,e11.4,a,e11.4,a)') "  Norm of the actual gradient vector: ", &
                  & gradnorm," (tol:",gthr,") -->  X "
   end if
   if (steplenmax .lt. dmaxthr)  then
      write(15,'(a,e11.4,a,e11.4,a)') "  Largest component of coord. step  : ", &
                  & steplenmax," (tol:",gthr,") -->  OK! "
   else 
      write(15,'(a,e11.4,a,e11.4,a)') "  Largest component of coord. step  : ", &
                  & steplenmax," (tol:",gthr,") -->  X "
   end if
   if (gradmax .lt. gmaxthr) then
      write(15,'(a,e11.4,a,e11.4,a)') "  Largest component of grad. vector : ", &
                  & gradmax," (tol:",gmaxthr,") -->  OK! "
   else 
      write(15,'(a,e11.4,a,e11.4,a)') "  Largest component of grad. vector : ", &
                  & gradmax," (tol:",gmaxthr,") -->  X "
   end if

!
!     Signal convergence if all four criterions are fulfilled!
!
   if ((steplen .lt. dthr) .and. (gradnorm .lt. gthr) .and. &
      &  (steplenmax .lt. dmaxthr) .and. (gradmax .lt. gmaxthr)) then
      write(*,*) "All convergence criterions are fulfilled!"
      write(*,*) "Convergence will therefore be signaled now!"
      write(*,*) 
      write(15,*) 
      write(15,*) "All convergence criterions are fulfilled!"
      write(15,*) "CONVERGENCE will therefore be signaled now!"
      write(15,*)
      exit
   end if
!
!     scale down the size of the step, if its larger than the tolerance DMAX
!
   if (steplen .gt. stepmax) then
      deltq_int=deltq_int*stepmax/steplen
   end if
!   deltq_int=deltq_int*0.1!05
!
!     Add step in internal coordinates to cartesian structure! 
!     Do this in an iterative process due to curvilinearity of 
!     the internal coordinate set.
!
   call int2step(xyz2,coord_int,deltq_int,deltq)

!
!     Convert step vector back to cartesian coordinates
!
!
!     Apply step to structure!
!
   xyz2=xyz2+deltq
   call xyz_2int(xyz2,coord_int2,natoms)

!   write(*,'(a,6e19.10)') "coord1",-coord_int+coord_int2
!   write(*,'(a,6e19.10)') "coord2",deltq_int
!
!     Write actual structure to file
!
   write(23,*) natoms
   write(23,*) "TS-optimization, step ",niter
   do i=1,natoms
      write(23,*) name(i),xyz2(:,i)*bohr
   end do
   write(15,*)  
end do
!
!     Signal if convergence couln't be reached and abort with error!
!
if (niter .gt. maxiter) then
   write(*,*) "Convergence for the TS-opt could't be reached! Raise the number of maximal"
   write(*,*) "steps by adding OPT_MAXITER or change the start structure!"
   call fatal
end if
!
!     estimate atom numbers from element symbols
!
do i=1,natoms
   call elem(name(i),atind(i))
end do
!
!     Write a molden file for inspectation of normal modes 
!
call calc_freq(hess,eigval_cart,eigvecs,.false.)
call g98fake("ts_opt.molden",natoms,atind,xyz2,eigval_cart,eigvecs,eigvecs)
write(*,*) "File with normal mode vibrations and IR spectrum was written to"
write(*,*) "ts_opt.molden. Open it with molden to inspect the normal modes!"
write(*,*)
write(15,*) "File with normal mode vibrations and IR spectrum was written to"
write(15,*) "ts_opt.molden. Open it with molden to inspect the normal modes!"
write(15,*) 

!
!     Write the final converged structure to file 
!
write(23,*) natoms
write(23,*) "TS-optimization: Final structure, converged!"
do i=1,natoms
   write(23,*) name(i),xyz2(:,i)*bohr
end do

write(*,*) "Succession of optimized structures written to 'ts_traj.xyz'"
write(*,*) "Optimized TS structure written to file 'ts_opt.xyz'."
open(unit=39,file="ts_opt.xyz",status="replace")

write(39,*) natoms
write(39,*) "TS-optimization: Final structure, converged!"
do i=1,natoms
   write(39,*) name(i),xyz2(:,i)*bohr
end do

close(39)

!
!     Overwrite global xyz coordinate array with optimized structure 
!
do i=1,natoms
   x(i)=xyz2(1,i)*bohr
   y(i)=xyz2(2,i)*bohr
   z(i)=xyz2(3,i)*bohr
end do



return
end subroutine opt_ts
