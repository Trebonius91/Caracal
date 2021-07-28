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
!     subroutine optimize_3evb: optimize coupling terms of 3x3 EVB
!
!     part of EVB
!
subroutine optimize_3evb(energies_qmdff,datapoints,&
           energies_result,fileenergy,filegeo)
use general
use evb_mod
implicit none
real(kind=8),dimension(datapoints,3)::energies_qmdff
real(kind=8),dimension(datapoints)::energies_reference,energies_result
real(kind=8),dimension(datapoints)::dq_dq1,dq_dq2
integer::datapoints,fileenergy_unit,nat
integer::i,j,k,m,par_num,maxstep,worse
real(kind=8),dimension(:),allocatable::beta_grad !the beta_vector(gradient)
real(kind=8),dimension(:,:),allocatable::alpha_hess !the alpha_matrix(hessian)
real(kind=8),dimension(:),allocatable::deltaA,deriv,num_step,upper,lower,sp
real(kind=8),dimension(:),allocatable::int_coord
real(kind=8)::Lambdas,DLambda
real(kind=8)::func
real(kind=8),dimension(3*n_one)::dq1,dq2
real(kind=8),dimension(3,n_one)::coord
real(kind=8)::chi,chi_old,e_evb,e_ref,stepsize,inc_chi
real(kind=8),dimension(3,3)::evb_mat !the well-known evb-matrix
logical::bad_run    ! decides, if a levenberg marquardt run should be broken off
! random number generation:
INTEGER SEED,t
INTEGER :: i_seed
INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
INTEGER, DIMENSION(1:8) :: dt_seed
real(kind=8)::r
! for multistart local search  
real(kind=8),dimension(:),allocatable::par_opt
real(kind=8)::glob_chi
integer::msls_step,p
character(len=60)::fileenergy,filegeo
logical::has_next
!
!     variables for lapack-dgesv-Subroutine
!
integer::N1,INFO
real(kind=8),dimension(:),allocatable::B
nat=natoms
fileenergy_unit=20
!
!     Set the number of optimizable parameters 
!
if (off_basis=="const") then
   par_num=3
else if (off_basis=="1g") then
   if (.not. full) then
      par_num=4
   else 
      par_num=6
   end if
else if (off_basis=="sd2") then
   if (.not. full) then
      par_num=12
   else
      par_num=14
   end if
end if
allocate(sp(par_num))
N1=par_num
allocate(B(N1))
!
!    other arrays
!
allocate(beta_grad(par_num))
allocate(alpha_hess(par_num,par_num))
allocate(deltaA(par_num))
allocate(deriv(par_num))
allocate(num_step(par_num),upper(par_num),lower(par_num))
!
!    if an EVB-DQ-calculation is desired, read in the needed coordinates
!    of the reactionpath and calculate the geometrical coordinate
!
if (evb_dq) then
   allocate(int_ts(3*nat-6))
   allocate(int_coord(3*nat-6))
   has_next=.true.
   open(unit=34,file=filegeo,status='old')
   do m=1,datapoints
      call next_geo(coord,natoms,34,has_next)
      call xyz_2int(coord,int_coord,nat)
      if (.not.has_next) exit
!     The geometrical coordinate dq=q-q_ts in internal coordinates
      int_ts=ts_coordinates
      do i=1,3*nat-6
         dq1(i)=int_coord(i)-int_ts(i)
         if (i .gt. 2*nat-3) then
            dq1(i)=sqrt((int_coord(i)-int_ts(i))*(int_coord(i)-int_ts(i))+&
                  &(int_coord(i+nat-3)-int_ts(i+nat-3))*&
                  &(int_coord(i+nat-3)-int_ts(i+nat-3)))
         end if
      end do
      int_ts=ts_coordinates2
      do i=1,3*nat-6
         dq2(i)=int_coord(i)-int_ts(i)
         if (i .gt. 2*nat-3) then
            dq2(i)=sqrt((int_coord(i)-int_ts(i))*(int_coord(i)-int_ts(i))+&
                  &(int_coord(i+nat-3)-int_ts(i+nat-3))*&
                  &(int_coord(i+nat-3)-int_ts(i+nat-3)))
         end if
      end do
      dq_dq1(m)=0
      dq_dq2(m)=0
      do i=1,3*n_one-6
         dq_dq1(m)=dq_dq1(m)+abs(dq1(i))/50
         dq_dq2(m)=dq_dq2(m)+abs(dq2(i))/50
      end do
   end do
   close(34)
end if

write(15,*) "########################################################"
write(15,*) "#           3x3-EVB : 1G COUPLINGTERM                  #"
write(15,*) "#  Doing the multi start local search alorithm in      #"
write(15,*) "#  combination with the Levenberg-Marquardt-Algorithm  #"
write(15,*) "#  to optimize the parameters of the off-diagonal term #"
write(15,*) "########################################################"
write(15,*)
write(15,*) "Your start-matrix looks like:"

if (off_basis == "const") then
   write(15,*) "(       E1             a                c          )"
   write(15,*) "(       a              E2               b          )"
   write(15,*) "(       c              b                E3         )"
end if
if (.not.evb_dq) then
   if (off_basis == "1g") then
      if (.not.full) then
         write(15,*) "(       E1          a*exp(-b*dE12^2)        0          )"
         write(15,*) "( a*exp(-b*dE12^2)       E2           c*exp(-d*dE23^2) )"
         write(15,*) "(       0           c*exp(-d*dE23^2)        E3         )"
      else
         write(15,*) "(       E1          a*exp(-b*dE12^2)  e*exp(-f*dE12)   )"
         write(15,*) "( a*exp(-b*dE12^2)        E2          c*exp(-d*dE23^2) )"
         write(15,*) "( e*exp(-f*dE23^2)  c*exp(-d*dE23^2)        E3         )"
      end if
   else if (off_basis == "sd2") then
      if (.not.full) then 
         write(15,*) "(                          a*exp(-b*dE12^2)                            )"
         write(15,*) "(       E1               +dE^2*c*exp(-d*dE12^2)          0             )"
         write(15,*) "(                        +dE^2*e*exp(-f*dE12^2)                        )"
         write(15,*) "(     a*exp(-b*dE12^2)                            g*exp(-h*dE23^2)     )"
         write(15,*) "( +dE^2*c*exp(-d*dE12^2)         E2             +dE^2*i*exp(-j*dE23^2) )"
         write(15,*) "( +dE^2*e*exp(-f*dE12^2)                        +dE^2*k*exp(-l*dE23^2) )"
         write(15,*) "(                          g*exp(-h*dE23^2)                            )"
         write(15,*) "(       0                +dE^2*i*exp(-j*dE23^2)          E3            )"
         write(15,*) "(                        +dE^2*k*exp(-l*dE23^2)                        )"
      else
         write(15,*) "(                          a*exp(-b*dE12^2)                            )"
         write(15,*) "(       E1               +dE^2*c*exp(-d*dE12^2)    m*exp(-n*dE13^2)    )"
         write(15,*) "(                        +dE^2*e*exp(-f*dE12^2)                        )"
         write(15,*) "(     a*exp(-b*dE12^2)                            g*exp(-h*dE23^2)     )"
         write(15,*) "( +dE^2*c*exp(-d*dE12^2)         E2             +dE^2*i*exp(-j*dE23^2) )"
         write(15,*) "( +dE^2*e*exp(-f*dE12^2)                        +dE^2*k*exp(-l*dE23^2) )"
         write(15,*) "(                          g*exp(-h*dE23^2)                            )"
         write(15,*) "(    m*exp(-n*dE13^2)    +dE^2*i*exp(-j*dE23^2)          E3            )"
         write(15,*) "(                        +dE^2*k*exp(-l*dE23^2)                        )"
     end if
   end if
else 
   if (off_basis == "1g") then
      write(15,*) "(       E1          a*exp(-b*dQ12^2)        0          )"
      write(15,*) "( a*exp(-b*dQ12^2)       E2           c*exp(-d*dQ23^2) )"
      write(15,*) "(       0           c*exp(-d*dQ23^2)        E3         )"
   else if (off_basis == "sd2") then
      write(15,*) "(                          a*exp(-b*dQ12^2)                            )"
      write(15,*) "(       E1               +dE^2*c*exp(-d*dQ12^2)          0             )"
      write(15,*) "(                        +dE^2*e*exp(-f*dQ12^2)                        )"
      write(15,*) "(     a*exp(-b*dQ12^2)                            g*exp(-h*dQ23^2)     )"
      write(15,*) "( +dE^2*c*exp(-d*dQ12^2)         E2             +dE^2*i*exp(-j*dQ23^2) )"
      write(15,*) "( +dE^2*e*exp(-f*dQ12^2)                        +dE^2*k*exp(-l*dQ23^2) )"
      write(15,*) "(                          g*exp(-h*dQ23^2)                            )"
      write(15,*) "(       0                +dE^2*i*exp(-j*dQ23^2)          E3            )"
      write(15,*) "(                        +dE^2*k*exp(-l*dQ23^2)                        )"
   end if
end if
write(15,*)
write(*,*)  " . ...- -... --.- -- -.. ..-. ..-."
open(unit=fileenergy_unit,file=fileenergy,status='unknown')
do i=1,datapoints
   read(fileenergy_unit,*) energies_reference(i)
end do
close(fileenergy_unit)

!
!   Generate random start numbers for all parameters
!   (taken from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html)
!

glob_chi=0d0
msls_step=0
do p=1,maxstart
   bad_run=.false.
   do i=1,par_num
      CALL RANDOM_SEED(size=i_seed)
      ALLOCATE(a_seed(1:i_seed))
      CALL RANDOM_SEED(get=a_seed)
      CALL DATE_AND_TIME(values=dt_seed)
      a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
      CALL RANDOM_SEED(put=a_seed)
      DEALLOCATE(a_seed)
! ----- Done setting up random seed -----

      CALL RANDOM_NUMBER(r)
!
!  Interval, in which the number should be
!
      r=r*(upper_bond-lower_bond)+lower_bond
      sp(i)=r
   end do

!
!     estimate the starting-point
!
   chi=0
   open(unit=123,file="test.plot",status='unknown')
   do i=1,datapoints
      call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                 &energies_qmdff(i,3),energies_reference(i),sp,&
                 &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))
      chi=chi+inc_chi
      write(123,*) E_evb
   end do
   write(15,*) "Start-parameters:"
   write(15,*) sp(1:par_num)
   write(15,*) "Fitness in the starting-point:",chi
   close(123)


   stepsize=diff_step !size of the numdiff-step (in fraction of the parameter itself)
   Lambdas=lm_par 
   maxstep=1

   do
      write(15,*) " "
      write(15,'(A,I3)') "STEP",maxstep
      write(15,*) " "
     !reset gradient and hessian   
      do i=1,par_num
         beta_grad(i)=0
         upper(i)=0
         lower(i)=0
      end do

      do i=1,par_num
         do j=1,par_num
            alpha_hess(i,j)=0
         end do
      end do

!
!    some parameters are very small,others are rather big:
!    make the num-diff step-size dependant of the parameter-value
!   
      do i=1,par_num
         num_step(i)=sp(i)*stepsize
      end do
!
!    construct the beta-vector (numerical first derivatives/gradient)
!
      do i=1,datapoints
         E_ref=energies_reference(i)
         do j=1,par_num
!    increment
            sp(j)=sp(j)+num_step(j)
            call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                    &energies_qmdff(i,3),energies_reference(i),sp,&
                    &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))
            upper(j)=E_evb!(E_ref-E_evb)*(E_ref-E_evb)
!    decrement
            sp(j)=sp(j)-2d0*num_step(j)
            call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                    &energies_qmdff(i,3),energies_reference(i),sp,&
                    &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))
            lower(j)=E_evb!(E_ref-E_evb)*(E_ref-E_evb)    
!    reset parameter and calculate numerical derivative  
            sp(j)=sp(j)+num_step(j)
            deriv(j)=(upper(j)-lower(j))/(2*num_step(j))
         end do
         call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                    &energies_qmdff(i,3),energies_reference(i),sp,&
                    &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))
         func=E_evb
!
!     gradient-vector
!
         do j=1,par_num
            beta_grad(j)=beta_grad(j)+((E_ref-func)*deriv(j))
         end do
!
!     hessian matrix
!
         do j=1,par_num
            do k=1,par_num
               alpha_hess(k,j)=alpha_hess(k,j)+deriv(j)*deriv(k)
            end do
         end do
      end do
!
!     Manipulate diagonal elements: add Lambda
!
      do j=1,par_num
         alpha_hess(j,j)=alpha_hess(j,j)+Lambdas*alpha_hess(j,j)
      end do

!
!     call the DGL-solving routine 
!
      call solving_lgs(par_num,alpha_hess,beta_grad,info,deltaA)  
      if (abs(info-0) >= 1E-8) then
         bad_run=.true.
      end if

      do i=1,par_num
         sp(i)=sp(i)+deltaA(i)
      end do
      chi_old=chi
      chi=0
!calculate the new fitness
      do i=1,datapoints
         call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                       &energies_qmdff(i,3),energies_reference(i),sp,&
                       &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))

         chi=chi+inc_chi
      end do
      if (chi >= chi_old) then
         Lambdas=Lambdas*lm_par_change
         do i=1,par_num
            sp(i)=sp(i)-deltaA(i)
         end do
         worse=1
      else
!      Lambda=Lambda/lm_par_change
         worse=0
      end if

!  prevents that parameters grows or sinks too much
!
      do i=1,par_num
         if (sp(i).ge.10000 .or. sp(i).le.1E-7) then
            goto 10
         end if
      end do
      goto 12
10    Lambdas=Lambdas*lm_par_change
      do i=1,par_num
         sp(i)=sp(i)-deltaA(i)
      end do
      worse=1
12    maxstep=maxstep+1
      if (abs(chi_old-chi) < lm_threshold .or. abs(chi_old-chi)/chi< lm_threshold/100 &
         .or. maxstep>optsteps .or. bad_run) exit
 
      if (worse==1) then
         chi=chi_old
      end if
      write(15,'(A,(g25.18))') "The current fitness:",chi
      write(15,*) "The current parameter-values:" 
      write(15,*) sp(1:par_num)
   end do
!
!   if the factorization has failed
!
   if (bad_run) then
      maxstep=optsteps
   end if

!
!     Multi start local search: only if optimization has converged: 
!     if the actual fitness value is better than the best before, replace 
!     it with the actual value
!
112 if (maxstep<optsteps) then
      if (glob_chi .eq. 0) then
         glob_chi=chi
         par_opt=sp
         write(*,*) "After Multi-Start-Local-Search step:",msls_step
         write(*,*) "Local optimization finished, new best fitness:",glob_chi
         write(15,*) "------------------------------"
         write(15,*) "Local optimization finished!"
         write(15,*) "First best fitness value:",glob_chi
         write(15,*) "------------------------------"
      else if (chi .le. glob_chi) then
         glob_chi=chi
         par_opt=sp
         write(*,*) "After Multi-Start-Local-Search step:",msls_step
         write(*,*) "Local optimization finished, new best fitness:",glob_chi
         write(15,*) "------------------------------"
         write(15,*) "Local optimization finished!"
         write(15,*) "New best fitness value:",glob_chi
         write(15,*) "------------------------------"
      else
         write(15,*) "------------------------------"
         write(15,*) "Local optimization finished!"
         write(15,*) "No better value found..."
         write(15,*) "------------------------------"
      end if
   else
      write(15,*) "The optimization didn´t converge in the given cycles."
      write(15,*) "The result is therefore ignored."
      write(15,*)
   end if
   msls_step=msls_step+1
end do
if (glob_chi .eq. 0) then
   write(*,*) "I´m sorry but the optimization could not archive a useful result..."
   write(*,*) "Please retry it or change the settings."
   call fatal
else
   write(*,*)  ".. .----. -- -.. --- -. . "
   write(*,*) "Multi start local optimization finished!"
   write(*,*) "The best fitness-value is:", glob_chi
   write(*,*) "The optimized parameters are:"
   write(*,*) par_opt
   write(*,*) "Look into evbopt.log for further details."
   write(15,*) "Multi start local optimization finished!"
   write(15,*) "The best fitness-value is:", glob_chi
   write(15,*) "The optimized parameters are:"
   write(15,*) par_opt
   write(15,*)
end if

!
!   At last, calculate the path energies 
!   for the optimized alpha and b_vec values
!
if (maxstart .gt. 1) then
   sp=par_opt
end if

do i=1,datapoints
   call energy_1g(energies_qmdff(i,1),energies_qmdff(i,2),&
                 &energies_qmdff(i,3),energies_reference(i),sp,&
                 &inc_chi,par_num,E_evb,dq_dq1(i),dq_dq2(i))
   energies_result(i)=E_evb
end do

!
!   deallocation
!

deallocate(beta_grad)
deallocate(alpha_hess)
deallocate(deltaA)
return

return
end subroutine optimize_3evb
