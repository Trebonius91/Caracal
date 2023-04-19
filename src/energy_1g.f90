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
!     subroutine energy_1g: calculates the fitness for one structure in the 
!     3x3-EVB-couplingmode
!
!     part of EVB
!  
subroutine energy_1g(E1,E2,E3,E_ref,sp,inc_chi,par_num,E_evb,dQ12,dQ23)
use evb_mod
implicit none
integer::par_num
real(kind=8)::E1,E2,E3,E_ref,inc_chi,dQ12,dQ23
real(kind=8)::dE12,dE13,dE23,E_evb
real(kind=8),dimension(14)::sp
!
!     variables for lapack-dsyev-Subroutine
!
character(len=1)::JOBZ,UPLO
integer::INFO,LDA,LWORK,N
real(kind=8),dimension(:,:),allocatable::A,A1
real(kind=8),dimension(:),allocatable::W,WORK
real(kind=8),dimension(3,3)::evb_mat
JOBZ='N' !only eigenvalues
UPLO='U' !upper triangle of a
N=3
LDA=N
INFO=0
LWORK=N*N-1
allocate(A(N,N))
allocate(W(N))
allocate(WORK(LWORK))

evb_mat(1,1)=E1
evb_mat(2,2)=E2
evb_mat(3,3)=E3

if (evb_dq) then
   if (off_basis=="1g") then
      evb_mat(1,2)=sp(1)*exp(-sp(2)*dQ12*dQ12)
      evb_mat(2,1)=evb_mat(1,2)
      evb_mat(2,3)=sp(3)*exp(-sp(4)*dQ23*dQ23)
      evb_mat(3,2)=evb_mat(2,3)
      evb_mat(1,3)=0
      evb_mat(3,1)=0
   else if (off_basis=="sd2") then
      evb_mat(1,2)=sp(1)*exp(-sp(2)*dQ12*dQ12)+sp(3)*dQ12*dQ12*&
                  &exp(-sp(4)*dQ12*dQ12)+sp(5)*dQ12*dQ12*exp(-sp(6)*dQ12*dQ12)
      evb_mat(2,1)=evb_mat(1,2)
      evb_mat(2,3)=sp(7)*exp(-sp(8)*dQ23*dQ23)+sp(9)*dQ23*dQ23*&
                  &exp(-sp(10)*dQ23*dQ23)+sp(11)*dQ23*dQ23*exp(-sp(12)*dQ23*dQ23)
      evb_mat(3,2)=evb_mat(2,3)
      evb_mat(1,3)=0
      evb_mat(3,1)=0
   end if
else
   dE12=abs(E1-E2)
   dE13=abs(E1-E3)
   dE23=abs(E2-E3)
   evb_mat(1,1)=E1
   evb_mat(2,2)=E2
   evb_mat(3,3)=E3
   if (off_basis=="1g") then 
      evb_mat(1,2)=sp(1)*exp(-sp(2)*dE12*dE12)
      evb_mat(2,1)=evb_mat(1,2)
      evb_mat(2,3)=sp(3)*exp(-sp(4)*dE23*dE23)
      evb_mat(3,2)=evb_mat(2,3)
      if (full.eqv. .true.) then
         evb_mat(1,3)=sp(5)*exp(-sp(6)*dE13*dE13)
         evb_mat(3,1)=evb_mat(1,3)
      else 
         evb_mat(1,3)=0
         evb_mat(3,1)=0
      end if
   else if (off_basis=="sd2") then
      evb_mat(1,2)=sp(1)*exp(-sp(2)*dE12*dE12)+sp(3)*dE12*dE12*&
                  &exp(-sp(4)*dE12*dE12)+sp(5)*dE12*dE12*exp(-sp(6)*dE12*dE12)
      evb_mat(2,1)=evb_mat(1,2)
      evb_mat(2,3)=sp(7)*exp(-sp(8)*dE23*dE23)+sp(9)*dE23*dE23*&
                  &exp(-sp(10)*dE23*dE23)+sp(11)*dE23*dE23*exp(-sp(12)*dE23*dE23)
      evb_mat(3,2)=evb_mat(2,3)
      if (full .eqv. .true.) then
         evb_mat(1,3)=sp(13)*exp(-sp(14)*dE13*dE13)
         evb_mat(3,1)=evb_mat(1,3)
      else
         evb_mat(1,3)=0
         evb_mat(3,1)=0
      end if   
   end if
end if
A=evb_mat
call DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
E_evb=W(1)
inc_chi=(E_ref-E_evb)*(E_ref-E_evb)
deallocate(A,W,WORK)
return
end subroutine energy_1g
