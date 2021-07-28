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
!     subroutine diagonalize: calls the LAPACK dsyev routine for 
!       diagonalizing matrix
!
!     part of EVB
!
subroutine diagonalize(Nd,A,W)
implicit none
character(len=1)::JOBZ,UPLO
integer::Nd,LDA,LWORK,INFO,i
real(kind=8)::A(Nd,Nd),W(Nd)
real(kind=8),dimension(:),allocatable::WORK
JOBZ='V'
UPLO='U'
LDA=Nd
INFO=0
LWORK=Nd*Nd-1
allocate(WORK(LWORK))
call dsyev(JOBZ,UPLO,Nd,A,LDA,W,WORK,LWORK,INFO)
if (info .ne. 0) then
   write(*,*) "The matrix diagonalization in diagonalize.f90 failed!"
   write(*,*) "The returned Lapack-INFO code was: ", INFO
   call fatal
end if
!
!    fill matrix with eigenvalues on diagonals
!
A=0.d0
do i=1,Nd
   A(i,i)=W(i)
end do

return
end subroutine diagonalize
