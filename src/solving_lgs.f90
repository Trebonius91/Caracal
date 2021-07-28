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
!     subroutine solving_lgs: solve a quadratic system of 
!     linear equations 
!
!     part of EVB
!
subroutine solving_lgs(N,A,B,INFO,deltaA)
!calls the lapack-subroutine for solving LGS
integer::N,NRHS,LDA,LDB,INFO
real(kind=8),dimension(N,N)::A
real(kind=8),dimension(N)::B,deltaA
integer,dimension(N)::IPIV
NRHS=1
LDA=N
LDB=N
call dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)

if (abs(INFO-0) >= 1E-8) then
   write(15,*) "Error, the factorization of the Hessian failed!"
   write(15,*) "Aborting the calculation ..."
   write(15,*) "Start new local optimization"
!   write(*,*) "ERROR: No convergence could be archieved!"
!   write(*,*) "Open Logfile 'evbopt.log' for details."
!   call fatal
end if
deltaA=B
return
end subroutine solving_lgs

