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
!     subroutine getpar: fill QMDFF hessian with values of bonded 
!                  energy terms
!
!     part of QMDFF
!
subroutine getpar(p,torsfit,hessatoms)
use qmdff
implicit none
logical::torsfit
real(kind=8)::p(*)
integer::i,j,k
integer::hessatoms(5,*)

k=0
do i=1,ntors
   if (tors(6,i).eq.0) cycle
   k=k+1
   p(k)=vtors(2,i)
   hessatoms(1,k)=tors(1,i)
   hessatoms(2,k)=tors(2,i)
   hessatoms(3,k)=tors(3,i)
   hessatoms(4,k)=tors(4,i)
   hessatoms(5,k)=4
end do
do i=1,nbond
   k=k+1
   p(k)=vbond(2,i)
   hessatoms(1,k)=bond(1,i)
   hessatoms(2,k)=bond(2,i)
   hessatoms(5,k)=2
end do
do i=1,nangl
   k=k+1
   p(k)=vangl(2,i)
   hessatoms(1,k)=angl(1,i)
   hessatoms(2,k)=angl(2,i)
   hessatoms(3,k)=angl(3,i)
   hessatoms(5,k)=3
end do

return
end subroutine getpar
