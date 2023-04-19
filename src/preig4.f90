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
!     subroutine preig4: print orbital energies and accupancies 
!          for Extended Hückel in QMDFF generation
!
!     part of QMDFF
!

subroutine preig4(IO,OCC,E,NORBS)
IMPLICIT REAL*4 (A-H,O-Z)
DIMENSION E(*), OCC(*)

write(io,'(/,10x,''eigenvalues'')')
N=6
NTIMES=NORBS/N
NREST=MOD(NORBS,N)
IF (NTIMES.EQ.0) NREST=NORBS

J=1
N2=N

DO K=1,NTIMES
   WRITE(IO,100)(I,I=J,N2)
   WRITE(IO,200)(occ(i),I=J,N2)
   WRITE(IO,300)(E(i),I=J,N2)
   J =J +N
   N2=N2+N
END DO

IF (NREST.GT.0.OR.NTIMES.EQ.0) THEN
   WRITE(IO,100)(I,I=J,J+NREST-1)
   WRITE(IO,200)(occ(i),I=J,J+NREST-1)
   WRITE(IO,300)(E(i),I=J,J+NREST-1)
END IF
 
 100  FORMAT(' #    : ',2X,6(4X,I4,2X))
 200  FORMAT(' occ. : ',2X,6(4X,F5.3,1X))
 300  FORMAT(' eps  : ',2X,6F10.5)

return
end subroutine preig4
