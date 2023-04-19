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
!     subroutine preig2: print frequencies of hessian matrices
!
!     part of QMDFF
!
subroutine preig2(IO,E,NORBS)
IMPLICIT REAL*8 (A-H,O-Z)
DIMENSION E(*)

N=6
NTIMES=NORBS/N
NREST=MOD(NORBS,N)
IF (NTIMES.EQ.0) NREST=NORBS
J=1
N2=N
DO K=1,NTIMES
   WRITE(10,300)J,N2,(E(i),I=J,N2)
   J =J +N
   N2=N2+N
END DO
IF (NREST.GT.0.OR.NTIMES.EQ.0) THEN
   WRITE(10,300)J,J+NREST-1,(E(i),I=J,J+NREST-1)
END IF

100  FORMAT(' #    : ',2X,6(4X,I4,2X))
300  FORMAT(' freq.',i4,'-',i4,':',2X,10F9.2)

return
end subroutine preig2
