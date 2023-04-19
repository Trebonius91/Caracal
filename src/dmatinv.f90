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
!     subroutine dmatinv: invert a quadratic matrix, needed 
!        during QMDFF fit
!
!     part of QMDFF
!

subroutine dmatinv (A,LDM,N,D)
IMPLICIT NONE
!
!     *** Deklaration globaler Variablen ***
!
INTEGER, INTENT(IN) :: LDM, N
REAL*8, INTENT(OUT)   :: D
REAL*8, INTENT(INOUT) :: A(LDM,*)
!
!     *** Deklaration lokaler Variablen ***
!
INTEGER             :: I, J, K, L(N), M(N)
REAL*8                :: BIGA, TEMP
REAL*8, PARAMETER     :: TOL = 1.0D-12


D = 1.0
DO K = 1, N
   L(K) = K
   M(K) = K
   BIGA = A(K,K)
   DO J = K, N
      DO I = K, N
         IF ( ABS(BIGA).LT.ABS(A(J,I)) ) THEN
            BIGA = A(J,I)
            L(K) = I
            M(K) = J
         END IF
      END DO
   END DO
   J = L(K)
   IF ( J.GT.K ) THEN
      DO I = 1, N
         TEMP = -A(I,K)
         A(I,K) = A(I,J)
         A(I,J) = TEMP
      END DO
   END IF
   I = M(K)
   IF ( I.GT.K ) THEN
      DO J = 1, N
         TEMP = -A(K,J)
         A(K,J) = A(I,J)
         A(I,J) = TEMP
      END DO
   END IF
   IF ( ABS(BIGA).LT.TOL ) THEN
      D = 0.0
      RETURN
   END IF
   DO I = 1, N
      IF ( I.NE.K ) A(K,I) = A(K,I)/(-BIGA)
   END DO
   DO I = 1, N
      DO J = 1, N
         IF ( I.NE.K ) THEN
            IF ( J.NE.K ) A(J,I) = A(K,I)*A(J,K) + A(J,I)
         END IF
      END DO
   END DO
   DO J = 1, N
      IF ( J.NE.K ) A(J,K) = A(J,K)/BIGA
   END DO
   D = MAX(-1.0D25,MIN(1.0D25,D))
   D = D*BIGA
   A(K,K) = 1.0/BIGA
END DO

K = N
!
!     write matrix A for output
!
DO
   K = K - 1
   IF ( K.LE.0 ) EXIT
   I = L(K)
   IF ( I.GT.K ) THEN
      DO J = 1, N
         TEMP = A(K,J)
         A(K,J) = -A(I,J)
         A(I,J) = TEMP
      END DO
   END IF
   J = M(K)
   IF ( J.GT.K ) THEN
      DO I = 1, N
         TEMP = A(I,K)
         A(I,K) = -A(I,J)
         A(I,J) = TEMP
      END DO
   END IF
END DO

return
end subroutine dmatinv
