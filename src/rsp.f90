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
!     subroutine rsp: diagonalize symmetric matrix for eigenvalues 
!     and eigenvectors
!
!     part of QMDFF
!
      SUBROUTINE RSP(A,N,MATZ,W,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(n*(n+1)/2),  W(n), Z(n,n)  
      DIMENSION :: FV1(2*n),FV2(2*n)
!******************************************************************
!
!   EISPACK DIAGONALIZATION ROUTINES: TO FIND THE EIGENVALUES AND
!           EIGENVECTORS (IF DESIRED) OF A REAL SYMMETRIC PACKED MATRIX.
! ON INPUT-      N  IS THE ORDER OF THE MATRIX  A,
!                A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
!                   PACKED MATRIX STORED ROW-WISE,
!             MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF ONLY
!                   EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
!                   ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND
!                   EIGENVECTORS.
! ON OUTPUT-     W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
!                Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
!
!******************************************************************
! THIS SUBROUTINE WAS CHOSEN AS BEING THE MOST RELIABLE. (JJPS)
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      SAVE FIRST, EPS, ETA, NV
      LOGICAL FIRST

      if(n.eq.1)then
         Z(1,1)=1.0d0 
         W(1)=A(1)
         return
      endif

      CALL EPSETA(EPS,ETA)
      
      NV=(N*(N+1))/2
      NM=N 
      CALL  TRED3(N,NV,A,W,FV1,FV2,EPS,EPS)
      IF (MATZ .NE. 0) GO TO 10
!     ********** FIND EIGENVALUES ONLY **********
      CALL  TQLRAT_a(N,W,FV2,IERR,EPS)
      GO TO 40
!     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   10 DO 30    I = 1, N
         DO 20    J = 1, N
            Z(J,I)=0.0D0
   20    CONTINUE
         Z(I,I)=1.0D0
   30 CONTINUE
      CALL  TQL2_a(NM,N,W,FV1,Z,IERR,EPS)
      IF (IERR .NE. 0) GO TO 40
      CALL  TRBAK3(NM,N,NV,A,N,Z)
!     ********** LAST CARD OF RSP **********
   40 RETURN
      END
      SUBROUTINE EPSETA(EPS,ETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
!     COMPUTE AND RETURN ETA, THE SMALLEST REPRESENTABLE NUMBER,
!     AND EPS IS THE SMALLEST NUMBER FOR WHICH 1+EPS.NE.1.
!
      ETA = 1.D0
   10 IF((ETA/2.D0).EQ.0.D0) GOTO 20
      IF(ETA.LT.1.D-38) GOTO 20
      ETA = ETA / 2.D0
      GOTO 10
   20 EPS = 1.D0
   30 IF((1.D0+(EPS/2.D0)).EQ.1.D0) GOTO 40
      IF(EPS.LT.1.D-17) GOTO 40
      EPS = EPS / 2.D0
      GOTO 30
   40 RETURN
      END
      SUBROUTINE TQL2_a(NM,N,D,E,Z,IERR,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!               ===== PROCESSED BY AUGMENT, VERSION 4N =====
!     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
!               ----- LOCAL VARIABLES -----
!               ----- GLOBAL VARIABLES -----
      DIMENSION D(*), E(*), Z(NM,*)
!               ----- SUPPORTING PACKAGE FUNCTIONS -----
!               ===== TRANSLATED PROGRAM =====
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2_a,
!     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
!     WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
!     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
!     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
!     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
!     FULL MATRIX TO TRIDIAGONAL FORM.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
!
!        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
!          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
!          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
!          THE IDENTITY MATRIX.
!
!      ON OUTPUT-
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
!          UNORDERED FOR INDICES 1,2,...,IERR-1,
!
!        E HAS BEEN DESTROYED,
!
!        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
!          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
!          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
!          EIGENVALUES,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!
      IERR = 0
      IF (N .EQ. 1) GO TO 160
      DO 10   I = 2, N
      E(I-1) = E(I)
   10   CONTINUE
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
      DO 110   L = 1, N
         J = 0
         H=EPS*(ABS (D(L))+ABS (E(L)))
         IF (B .LT. H) B=H
!     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 20   M = L, N
            IF (ABS (E(M)).LE.B)  GO TO 30
!     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP **********
   20    CONTINUE
   30    IF (M .EQ. L) GO TO 100
   40    IF (J .EQ. 30) GO TO 150
         J = J + 1
!     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P=(D(L1)-G)/(2.0D0*E(L))
         R=SQRT (P*P+1.0D0)
         D(L)=E(L)/(P+SIGN (R,P))
         H = G - D(L)
         DO 50   I = L1, N
         D(I) = D(I) - H
   50    CONTINUE
         F = F + H
!     ********** QL TRANSFORMATION **********
         P = D(M)
         C=1.0D0
         S=0.0D0
         MML = M - L
!     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 90   II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS (P).LT.ABS (E(I)))  GO TO 60
            C = E(I) / P
            R=SQRT (C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C=1.0D0/R
            GO TO 70
   60       C = P / E(I)
            R=SQRT (C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S=1.0D0/R
            C = C * S
   70       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     ********** FORM VECTOR **********
            DO 80   K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
   80       CONTINUE
   90    CONTINUE
         E(L) = S * P
         D(L) = C * P
         IF (ABS (E(L)).GT.B)  GO TO 40
  100    D(L) = D(L) + F
  110 CONTINUE
!     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 140   II = 2, N
         I = II - 1
         K = I
         P = D(I)
         DO 120   J = II, N
            IF (D(J) .GE. P) GO TO 120
            K = J
            P = D(J)
  120    CONTINUE
         IF (K .EQ. I) GO TO 140
         D(K) = D(I)
         D(I) = P
         DO 130   J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  130    CONTINUE
  140 CONTINUE
      GO TO 160
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
  150 IERR = L
  160 RETURN
!     ********** LAST CARD OF TQL2_a **********
      END
      SUBROUTINE TQLRAT_a(N,D,E2,IERR,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!               ===== PROCESSED BY AUGMENT, VERSION 4N =====
!     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
!               ----- LOCAL VARIABLES -----
!               ----- GLOBAL VARIABLES -----
      DIMENSION D(*), E2(*)
!               ----- SUPPORTING PACKAGE FUNCTIONS -----
!               ===== TRANSLATED PROGRAM =====
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT_a,
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT-
!
!        N IS THE ORDER OF THE MATRIX,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
!
!        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT-
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES,
!
!        E2 HAS BEEN DESTROYED,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!
      IERR = 0
      IF (N .EQ. 1) GO TO 140
      DO 10   I = 2, N
      E2(I-1) = E2(I)
   10 CONTINUE
      F=0.0D0
      B=0.0D0
      E2(N)=0.0D0
      DO 120   L = 1, N
         J = 0
         H=EPS*(ABS (D(L))+SQRT (E2(L)))
         IF (B .GT. H) GO TO 20
         B = H
         C = B * B
!     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
   20    DO 30   M = L, N
            IF (E2(M) .LE. C) GO TO 40
!     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP **********
   30    CONTINUE
   40    IF (M .EQ. L) GO TO 80
   50    IF (J .EQ. 30) GO TO 130
         J = J + 1
!     ********** FORM SHIFT **********
         L1 = L + 1
         S=SQRT (E2(L))
         G = D(L)
         P=(D(L1)-G)/(2.0D0*S)
         R=SQRT (P*P+1.0D0)
         D(L)=S/(P+SIGN (R,P))
         H = G - D(L)
         DO 60   I = L1, N
         D(I) = D(I) - H
   60    CONTINUE
         F = F + H
!     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G.EQ.0.0D0) G=B
         H = G
         S=0.0D0
         MML = M - L
!     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 70   II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G.EQ.0.0D0) G=B
            H = G * P / R
   70    CONTINUE
         E2(L) = S * G
         D(L) = H
!     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H.EQ.0.0D0)  GO TO 80
         IF (ABS (E2(L)).LE.ABS (C/H))  GO TO 80
         E2(L) = H * E2(L)
         IF (E2(L).NE.0.0D0)  GO TO 50
   80    P = D(L) + F
!     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 100
!     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 90   II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 110
            D(I) = D(I-1)
   90    CONTINUE
  100    I = 1
  110    D(I) = P
  120 CONTINUE
      GO TO 140
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
  130 IERR = L
  140 RETURN
!     ********** LAST CARD OF TQLRAT_a **********
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!               ===== PROCESSED BY AUGMENT, VERSION 4N =====
!     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
!               ----- LOCAL VARIABLES -----
!               ----- GLOBAL VARIABLES -----
      DIMENSION A(*), Z(NM,*)
!               ===== TRANSLATED PROGRAM =====
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
!
!     ON INPUT-
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT,
!
!        N IS THE ORDER OF THE MATRIX,
!
!        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
!          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
!          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
!          N*(N+1)/2 POSITIONS,
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
!
!        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!          IN ITS FIRST M COLUMNS.
!
!     ON OUTPUT-
!
!        Z CONTAINS THE TRANSFORMED EIGENVECTORS
!          IN ITS FIRST M COLUMNS.
!
!     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
      IF (M .EQ. 0) GO TO 50
      IF (N .EQ. 1) GO TO 50
      DO 40   I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H.EQ.0.0D0)  GO TO 40
         DO 30   J = 1, M
            S=0.0D0
            IK = IZ
            DO 10   K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
   10       CONTINUE
!     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
            S = (S / H) / H
            IK = IZ
            DO 20   K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
   50 RETURN
!     ********** LAST CARD OF TRBAK3 **********
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2,EPS,ETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!               ===== PROCESSED BY AUGMENT, VERSION 4N =====
!     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
!               ----- LOCAL VARIABLES -----
!               ----- GLOBAL VARIABLES -----
      DIMENSION A(*), D(*), E(*), E2(*)
!               ----- SUPPORTING PACKAGE FUNCTIONS -----
!               ===== TRANSLATED PROGRAM =====
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
!     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
!     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT-
!
!        N IS THE ORDER OF THE MATRIX,
!
!        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
!          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
!
!        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
!          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
!          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
!
!     ON OUTPUT-
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
!          TRANSFORMATIONS USED IN THE REDUCTION,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 100   II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = ( I * L ) / 2
         H=0.0D0
         SCALE=0.0D0
         DO 10   K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE=SCALE+ABS( D(K) )
   10    CONTINUE
         IF ( SCALE.NE.0.D0 ) GO TO 20
         E(I)=0.0D0
         E2(I)=0.0D0
         GO TO 90
   20    DO 30   K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
   30    CONTINUE
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G=-SIGN (SQRT (H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 90
         F=0.0D0
         DO 70   J = 1, L
            G=0.0D0
            JK = (J * (J-1)) / 2
!     ********** FORM ELEMENT OF A*U **********
            K = 0
   40       K = K + 1
            JK = JK + 1
            G = G + A(JK) * D(K)
            IF ( K .LT. J ) GO TO 40
            IF ( K .EQ. L ) GO TO 60
   50       JK = JK + K
            K = K + 1
            G = G + A(JK) * D(K)
            IF ( K .LT. L ) GO TO 50
!     ********** FORM ELEMENT OF P **********
   60       CONTINUE
            E(J) = G / H
            F = F + E(J) * D(J)
   70    CONTINUE
!
         HH = F / (H + H)
         JK = 0
!     ********** FORM REDUCED A **********
         DO 80   J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
            DO 85   K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
   85    CONTINUE
   80    CONTINUE
         D(I) = A(IZ+1)
   90    CONTINUE
         A(IZ+1)=SCALE*SQRT (H)
  100 CONTINUE
!
      RETURN
!     ********** LAST CARD OF TRED3 **********
      END
