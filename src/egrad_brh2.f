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
!   Potential energy routine taken from the POTLIB library:
!   R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, A. W. Jasper, D. G. Truhlar, T. C. Allison, 
!   A. F. Wagner, B. C. Garrett, J. Espinosa-García, and J. C. Corchado, POTLIB, 
!   http://comp.chem.umn.edu/potlib. 

      subroutine egrad_brh2(q,Natoms,Nbeads,V,dVdq,info)
      integer, intent(in) :: Natoms
      integer, intent(in) :: Nbeads
      double precision, intent(in) :: q(3,Natoms,Nbeads)
      double precision, intent(out) :: V(Nbeads)
      double precision, intent(out) :: dVdq(3,Natoms,Nbeads)

      double precision :: R(3), dVdr(3)
      double precision :: xAB, yAB, zAB, rAB
      double precision :: xAC, yAC, zAC, rAC
      double precision :: xBC, yBC, zBC, rBC
      integer k, info
        info = 0
      do k = 1, Nbeads

        xAB = q(1,2,k) - q(1,1,k)
        yAB = q(2,2,k) - q(2,1,k)
        zAB = q(3,2,k) - q(3,1,k)
        rAB = sqrt(xAB * xAB + yAB * yAB + zAB * zAB)
        R(1) = rAB

        xAC = q(1,1,k) - q(1,3,k)
        yAC = q(2,1,k) - q(2,3,k)
        zAC = q(3,1,k) - q(3,3,k)
        rAC = sqrt(xAC * xAC + yAC * yAC + zAC * zAC)
        R(2) = rAC

        xBC = q(1,3,k) - q(1,2,k)
        yBC = q(2,3,k) - q(2,2,k)
        zBC = q(3,3,k) - q(3,2,k)
        rBC = sqrt(xBC * xBC + yBC * yBC + zBC * zBC)
        R(3) = rBC
        call pot_brh2(R, V(k), dVdr)
        dVdq(1,1,k) = dVdr(2) * xAC / rAC - dVdr(1) * xAB / rAB
        dVdq(2,1,k) = dVdr(2) * yAC / rAC - dVdr(1) * yAB / rAB
        dVdq(3,1,k) = dVdr(2) * zAC / rAC - dVdr(1) * zAB / rAB

        dVdq(1,2,k) = dVdr(1) * xAB / rAB - dVdr(3) * xBC / rBC
        dVdq(2,2,k) = dVdr(1) * yAB / rAB - dVdr(3) * yBC / rBC
        dVdq(3,2,k) = dVdr(1) * zAB / rAB - dVdr(3) * zBC / rBC

        dVdq(1,3,k) = dVdr(3) * xBC / rBC - dVdr(2) * xAC / rAC
        dVdq(2,3,k) = dVdr(3) * yBC / rBC - dVdr(2) * yAC / rAC
        dVdq(3,3,k) = dVdr(3) * zBC / rBC - dVdr(2) * zAC / rAC
      end do
c      stop "in pooott"
      end subroutine egrad_brh2




      SUBROUTINE initialize_brh2
C
C   System:          BrH2
C   Functional form: Diatomics-in-molecules plus three-center term
C   Common name:     BrH2 DIM-3C
C   Interface:       3-1S
C   References:      D. C. Clary
C                    Chem. Phys. 71, 117 (1982)
C                    I. Last and M. Baer
C                    in Potential Energy Surfaces and Dynamics
C                    edited by D. G. Truhlar (Plenum, New York, 1981)
C                    p. 519
C
C   Number of bodies: 3
C   Number of derivatives: 1
C   Number of electronic surfaces: 1
C
C   Calling Sequence: 
C      PREPOT - initializes the potential's variables and
C               must be called once before any calls to POT
C      POT    - driver for the evaluation of the energy and the derivatives 
C               of the energy with respect to the coordinates for a given 
C               geometric configuration
C
C   Units: 
C      energies    - hartrees
C      coordinates - bohr
C      derivatives - hartrees/bohr
C
C   Surfaces: 
C      ground electronic state
C
C   Zero of energy: 
C      The classical potential energy is set equal to zero for the Br
C      infinitely far from the H2 diatomic and R(H2) set equal to the
C      H2 equilibrium diatomic value.
C
C   Parameters:
C      Set in data statements
C
C   Coordinates:
C      Internal, Definition: R(1) = R(Br-H)
C                            R(2) = R(H-H)
C                            R(3) = R(Br-H)
C
C   Common Blocks (used between the calling program and this potential):
C      /PT1CM/ R(3), ENERGY, DEDR(3)
C        passes the coordinates, ground state electronic energy, and 
C        derivatives of the ground electronic state energy with respect 
C        to the coordinates.
C      /PT4CM/ IPTPRT, IDUM(19)
C        passes the FORTRAN unit number used for output from the potential
C      /PT5CM/ EASYAB, EASYBC, EASYAC
C        passes the energy in the three asymptotic valleys for an A + BC system.
C        The energy in the AB valley, EASYAB, is equal to the energy of the 
C        C atom "infinitely" far from the AB diatomic and R(AB) set equal to 
C        Re(AB), the equilibrium bond length for the AB diatomic.  
C        In this potential the AB valley represents H infinitely far from
C        the BrH diatomic and R(BrH) equal to Re(BrH).  Similarly, the terms
C        EASYBC and EASYAC represent the energies in the H2 and the HBr 
C        valleys, respectively.
C
C   Default Parameter Values:
C      Variable      Default value
C      IPTPRT             6
C
C*****
C                    
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c         COMMON /PT1CM/ R(3), ENERGY, DEDR(3)
         COMMON /PT4CM_brh2/ IPTPRT, IDUM(19)
         COMMON /PT5CM_brh2/ EASYAB, EASYBC, EASYAC
C
C
      LOGICAL LCOL
      DIMENSION H(10),DH1(10),DH2(10),DH3(10),U(4,4),E(4),SCR1(4),
     1           SCR2(4),R(3),dedr(3)
      COMMON /POTCOM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,
     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P
      save
C
C   Echo potential parameters to FORTRAN unit IPTPRT
c      WRITE (IPTPRT, 600) RHH,AHH,DHH,BHH,RHX,AHX,DHX,BHX
c      WRITE (IPTPRT, 610) ETAHH,ETA3S,ETA1P,ETA3P,
c     *                    XIH,XIX,G,ALFW
C
      EASYAB = DHX
      EASYBC = DHH
      EASYAC = DHX
C
      XIB = 0.5D0*(XIH + XIX)
      T = XIB*RHX
      C3 = 1.D0/3.D0
      DENOM = 2.D0*RHX*(1.D0+T*(1.D0+C3*T))*EXP(-T)
      G = G/RHX
      RT34 = 0.25D0*SQRT(3.D0)
      RT38 = 0.5D0*RT34
C
600   FORMAT(/,1X,'*****','Potential Energy Surface',1X,'*****',
     *      //,1X,T5,'BrH2 DIM-3C potential energy surface',
     *       //,1X,T5,'Potential energy surface parameters:',
     *       /,1X,T5,'HH parameters (atomic units):',
     *       /,1X,T10,'RHH',T16,'=',T18,1PE13.5,
     *            T40,'AHH',T46,'=',T48,1PE13.5,
     *       /,1X,T10,'DHH',T16,'=',T18,1PE13.5,
     *            T40,'BHH',T46,'=',T48,1PE13.5,
     *       /,1X,T5,'HBr parameters (atomic units):',
     *       /,1X,T10,'RHBr',T16,'=',T18,1PE13.5,
     *            T40,'AHBr',T46,'=',T48,1PE13.5,
     *       /,1X,T10,'DHBr',T16,'=',T18,1PE13.5,
     *            T40,'BHBr',T46,'=',T48,1PE13.5)
610   FORMAT(/,1X,T5,'Other parameters (atomic units):',
     *       /,1X,T10,'ETAHH',T16,'=',T18,1PE13.5,
     *            T40,'ETA3S',T46,'=',T48,1PE13.5,
     *       /,1X,T10,'ETA1P',T16,'=',T18,1PE13.5,
     *            T40,'ETA3P',T46,'=',T48,1PE13.5,
     * /,1X,T10,'XIH',T16,'=',T18,1PE13.5,T40,'XIX',T46,'=',T48,1PE13.5,
     * /,1X,T10,'G',T16,'=',T18,1PE13.5,T40,'ALFW',T46,'=',T48,1PE13.5,
     * //,1X,'*****')
C
      RETURN
c      end 
C   Changed to new subroutine layout with dummy parameters for calling
      entry POT_brh2(R,energy,dedr)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      real(kind=8)::R(3),dedr(3),energy
c      COMMON /PT4CM/ IPTPRT, IDUM(19)
c      COMMON /PT5CM/ EASYAB, EASYBC, EASYAC

C
C   The potential routines define R1=R(Br-H), R2=R(H-Br), and R3=R(H-H);
C   rearrange the input coordinates to match these definitions
C
      R1 = R(1)
      R2 = R(3)
      R3 = R(2)
C
      IF(R1.GT.R2+R3) R1 = R2+R3
      IF(R2.GT.R1+R3) R2 = R1+R3
      IF(R3.GT.R1+R2) R3 = R1+R2
      R1S = R1*R1
      R2S = R2*R2
      R3S = R3*R3
      R12 = R1*R2
      R13 = R1*R3
      R23 = R2*R3
      T3 = (R1S+R2S-R3S)*0.5D0
      T2 = (R1S+R3S-R2S)*0.5D0
      T1 = (R2S+R3S-R1S)*0.5D0
      CSA = T2/R13
      IF(ABS(CSA) .GT. ONE) CSA = SIGN(ONE,CSA)
      CSA2 = CSA*CSA
      SNA2 = 1.D0-CSA2
      LCOL = SNA2 .LT. EPS
      SNA = SQRT(SNA2)
      T22 = 2.D0*CSA
      T11 = 0.D0
      IF(.NOT.LCOL) T11 = 2.D0*(SNA2-CSA2)/SNA
      SN2A = T22*SNA
      T = T3/(R1*R13)
      DCSA21 = T22*T
      DSN2A1 = T11*T
      T = -R2/R13
      DCSA22 = T22*T
      DSN2A2 = T11*T
      T = T1/(R13*R3)
      DCSA23 = T22*T
      DSN2A3 = T11*T
      CSB = T1/R23
      IF(ABS(CSB).GT.ONE) CSB = SIGN(ONE,CSB)
      CSB2 = CSB*CSB
      SNB2 = 1.D0-CSB2
      SNB = SQRT(SNB2)
      T11 = 0.D0
      IF(.NOT.LCOL) T11 = 2.D0*(SNB2-CSB2)/SNB
      T22 = 2.D0*CSB
      SN2B = T22*SNB
      T = -R1/R23
      DCSB21 = T22*T
      DSN2B1 = T11*T
      T = T3/(R2*R23)
      DCSB22 = T22*T
      DSN2B2 = T11*T
      T = T2/(R23*R3)
      DCSB23 = T22*T
      DSN2B3 = T11*T
      T = T3/R12
      CSG2 = T*T
      T = 2.D0*T
      DCSG21 = T*T2/(R1*R12)
      DCSG22 = T*T1/(R12*R2)
      DCSG23 = -T*R3/R12
C
C  DIATOMIC CURVES
C    HH
      RDIF = R3 - RHH
      EX1 = EXP(-AHH*RDIF)
      RDIF2 = RDIF*RDIF
      EX2 = EXP(-BHH*RDIF*RDIF2)
      T1 = DHH*EX1*EX2
      V1HH = T1*(EX1-2.D0)
      V3HH = ETAHH*T1*(EX1+2.D0)
      T1 = 2.D0*AHH*T1
      T2 = 3.D0*BHH*RDIF2
      DV1HH = T1*(1.D0-EX1) - T2*V1HH
      DV3HH = -T1*ETAHH*(1.D0+EX1) - T2*V3HH
C    HX
      CALL VHX(R1,V1S1,V3S1,V1P1,V3P1,DV1S1,DV3S1,DV1P1,DV3P1)
      CALL VHX(R2,V1S2,V3S2,V1P2,V3P2,DV1S2,DV3S2,DV1P2,DV3P2)
C
      S11 = V1S1 + 3.D0*V3S1
      S12 = V1S2 + 3.D0*V3S2
      S21 = 3.D0*V1S1 + V3S1
      S22 = 3.D0*V1S2 + V3S2
      S31 = V1S1 - V3S1
      S32 = V1S2 - V3S2
      P11 = V1P1 + 3.D0*V3P1
      P12 = V1P2 + 3.D0*V3P2
      P21 = 3.D0*V1P1 + V3P1
      P22 = 3.D0*V1P2 + V3P2
      P31 = V1P1 - V3P1
      P32 = V1P2 - V3P2
C
      DS11 = DV1S1 + 3.D0*DV3S1
      DS12 = DV1S2 + 3.D0*DV3S2
      DS21 = 3.D0*DV1S1 + DV3S1
      DS22 = 3.D0*DV1S2 + DV3S2
      DS31 = DV1S1 - DV3S1
      DS32 = DV1S2 - DV3S2
      DP11 = DV1P1 + 3.D0*DV3P1
      DP12 = DV1P2 + 3.D0*DV3P2
      DP21 = 3.D0*DV1P1 + DV3P1
      DP22 = 3.D0*DV1P2 + DV3P2
      DP31 = DV1P1 - DV3P1
      DP32 = DV1P2 - DV3P2
C
C    CONSTRUCT 4X4 HAMILTONIAN PUT IN PACKED ARRAY
C
      H(1) = V1HH + 0.25D0*(S11*CSA2 + P11*SNA2 + S12*CSB2 + P12*SNB2)
      H(3) = V1HH + 0.25D0*(S11*SNA2 + P11*CSA2 + S12*SNB2 + P12*CSB2)
      H(6) = V3HH + 0.25D0*(S21*CSA2 + P21*SNA2 + S22*CSB2 + P22*SNB2)
      H(10) = V3HH + 0.25D0*(S21*SNA2 + P21*CSA2 + S22*SNB2 + P22*CSB2)
      T11 = S11 - P11
      T12 = S12 - P12
      H(2) = 0.D0
      IF(.NOT.LCOL) H(2) = 0.125D0*(T11*SN2A - T12*SN2B)
      H(4) = RT34*(S31*CSA2 + P31*SNA2 - S32*CSB2 - P32*SNB2)
      T31 = S31 - P31
      T32 = S32 - P32
      H(5) = 0.D0
      IF(.NOT.LCOL) H(5) = RT38*(T31*SN2A + T32*SN2B)
      H(7) = H(5)
      H(8) = RT34*(S31*SNA2 + P31*CSA2 - S32*SNB2 - P32*CSB2)
      T21 = S21 - P21
      T22 = S22 - P22
      H(9) = 0.D0
      IF(.NOT.LCOL) H(9) = 0.125D0*(T21*SN2A - T22*SN2B)
C
C   CONSTRUCT DERIVATIVE OF HAMILTONIAN MATRIX
      T = T11*DCSA21 + T12*DCSB21
      DH1(1) = 0.25D0*(DS11*CSA2 + DP11*SNA2 + T)
      DH1(2) = 0.125D0*((DS11-DP11)*SN2A + T11*DSN2A1 - T12*DSN2B1)
      DH1(3) = 0.25D0*(DS11*SNA2 + DP11*CSA2 - T)
      T = T11*DCSA22 + T12*DCSB22
      DH2(1) = 0.25D0*(DS12*CSB2 + DP12*SNB2 + T)
      DH2(2) = 0.125D0*(-(DS12-DP12)*SN2B + T11*DSN2A2 - T12*DSN2B2)
      DH2(3) = 0.25D0*(DS12*SNB2 + DP12*CSB2 - T)
      T = 0.25D0*(T11*DCSA23 + T12*DCSB23)
      DH3(1) = DV1HH + T
      DH3(2) = 0.125D0*(T11*DSN2A3 - T12*DSN2B3)
      DH3(3) = DV1HH - T
      T = T21*DCSA21 + T22*DCSB21
      DH1(6) = 0.25D0*(DS21*CSA2 + DP21*SNA2 + T)
      DH1(9) = 0.125D0*((DS21-DP21)*SN2A + T21*DSN2A1 - T22*DSN2B1)
      DH1(10) = 0.25D0*(DS21*SNA2 + DP21*CSA2 - T)
      T = T21*DCSA22 + T22*DCSB22
      DH2(6) = 0.25D0*(DS22*CSB2 + DP22*SNB2 + T)
      DH2(9) = 0.125D0*((DP22-DS22)*SN2B + T21*DSN2A2 - T22*DSN2B2)
      DH2(10) = 0.25D0*(DS22*SNB2 + DP22*CSB2 - T)
      T = 0.25D0*(T21*DCSA23 + T22*DCSB23)
      DH3(6) = DV3HH + T
      DH3(9) = 0.125D0*(T21*DSN2A3 - T22*DSN2B3)
      DH3(10) = DV3HH - T
      T = T31*DCSA21 - T32*DCSB21
      DH1(4) = RT34*(DS31*CSA2 + DP31*SNA2 + T)
      DH1(5) = RT38*((DS31-DP31)*SN2A + T31*DSN2A1 + T32*DSN2B1)
      DH1(7) = DH1(5)
      DH1(8) = RT34*(DS31*SNA2 + DP31*CSA2 - T)
      T = T31*DCSA22 - T32*DCSB22
      DH2(4) = RT34*(-DS32*CSB2 - DP32*SNB2 + T)
      DH2(5) = RT38*((DS32-DP32)*SN2B + T31*DSN2A2 + T32*DSN2B2)
      DH2(7) = DH2(5)
      DH2(8) = RT34*(-DS32*SNB2 - DP32*CSB2 - T)
      T = T31*DCSA23 - T32*DCSB23
      DH3(4) = RT34*T
      DH3(5) = RT38*(T31*DSN2A3 + T32*DSN2B3)
      DH3(7) = DH3(5)
      DH3(8) = -DH3(4)
C
C    DIAGONALIZE H PUT EIGENVECTORS INTO U, EIGENVALUES INTO E
C    USE EISPACK ROUTINE RSP
C
      CALL RSP(4,4,10,H,E,1,U,SCR1,SCR2,IERR)
      IF(IERR.NE.0) THEN
         WRITE(IPTPRT,6000) R1,R2,R3,H
         WRITE(IPTPRT,6001) V1HH,V3HH,V1S1,V3S1,V1S2,
     *                  V3S2,V1P1,V3P1,V1P2,V3P2
C
         STOP 'POT 3'
      END IF
C    EVALUATE DERIVATIVES OF LOWEST ROOT
      D1 = 0.D0
      D2 = 0.D0
      D3 = 0.D0
      DO 20 L = 1,4
      T1 = U(L,1)
      LL = (L*(L-1))/2
      DO 20 K = 1,4
      T2 = T1*U(K,1)
      KK = (K*(K-1))/2
      IF(L.GE.K) INDEX = LL + K
      IF(K.GT.L) INDEX = KK + L
      D1 = D1 + T2*DH1(INDEX)
      D2 = D2 + T2*DH2(INDEX)
      D3 = D3 + T2*DH3(INDEX)
   20 CONTINUE
C
C   SUPPLEMENTARY TERM
      T = XIH*R3
      EX1 = EXP(-T)
      SHH = (1.D0+T*(1.D0+C3*T))*EX1
      DSHH = -XIH*T*(1.D0+T)*C3*EX1
      T = XIB*R1
      EX1 = EXP(-T)
      SHX1 = R1*(1.D0+T*(1.D0+C3*T))*EX1/DENOM
      DSHX1 = (1.D0+T*(1.D0-C3*T*T))*EX1/DENOM
      T = XIB*R2
      EX1 = EXP(-T)
      SHX2 = R2*(1.D0+T*(1.D0+C3*T))*EX1/DENOM
      DSHX2 = (1.D0+T*(1.D0-C3*T*T))*EX1/DENOM
      RDIF = R1 - R2
      T1 = G*EXP(-ALFW*RDIF*RDIF)
      T2 = SHH*(SHX1+SHX2) + SHX1*SHX2
      ENERGY = E(1) + T2*T1*CSG2 + DHH
      T3 = 2.D0*ALFW*RDIF*CSG2
      DEDR(1) = D1 + (DSHX1*(SHH+SHX2)*CSG2 + T2*(DCSG21-T3))*T1
      DEDR(3) = D2 + (DSHX2*(SHH+SHX1)*CSG2 + T2*(DCSG22+T3))*T1
      DEDR(2) = D3 + (DSHH*(SHX1+SHX2)*CSG2 + T2*DCSG23)*T1
C
6000  FORMAT(/,2X,'Error returned by the subprogram RSP in POT', 
     *       /,2X,'R = ',T10,3(1PE13.5,1X),
     *       /,2X,'H = ',(T10,4(1PE13.5,1X)))
6001  FORMAT(2X,'V1HH,V3HH,V1S1,V3S1,V1S2,V3S2,V1P1,V3P1,V1P2,V3P2=',
     *       /,2X,(T10,4(1PE13.5,1X)))
C
      RETURN
      END
      SUBROUTINE VHX(R,V1S,V3S,V1P,V3P,DV1S,DV3S,DV1P,DV3P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /POTCOM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,
     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P
      RDIF = R - RHX
      RDIF2 = RDIF*RDIF
      EX1 = EXP(-AHX*RDIF)
      EX2 = EXP(-BHX*RDIF2*RDIF)
      T1 = DHX*EX1*EX2
      V1S = T1*(EX1-2.D0)
      V1 = T1*(EX1+2.D0)
      V3S = ETA3S*V1
      V1P = ETA1P*V1
      V3P = ETA3P*V1
      T1 = 2.D0*AHX*T1
      T2 = 3.D0*BHX*RDIF2
      DV1S = T1*(1.D0-EX1) - T2*V1S
      DV1 = -T1*(1.D0+EX1) - T2*V1
      DV3S = ETA3S*DV1
      DV1P = ETA1P*DV1
      DV3P = ETA3P*DV1
      RETURN
      END
C
      BLOCK DATA PTPARM
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PT4CM_brh2/  IPTPRT, IDUM(19)
      COMMON /POTCOM/ EPS,ONE,RHH,AHH,DHH,BHH,XIH,XIX,G,ALFW,RHX,AHX,
     &                DHX,BHX,ETAHH,ETA3S,ETA1P,ETA3P
        DATA EPS /1.D-6/, ONE /1.0D0/
        DATA RHH,AHH,DHH,BHH /1.4016D0, 1.0291D0, 0.17447D0, 0.018D0/
        DATA XIH,XIX,G,ALFW /1.D0, 1.6D0, 0.22D0, 1.0D0/
        DATA IPTPRT /6/
        DATA RHX /2.673D0/
        DATA AHX /0.957D0/
        DATA DHX /0.1439D0/
        DATA BHX /0.012D0/
        DATA ETAHH /0.393764D0/
        DATA ETA3S /0.322D0/
        DATA ETA1P /0.20286D0/
        DATA ETA3P /0.1771D0/
      END
C
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2,IERR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(NV),W(N),Z(NM,N),FV1(N),FV2(N)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC PACKED MATRIX.
C
C     ON INPUT-
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT,
C
C        N  IS THE ORDER OF THE MATRIX  A,
C
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE
C        LESS THAN  N*(N+1)/2,
C
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C        PACKED MATRIX STORED ROW-WISE,
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT-
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DATA ZERO,ONE/0.0D0,1.0D0/
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50
C
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 DO 40 I = 1, N
C
         DO 30 J = 1, N
            Z(J,I) = ZERO
   30    CONTINUE
C
         Z(I,I) = ONE
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
C     ********** LAST CARD OF RSP **********
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION A(NV),D(N),E(N),E2(N)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DATA ZERO/0.0D0/
      DO  300 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I * L) / 2
         H = ZERO
         SCALE = ZERO
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
  120    CONTINUE
C
         IF (SCALE .NE. ZERO) GO TO 140
  130    E(I) = ZERO
         E2(I) = ZERO
         GO TO 290
C
  140    DO 150 K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 290
         F = ZERO
C
         DO 240 J = 1, L
            G = ZERO
            JK = (J * (J-1)) / 2
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, L
               JK = JK + 1
               IF (K .GT. J) JK = JK + K - 2
               G = G + A(JK) * D(K)
  180       CONTINUE
C     ********** FORM ELEMENT OF P **********
            E(J) = G / H
            F = F + E(J) * D(J)
  240    CONTINUE
C
         HH = F / (H + H)
         JK = 0
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
  260    CONTINUE
C
  290    D(I) = A(IZ+1)
         A(IZ+1) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED3 **********
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(N),E2(N)
      DOUBLE PRECISION MACHEP
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C                **********
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
      MACHEP = TWO**(-37)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F = ZERO
      B = ZERO
      E2(N) = ZERO
C
      DO 290 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + SQRT(E2(L)))
         IF (B .GT. H) GO TO 105
         B = H
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105    DO 110 M = L, N
            IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 210
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         S = SQRT(E2(L))
         G = D(L)
         P = (D(L1) - G) / (TWO * S)
         R = SQRT(P*P+ONE)
         D(L) = S / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G .EQ. ZERO) G = B
         H = G
         S = ZERO
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G .EQ. ZERO) G = B
            H = G * P / R
  200    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H .EQ. ZERO) GO TO 210
         IF (ABS(E2(L)) .LE. ABS(C/H)) GO TO 210
         E2(L) = H * E2(L)
         IF (E2(L) .NE. ZERO) GO TO 130
  210    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 270
            D(I) = D(I-1)
  230    CONTINUE
C
  250    I = 1
  270    D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(N),E(N),Z(NM,N)
      DOUBLE PRECISION MACHEP
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C                **********
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
      MACHEP = TWO**(-37)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = ZERO
      B = ZERO
      E(N) = ZERO
C
      DO 240 L = 1, N
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 110 M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (TWO * E(L))
         R = SQRT(P*P+ONE)
         D(L) = E(L) / (P + SIGN(R,P))
         H = G - D(L)
C
         DO 140 I = L1, N
  140    D(I) = D(I) - H
C
         F = F + H
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C = ONE
         S = ZERO
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 200 II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+ONE)
            E(I+1) = S * P * R
            S = C / R
            C = ONE / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+ONE)
            E(I+1) = S * E(I) * R
            S = ONE / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NV),Z(NM,M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      DATA ZERO/0.0D0/
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H .EQ. ZERO) GO TO 140
C
         DO 130 J = 1, M
            S = ZERO
            IK = IZ
C
            DO 110 K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
  110       CONTINUE
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
            S = (S / H) / H
            IK = IZ
C
            DO 120 K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
C     ********** LAST CARD OF TRBAK3 **********
      END
