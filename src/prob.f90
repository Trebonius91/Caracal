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
!     Subroutine prob: Calculate probability distribution (histogram)
!     from umbrella samplings for wham subroutine
!
!     part of EVB
!

!*********************SUBROUTINE PROB BEGINS HERE************************************
!====================================================================================
!     FOR A GIVEN WINDOW, THIS SUBROUTINE RETURNS:
!     (a) THE NUMBER OF SAMPLED SNAPSHOTS (NI)
!     (b) THE NORMALIZED PROBABILITY DISTRIBUTION (HIST)
!     (c) THE BIASING POTENTIAL AROUND THE WINDOW CENTER (V)
!
!
subroutine prob(XI_L,DXI,nbins,XI_I,L,DT,TCUT,KS,NSAMP,F,NI,HIST,v_res)
use general
use evb_mod
implicit none
integer::I,J,nbins,NC,NAT,BIN,NSAMP,T_TOT,NMAX,NI
parameter (NMAX=5)
real(kind=8)::L,XI_L,XI_U,DXI,DT,TCUT,FAC,XI_I,KS
real(kind=8)::POS(3,NMAX),TMP,DIST
real(kind=8)::HIST(nbins),v_res(nbins),TMP2
CHARACTER(len=3)::F
data FAC/1.889725988578923D0/
!
!     OPEN FILES 1.7.xyz, 1.8.xyz,........,4.5.xyz (based on F)
!
!OPEN(10,FILE=F//".xyz",STATUS="OLD")
!
!     THE OUTPUT NORMALIZED PROBABILITY DISTRIBUTION WILL STORED IN
!     THE FILES PROB_1.7.dat, PROB_1.8.dat, AND SO ON.
!
!OPEN(11,FILE="PROB_"//F//".dat",STATUS="UNKNOWN")
!
!     INITIAL HISTOGRAM VALUES TO ZERO
!
DO I = 1, nbins
   HIST(I) = 0.D0
END DO
NC = 0 ! INITIALIZE TOTAL NUMBER OF SNAPSHOTS TO ZERO
NI = 0 ! INITIALIZE TOTAL NUMBER OF SAMPLED SNAPSHOTS TO ZERO
!
!     READING SNAPSHOTS UNTIL END OF FILE IS REACHED
!
DO WHILE(.TRUE.)
   READ(10,*,ERR=30,END=40)NAT
   READ(10,*,ERR=30,END=40)
   DO I = 1, NAT
      READ(10,*,ERR=30,END=40)F,(POS(J,I),J=1,3),TMP,TMP,TMP
   END DO
   NC = NC + 1
   T_TOT=DT*DBLE(NC) ! TOTAL SIMULATION TIME UP TILL NOW
   IF(MOD(NC,NSAMP).NE.0)GOTO 50 ! SAMPLING
   IF(T_TOT.LT.TCUT)GOTO 50 ! DISCARD PRE-EQULIBRATION SNAPSHOTS
   NI = NI + 1
   DIST = 0.D0 ! INITIALIZE C-Cl DISTANCE TO ZERO
   DO J = 1, 3 ! COMPUTE DIST (IN ANGSTROM) SUBJECT TO PBC
      TMP = POS(J,1)-POS(J,2)
      TMP = TMP - L*DNINT(TMP/L)
      DIST = DIST + TMP*TMP
   END DO
   DIST=DSQRT(DIST)
   DIST=DIST-XI_L ! SHIFT DISTANCE RELATIVE LOWER BOUND
   BIN = 1 + DINT(DIST/DXI) ! LOCATE BIN
   IF(BIN.GT.nbins)GOTO 50 ! DISCARD POINTS GREATER THAN XI_L
   HIST(BIN)=HIST(BIN) + 1.D0 ! PLACE IN APPROPRIATE BIN
50 CONTINUE
END DO
30 CONTINUE
WRITE(*,*) "ERROR IN READING INPUT FILE"

40 CONTINUE
CLOSE(10)
!
!     COMPUTE NORMALIZATION FACTOR
!
TMP = 0
DO I = 1, N
   TMP = TMP + HIST(I)
END DO
!
!     NOTE THAT, IN PRINCIPLE, TMP=DBLE(NI)
!     COMPUTE NORMALIZED DISTRIBUTION AND BIASING POTENTIAL
!
WRITE(11,'(A)') "# COORDINATE POTENTIAL PROBABILITY"
DO I = 1, nbins
   TMP2 = XI_L + DXI*(DBLE(I)-0.5D0) ! REACTION COORDINATE IN ANGSTROM
   HIST(I) = HIST(I)/TMP             ! NORMALIZED PROBALITY AT TMP2
   v_res(I) = KS*(TMP2-XI_I)**2          ! BIASING WINDOW POTENTIAL
   v_res(I) = v_res(I)*FAC**2                ! CONVERT POTENTIAL TO A.U.
   WRITE(11,'(3F12.6)')TMP2,HIST(I),v_res(I)
END DO
CLOSE(11)

RETURN
end subroutine prob

