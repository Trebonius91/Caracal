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
!     fermismear: calculate fermi smearing: floating point 
!     occupations of border orbitals for EHT calculations
!
!     part of QMDFF
!
subroutine fermismear(prt,NORBS,NEL,T,eig,occ)
IMPLICIT NONE
integer norbs
integer nel
real*4  eig(norbs)
real*4  occ(norbs)
real*8  t
LOGICAL PRT

real*8 boltz,bkt,occt,e_fermi,total_number,thr
real*8 total_dfermi,dfermifunct,fermifunct,s,change_fermi

PARAMETER (BOLTZ = 3.166808578545117E-06*27.2112)
integer ncycle,i,j,m,k,i1,i2
!
!     S. Grimme, University Bonn, Dec. 2012      
!
BKT = BOLTZ*T

E_FERMI = 0.5*(EIG(NEL)+EIG(NEL+1))
OCCT=NEL

NCYCLE = 0
!
!     older form of do-loop...
!
   10 TOTAL_NUMBER = 0.0
   TOTAL_DFERMI = 0.0
   NCYCLE = NCYCLE+1
   DO I = 1, NORBS
      FERMIFUNCT = 0.0
      if((EIG(I)-E_FERMI)/BKT.lt.50) then
         FERMIFUNCT = 1.0/(EXP((EIG(I)-E_FERMI)/BKT)+1.0)
         DFERMIFUNCT = EXP((EIG(I)-E_FERMI)/BKT) / &
          & (BKT*(EXP((EIG(I)-E_FERMI)/BKT)+1.0)**2)
      ELSE
         DFERMIFUNCT = 0.0
      END IF
      OCC(I) = FERMIFUNCT
      TOTAL_NUMBER = TOTAL_NUMBER + FERMIFUNCT
      TOTAL_DFERMI = TOTAL_DFERMI + DFERMIFUNCT
   END DO
   CHANGE_FERMI = (OCCT-TOTAL_NUMBER)/TOTAL_DFERMI
   E_FERMI = E_FERMI+CHANGE_FERMI
IF(ABS(OCCT-TOTAL_NUMBER).GT.1.0E-8.AND.NCYCLE.LT.200) GOTO 10

IF (PRT) THEN
WRITE(*,'(''temp, E(Fermi) in cycles '',2F10.4,i5)') &
     & T,E_FERMI,NCYCLE
ENDIF

RETURN
end subroutine fermismear
