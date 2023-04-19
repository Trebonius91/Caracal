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
!    subroutine cm5mod: original subroutine for calculation of CM5 charges
!
!    part of QMDFF
!
subroutine cm5mod(NAT,IZ,CHIR,Q,CCM5,DHIRX,DHIRY,DHIRZ, &
  & DCM5X,DCM5Y,DCM5Z)
IMPLICIT REAL*8 (A-H,O-Z)
PARAMETER (MZ=118)
DIMENSION Q(3,*),IZ(*),CHIR(*),CCM5(*)
DIMENSION RAD(MZ),A0(MZ),D(MZ,MZ)


!
!     COVALENT RADII of all elements
!     based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!     in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!     edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!     corrected Nov. 17, 2010 for the 92nd edition.
!
rad = (/ &
 & 0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0, &
 & 0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0, &
 & 2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0, &
 & 1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0, &
 & 2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0, &
 & 1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0, &
 & 2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0, &
 & 1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0, &
 & 1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0, &
 & 1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0, &
 & 1.84D0,1.83D0,1.80D0,1.80D0,1.73D0,1.68D0,1.68D0,1.68D0,1.65D0, &
 & 1.67D0,1.73D0,1.76D0,1.61D0,1.57D0,1.49D0,1.43D0,1.41D0,1.34D0, & 
 & 1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,1.75D0,1.65D0, &
 & 1.57D0 /)


!
!     CM5 MODEL PARAMETERS
!
DO I=1,MZ
   A0(I)=0.D0
   DO J=1,MZ
      D(I,J)=0.D0
   END DO
END DO

!
!     definition of atomwise parameters
!

a0= (/ &
 &  0.0056,-0.1543, 0.0000, 0.0333,-0.1030,-0.0446,-0.1072,-0.0802, &
 & -0.0629,-0.1088, 0.0184, 0.0000,-0.0726,-0.0790,-0.0756,-0.0565, &
 & -0.0444,-0.0767, 0.0130, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,-0.0512,-0.0557, &
 & -0.0533,-0.0399,-0.0313,-0.0541, 0.0092, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 & -0.0361,-0.0393,-0.0376,-0.0281,-0.0220,-0.0381, 0.0065, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 & -0.0255,-0.0277,-0.0265,-0.0198,-0.0155,-0.0269, 0.0046, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 &  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, &
 & -0.0179,-0.0195,-0.0187,-0.0140,-0.0110,-0.0189 /)

DO K1=1,MZ
   DO K2=K1+1,MZ
      D(K1,K2)=A0(K1)-A0(K2)
   END DO
END DO
!
!     PAIRWISE PARAMETERS
!
D( 1, 6)= 0.0502
D( 1, 7)= 0.1747
D( 1, 8)= 0.1671
D( 6, 7)= 0.0556
D( 6, 8)= 0.0234
D( 7, 8)=-0.0346

DO I=1,MZ
   DO J=I+1,MZ
      D(J,I)=-D(I,J)
   END DO
END DO
!     ALPHA
ALP=2.4740

!     C-COEFFICIENT: 0.7050   ! ALREADY INCLUDED IN A0


DO K=1,NAT
   CCM5(K)=CHIR(K)
   DO K1=1,NAT
      IF (IZ(K).NE.IZ(K1)) THEN
         DIS=DSQRT((Q(1,K)-Q(1,K1))**2+(Q(2,K)-Q(2,K1))**2+ &
           & (Q(3,K)-Q(3,K1))**2)
         BKK=DEXP(-ALP*(DIS-RAD(IZ(K))-RAD(IZ(K1))))
         CCM5(K)=CCM5(K)+BKK*D(IZ(K),IZ(K1))
      END IF
   END DO
END DO
DHIRX=0.D0
DHIRY=0.D0
DHIRZ=0.D0
DCM5X=0.D0
DCM5Y=0.D0
DCM5Z=0.D0
DO J=1,NAT
   DHIRX=DHIRX+Q(1,J)*CHIR(J)*4.803242D0
   DHIRY=DHIRY+Q(2,J)*CHIR(J)*4.803242D0
   DHIRZ=DHIRZ+Q(3,J)*CHIR(J)*4.803242D0
   DCM5X=DCM5X+Q(1,J)*CCM5(J)*4.803242D0
   DCM5Y=DCM5Y+Q(2,J)*CCM5(J)*4.803242D0
   DCM5Z=DCM5Z+Q(3,J)*CCM5(J)*4.803242D0
END DO

RETURN
end subroutine

