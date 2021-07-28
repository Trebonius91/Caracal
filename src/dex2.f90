!
!     function dex2: formate of basis function array in setsto3/4
!     for EHT calculations
!
!     part of QMDFF
!
double precision function dex2(m)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
IF(M .LT. 2) THEN
   DEX2=1
ELSE
   DEX2=1
   do I=1,M,2
      DEX2=DEX2*I
   end do
ENDIF
RETURN
end function dex2
