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
!     subroutine thermocal: calculate partition function and with it 
!           thermodynamical properties
!
!     part of QMDFF
!

SUBROUTINE THERMOCAL(A,B,C,avmom, &
   &       LINEAR,ATOM,SYM,WT,VIBS,NVIBS,ESCF, &
   &       tempinput,sthr,et,g298,ts,zp)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
DIMENSION VIBS(*)
LOGICAL::LINEAR,ATOM,pr
!
!     THE FOLLOWING CONSTANTS ARE NOW DEFINED:
!          PI  = CIRCUMFERENCE TO DIAMETER OF A CIRCLE
!          R   = GAS CONSTANT IN CALORIES/MOLE
!          H   = PLANCK'S CONSTANT IN ERG-SECONDS
!          AK  = BOLTZMANN CONSTANT IN ERG/DEGREE
!          AC  = SPEED OF LIGHT IN CM/SEC
!
PI =3.14159D0 
R=1.98726D0
H=6.626176D-27
AK=1.3807D-16
AC=2.99792458D+10
siK=1.38066D-23
siH=6.626076D-34
siNA=6.022137D+23

pr=.false.

T=tempinput 
!
!   ***   INITIALISE SOME VARIABLES   ***
!
C1=H*AC/AK/T
QV=1.0D0
HV=0.0D0
E0=0.0D0
CPV=0.0D0
SV=0.0D0
ZP =0.0D0
!
!   ***   CONSTRUCT THE FREQUENCY DEPENDENT PARTS OF PARTITION FUNCTION
!
DO 50 I=1,NVIBS
   WI=VIBS(I)
   ZP=ZP+WI*1.4295718D-3
   EWJ=EXP(-WI*C1)
   QV=QV/(1-EWJ)
   HV=HV+WI*EWJ/(1.-EWJ)
   E0=E0+WI
   CPV=CPV+WI*WI*EWJ/(1.-EWJ)/(1.-EWJ)
!
!     replace low-lying vibs for S by rotor approx.
!              E in J (wi in cm-1)
!
   e=(1.d-14+wi)*0.01196266D+3/siNA
!
!     mom of intertia corresponding to the rotor with frequency wi         
!   
   xxmom=(siH/(2.0d0*PI))**2.0d0/(2.0d0*e)
!
!     write(*,*) xxmom,avmom,xxmom*avmom/(xxmom+avmom)
!     this reduced moment limits the rotational moment of
!     inertia for this vib to that of the total molecule
!     rotation/3
!
   xmom=xxmom*avmom/(xxmom+avmom)
!
!     free rotor entropy
!     Cramer, page 328 for one degree of freedom or
!     http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
!
   SV3=R*(0.5d0+LOG(SQRT(8.*PI**3.*XMOM*siK*T)/siH))
!
!     harm. osc. entropy
!
   if(wi.gt.0)then
      SV1=LOG(1.0D0-EWJ)
      SV2=WI*EWJ/(1.0d0-EWJ)
   else
      SV1=0
      SV2=0
   end if        
!
!     fermi weigthing
!     wofrot=1./(1.+exp( (wi-sthr)/20.0 ) )
!     Head-Gordon weighting
!
   wofrot=1-1./(1.+(sthr/(wi+1.d-14))**4.0)
   if (sthr.lt.0) wofrot=0.0d0
   SV=SV+(1.-wofrot)*(SV2*R*C1-R*SV1)+wofrot*SV3
   if (wi.lt.300.and.pr) & 
      &    write(10,'(''omega='',F6.2,5x,''Svib='',F10.3, &
      &       5x,''Srot='',F10.3,5x,''Sused='',F10.3)') &
      &       wi,(SV2*R*C1-R*SV1)*t/1000.,SV3*t/1000., &
      &       t*((1.-wofrot)*(SV2*R*C1-R*SV1)+wofrot*SV3)/1000.
50    CONTINUE
!
!   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
!
HV=HV*R*H*AC/AK
E0=E0*1.4295D0
CPV=CPV*R*C1*C1
!   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST LINEAR MOLECULES
IF (ATOM       ) GOTO 70
IF (.NOT.LINEAR) GOTO 60
QR=1/(C1*A*SYM)
HR=R*T
CPR=R
SR=R*(LOG(T*AK/(H*AC*A*SYM)))+R
GOTO 70
60  QR=SQRT(PI/(A*B*C*C1*C1*C1))/SYM
HR=3.0D0*R*T/2.0D0
CPR=3.0D0*R/2.0D0
SR=0.5D0*R*(3.D0*LOG(T*AK/(H*AC)) &
    &  -2.D0*LOG(SYM)+LOG(PI/(A*B*C))+3.D0)
70    CONTINUE
!
!   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
!
QINT=QV*QR
HINT=HV+HR
CPINT=CPV+CPR
SINT=SV+SR
!
!   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
!
QTR=(SQRT(2.D0*PI*WT*T*AK*1.6606D-24)/H)**3
!
!     this is 3/2RT+PV=5/2RT
!
HTR=5.0D0*R*T/2.0D0
CPTR=5.0D0*R/2.0D0
STR=2.2868D0*(5.0D0*LOG10(T)+3.0D0*LOG10(WT))-2.3135D0
!
!   ***   CONSTRUCT TOTALS   ***
!
CPTOT=CPTR+CPINT
STOT=STR+SINT
HTOT=HTR+HINT

autokcal=627.509541
autokj  =2625.49
caltoj  =autokj/autokcal

!if (pr)then
   write(10,*)
   WRITE(10,'('' Calculated thermodynamic properties: '')')
   write(10,*) "-------------------------------------------", &
           & "------------------------" 
   WRITE(10,'(''   T (K)   part. function  '', &
      & ''    enthalpy   heat capacity  entropy'')')
   WRITE(10,'(  ''                             '', &
      & ''   cal/mole    cal/K/mol   cal/K/mol'')')
   write(10,*) "-------------------------------------------", &
           & "------------------------"
   WRITE(10,'(f7.2,''  vib.'',G15.4 &
      & , 2X,3F11.5        )')T,QV,  HV,  CPV,  SV
   WRITE(10,'(7X,''  rot.'',G13.3 &
      &  ,2X,3F11.3        )')      QR,  HR,  CPR,  SR
   WRITE(10,'(7X,''  int.'',G13.3 &
      &  ,2X,3F11.3        )')      QINT,HINT,CPINT,SINT
   WRITE(10,'(7X,''  tra.'',G13.3 &
      &  , 2X,3F11.3)') QTR, HTR, CPTR, STR
   WRITE(10,'(7X,''  tot.'',13X,3x  ,F11.4,3F11.4)') &
      &   HTOT,CPTOT,STOT
   write(10,*) "-------------------------------------------", &
           & "------------------------"
!end if

ht=htot/1000.
et=ht+zp
ts=stot*t/1000.

g298=et-ts


!if (pr) then
   write(10,*)
   write(10,*) "-------------------------------------------", &
           & "----------------"
   WRITE(10,'(20x,''[au]       [kcal/mol]        [kJ/mol]'')')
   write(10,*) "-------------------------------------------", &
           & "----------------"
   WRITE(10,'('' ZPVE        '',3G16.6)')zp/autokcal,zp,zp*caltoj
   WRITE(10,'('' H(0)-H(T)+PV'',3G16.6)')ht/autokcal,ht,ht*caltoj
   write(10,*) "--------------------------------------------", &
       &     "----------------"
   WRITE(10,'('' H(T)        '',3G16.6)')et/autokcal,et,et*caltoj
   WRITE(10,'('' T*S         '',3G16.6)')ts/autokcal,ts,ts*caltoj
   WRITE(10,'('' G(T)        '',3G16.6)')g298/autokcal,g298,g298*caltoj
   write(10,*) "--------------------------------------------", &
       &     "----------------"

!end if

return
end subroutine thermocal
