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

!
!     setsto4: STO-3G parameters in a four gaussian expansion for 
!      EHT basis set
!
!     part of QMDFF
!
subroutine setsto4(nprim,n,l,zeta,expo,cont)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
dimension expo(*),cont(*)
dimension ALLC(4,5,3),ALLZ(4,5,3)
integer::iam,nprim
real(kind=8)::DEX(-1:96)

nprim=4

PI=4.D0*ATAN(1.D0)
DO I=-1,10
   DEX(I)=DEX2(I)
ENDDO   
!     1s
ALLZ(1,1,1) =5.216844534
ALLZ(2,1,1) =9.546182760D-1
ALLZ(3,1,1) =2.652034102D-1
ALLZ(4,1,1) =8.801862774D-2
ALLC(1,1,1) =5.675242080D-2
ALLC(2,1,1) =2.601413550D-1
ALLC(3,1,1) =5.328461143D-1
ALLC(4,1,1) =2.916254405D-1
!     2s
ALLZ(1,2,1) =1.161525551D1
ALLZ(2,2,1) =2.000243111
ALLZ(3,2,1) =1.607280687D-1
ALLZ(4,2,1) =6.125744532D-2
ALLC(1,2,1) =-1.198411747D-2
ALLC(2,2,1) =-5.472052539D-2
ALLC(3,2,1) =5.805587176D-1
ALLC(4,2,1) =4.770079976D-1
!     2p
ALLZ(1,2,2) =1.798260992
ALLZ(2,2,2) =4.662622228D-1
ALLZ(3,2,2) =1.643718620D-1
ALLZ(4,2,2) =6.543927065D-2
ALLC(1,2,2) =5.713170255D-2
ALLC(2,2,2) =2.857455515D-1
ALLC(3,2,2) =5.517873105D-1
ALLC(4,2,2) =2.632314924D-1
!     3s
ALLZ(1,3,1) =1.513265591
ALLZ(2,3,1) =4.262497508D-1
ALLZ(3,3,1) =7.643320863D-2
ALLZ(4,3,1) =3.760545063D-2
ALLC(1,3,1) =-3.295496352D-2
ALLC(2,3,1) =-1.724516959D-1
ALLC(3,3,1) =7.518511194D-1
ALLC(4,3,1) =3.589627317D-1
!     3p
ALLZ(1,3,2) =1.853180239
ALLZ(2,3,2) =1.915075719D-1
ALLZ(3,3,2) =8.655487938D-2
ALLZ(4,3,2) =4.184253862D-2
ALLC(1,3,2) =-1.434249391D-2
ALLC(2,3,2) =2.755177580D-1
ALLC(3,3,2) =5.846750879D-1
ALLC(4,3,2) =2.144986514D-1
!     3p
ALLZ(1,4,1) =3.242212833D-1
ALLZ(2,4,1) =1.663217177D-1
ALLZ(3,4,1) =5.081097451D-2
ALLZ(4,4,1) =2.829066600D-2
ALLC(1,4,1) =-1.120682822D-1
ALLC(2,4,1) =-2.845426863D-1
ALLC(3,4,1) =8.909873788D-1
ALLC(4,4,1) =3.517811205D-1
!     4p
ALLZ(1,4,2) =1.492607880
ALLZ(2,4,2) =4.327619272D-1
ALLZ(3,4,2) =7.553156064D-2
ALLZ(4,4,2) =3.706272183D-2
ALLC(1,4,2) =-6.035216774D-3
ALLC(2,4,2) =-6.013310874D-2
ALLC(3,4,2) =6.451518200D-1
ALLC(4,4,2) =4.117923820D-1
!     5s
ALLZ(1,5,1) =8.602284252D-1
ALLZ(2,5,1) =1.189050200D-1
ALLZ(3,5,1) =3.446076176D-2
ALLZ(4,5,1) =1.974798796D-2
ALLC(1,5,1) =1.103657561D-2
ALLC(2,5,1) =-5.606519023D-1
ALLC(3,5,1) =1.179429987
ALLC(4,5,1) =1.734974376D-1
!     5p
ALLZ(1,5,2) =3.962838833D-1
ALLZ(2,5,2) =1.838858552D-1
ALLZ(3,5,2) =4.943555157D-2
ALLZ(4,5,2) =2.750222273D-2
ALLC(1,5,2) =-1.801459207D-2
ALLC(2,5,2) =-1.360777372D-1
ALLC(3,5,2) =7.533973719D-1
ALLC(4,5,2) =3.409304859D-1
!     3d
ALLZ(1,3,3) =9.185846715D-1
ALLZ(2,3,3) =2.920461109D-1
ALLZ(3,3,3) =1.187568890D-1
ALLZ(4,3,3) =5.286755896D-2
ALLC(1,3,3) =5.799057705D-2
ALLC(2,3,3) =3.045581349D-1
ALLC(3,3,3) =5.601358038D-1
ALLC(4,3,3) =2.432423313D-1
!     4d
ALLZ(1,4,3) =1.995825422
ALLZ(2,4,3) =1.823461280D-1
ALLZ(3,4,3) =8.197240896D-2
ALLZ(4,4,3) =4.000634951D-2
ALLC(1,4,3) =-2.816702620D-3
ALLC(2,4,3) =2.177095871D-1
ALLC(3,4,3) =6.058047348D-1
ALLC(4,4,3) =2.717811257D-1
!     5d
ALLZ(1,5,3) =4.230617826D-1
ALLZ(2,5,3) =8.293863702D-2
ALLZ(3,5,3) =4.590326388D-2
ALLZ(4,5,3) =2.628744797D-2
ALLC(1,5,3) =-2.421626009D-2
ALLC(2,5,3) =3.937644956D-1
ALLC(3,5,3) =5.489520286D-1
ALLC(4,5,3) =1.190436963D-1

iam=l-1
DO J=1,nprim
   cont(J)=ALLC(J,N,L)
   expo(J)=ALLZ(J,N,L)*zeta**2
   XNORM=(2.D0*EXPO(J)/PI)**0.75D0*(4.D0*EXPO(J))**(IAM/2.D0)/&
     &  SQRT(DEX(2*IAM-1))
   cont(j)=cont(j)*xnorm                
END DO      

end subroutine setsto4
