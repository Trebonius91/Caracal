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
!     setsto3: Define GTO expansion values for 3 gaussians for EHT basis
!
!     part of QMDFF
!
subroutine setsto3(nprim,n,l,zeta,expo,cont)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
dimension expo(*),cont(*)
dimension ALLC(3,5,3),ALLZ(3,5,3)
integer::iam,nprim
real(kind=8)::DEX(-1:96)

nprim=3

PI=4.D0*ATAN(1.D0)
DO I=-1,10
   DEX(I)=DEX2(I)
ENDDO   
!     1s: first three: alpha, second three d
ALLZ(1,1,1) =2.227660584D00 
ALLZ(2,1,1) =4.057711562D-01
ALLZ(2,1,1) =4.057711562D-01
ALLZ(3,1,1) =1.098175104D-01
ALLC(1,1,1) =1.543289673D-01
ALLC(2,1,1) =5.353281423D-01
ALLC(3,1,1) =4.446345422D-01
!     2s
ALLZ(1,2,1) =2.581578398D00
ALLZ(2,2,1) =1.567622104D-01
ALLZ(3,2,1) =6.018332272D-02
ALLC(1,2,1) =-5.994474934D-02
ALLC(2,2,1) =5.960385398D-01
ALLC(3,2,1) =4.581786291D-01
!     2p
ALLZ(1,2,2) =9.192379002D-01
ALLZ(2,2,2) =2.359194503D-01
ALLZ(3,2,2) =8.009805746D-02
ALLC(1,2,2) =1.623948553D-01
ALLC(2,2,2) =5.661708862D-01
ALLC(3,2,2) =4.223071752D-01
!     3s
ALLZ(1,3,1) =5.641487709D-01
ALLZ(2,3,1) =6.924421391D-02
ALLZ(3,3,1) =3.269529097D-02
ALLC(1,3,1) =-1.782577972D-01
ALLC(2,3,1) =8.612761663D-01
ALLC(3,3,1) =2.261841969D-01
!     3p
ALLZ(1,3,2) =2.692880368D00
ALLZ(2,3,2) =1.489359592D-01
ALLZ(3,3,2) =5.739585040D-02
ALLC(1,3,2) =-1.061945788D-02
ALLC(2,3,2) =5.218564264D-01
ALLC(3,3,2) =5.450015143D-01
!     4s
ALLZ(1,4,1) =2.267938753D-01
ALLZ(2,4,1) =4.448178019D-02
ALLZ(3,4,1) =2.195294664D-02
ALLC(1,4,1) =-3.349048323D-01
ALLC(2,4,1) =1.056744667D00
ALLC(3,4,1) =1.256661680D-01
!     4p
ALLZ(1,4,2) =4.859692220D-01
ALLZ(2,4,2) =7.430216918D-02
ALLZ(3,4,2) =3.653340923D-02
ALLC(1,4,2) =-6.147823411D-02
ALLC(2,4,2) =6.604172234D-01
ALLC(3,4,2) =3.932639495D-01
!     5s
ALLZ(1,5,1) =1.080198458D-01
ALLZ(2,5,1) =4.408119382D-02
ALLZ(3,5,1) =2.610811810D-02
ALLC(1,5,1) =-6.617401158D-01
ALLC(2,5,1) =7.467595004D-01
ALLC(3,5,1) =7.146490945D-01
!     5p
ALLZ(1,5,2) =2.127482317D-01
ALLZ(2,5,2) =4.729648620D-02
ALLZ(3,5,2) =2.604865324D-02
ALLC(1,5,2) =-1.389529695D-01
ALLC(2,5,2) =8.076691064D-01
ALLC(3,5,2) =2.726029342D-01
!     3d
ALLZ(1,3,3) =5.229112225D-1
ALLZ(2,3,3) =1.639595876D-1
ALLZ(3,3,3) =6.386630021D-2
ALLC(1,3,3) =1.686596060D-1
ALLC(2,3,3) =5.847984817D-1
ALLC(3,3,3) =4.056779520D-1             
!     4d
ALLZ(1,4,3) =1.777717219D-1
ALLZ(2,4,3) =8.040647350D-2
ALLZ(3,4,3) =3.949855551D-2
ALLC(1,4,3) =2.308552718D-1
ALLC(2,4,3) =6.042409177D-1
ALLC(3,4,3) =2.595768926D-1
!     5d
ALLZ(1,5,3) =4.913352950D-1
ALLZ(2,5,3) =7.329090601D-2
ALLZ(3,5,3) =3.594209290D-2
ALLC(1,5,3) =-2.010175008D-2
ALLC(2,5,3) =5.899370608D-1
ALLC(3,5,3) =4.658445960D-1

iam=l-1
!
!     set parameters into functional expression
!
DO J=1,nprim
   cont(J)=ALLC(J,N,L)
   expo(J)=ALLZ(J,N,L)*zeta**2
   XNORM=(2.D0*EXPO(J)/PI)**0.75D0*(4.D0*EXPO(J))**(IAM/2.D0)/ & 
     &  SQRT(DEX(2*IAM-1))
   cont(j)=cont(j)*xnorm                
ENDDO      

end subroutine setsto3
