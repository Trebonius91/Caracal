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
!   Potential energy routine taken from the POTLIB library:
!   R. J. Duchovic, Y. L. Volobuev, G. C. Lynch, A. W. Jasper, D. G. Truhlar, T. C. Allison, 
!   A. F. Wagner, B. C. Garrett, J. Espinosa-García, and J. C. Corchado, POTLIB, 
!   http://comp.chem.umn.edu/potlib. 

      subroutine egrad_oh3(q,Natoms,Nbeads,V,dVdq,info)
      integer, intent(in) :: Natoms
      integer, intent(in) :: Nbeads
      double precision, intent(in) :: q(3,Natoms,Nbeads)
      double precision, intent(out) :: V(Nbeads)
      double precision, intent(out) :: dVdq(3,Natoms,Nbeads)

      double precision :: dVdr_r(6),dVdr_l(6),h_test(6,6)
      double precision :: R(6), dVdr(6)
      double precision :: xOH1, yOH1, zOH1, rOH1
      double precision :: xOH2, yOH2, zOH2, rOH2
      double precision :: xOH3, yOH3, zOH3, rOH3
      double precision :: xH1H2, yH1H2, zH1H2, rH1H2
      double precision :: xH1H3, yH1H3, zH1H3, rH1H3
      double precision :: xH2H3, yH2H3, zH2H3, rH2H3
      integer k, info
        info = 0
      do k = 1, Nbeads
        
        xOH1 = q(1,2,k) - q(1,1,k)
        yOH1 = q(2,2,k) - q(2,1,k)
        zOH1 = q(3,2,k) - q(3,1,k)
        rOH1 = sqrt(xOH1 * xOH1 + yOH1 * yOH1 + zOH1 * zOH1)
        R(1) = rOH1

        xOH2 = q(1,3,k) - q(1,1,k)
        yOH2 = q(2,3,k) - q(2,1,k)
        zOH2 = q(3,3,k) - q(3,1,k)
        rOH2 = sqrt(xOH2 * xOH2 + yOH2 * yOH2 + zOH2 * zOH2)
        R(2) = rOH2

        xOH3 = q(1,4,k) - q(1,1,k)
        yOH3 = q(2,4,k) - q(2,1,k)
        zOH3 = q(3,4,k) - q(3,1,k)
        rOH3 = sqrt(xOH3 * xOH3 + yOH3 * yOH3 + zOH3 * zOH3)
        R(3) = rOH3

        xH1H2 = q(1,3,k) - q(1,2,k)
        yH1H2 = q(2,3,k) - q(2,2,k)
        zH1H2 = q(3,3,k) - q(3,2,k)
        rH1H2 = sqrt(xH1H2 * xH1H2 + yH1H2 * yH1H2 + zH1H2 * zH1H2)
        R(4) = rH1H2

        xH1H3 = q(1,4,k) - q(1,2,k)
        yH1H3 = q(2,4,k) - q(2,2,k)
        zH1H3 = q(3,4,k) - q(3,2,k)
        rH1H3 = sqrt(xH1H3 * xH1H3 + yH1H3 * yH1H3 + zH1H3 * zH1H3)
        R(5) = rH1H3

        xH2H3 = q(1,4,k) - q(1,3,k)
        yH2H3 = q(2,4,k) - q(2,3,k)
        zH2H3 = q(3,4,k) - q(3,3,k)
        rH2H3 = sqrt(xH2H3 * xH2H3 + yH2H3 * yH2H3 + zH2H3 * zH2H3)
        R(6) = rH2H3
c
c       TEST: calculate internal hessian
c        R=R*0.52917721092
c        do i=1,6
c           write(*,*) i,R(i)
c        end do
c        h_test=0.d0
c        step=0.00001d0
c        do i=1,6
c           R(i)=R(i)+step
c
c           call pot_oh3(R, V(k), dVdr_r)
c
c           R(i)=R(i)-2.*step
c
c           call pot_oh3(R, V(k), dVdr_l)
c
c           do j=1,6
c              h_test(i,j)=(dVdr_r(j)-dVdr_l(j))/(2.*step) 
c           end do
c
c
c         end do 
c         write(*,*) "TEST: internal hessian:"
c         do i=1,6
c            write(*,'(6e16.7)') h_test(i,:)
c         end do
c         stop "Jjfggdd"
c       END TEST
        call pot_oh3(R, V(k), dVdr)
c        write(*,*) "Internal gradeint from OH3"
c        write(*,*) V
c        write(*,*) dVdR
c        do i=1,6
c           write(*,*) i,dVdR(i)
c        end do


        dVdq=0.d0
c        stop "Jjhgg"
c
c       The first atom (O)
c
        dVdq(1,1,k) = dVdq(1,1,k) - dVdr(1) * xOH1 / rOH1
        dVdq(2,1,k) = dVdq(2,1,k) - dVdr(1) * yOH1 / rOH1
        dVdq(3,1,k) = dVdq(3,1,k) - dVdr(1) * zOH1 / rOH1

        dVdq(1,1,k) = dVdq(1,1,k) - dVdr(2) * xOH2 / rOH2
        dVdq(2,1,k) = dVdq(2,1,k) - dVdr(2) * yOH2 / rOH2
        dVdq(3,1,k) = dVdq(3,1,k) - dVdr(2) * zOH2 / rOH2

        dVdq(1,1,k) = dVdq(1,1,k) - dVdr(3) * xOH3 / rOH3
        dVdq(2,1,k) = dVdq(2,1,k) - dVdr(3) * yOH3 / rOH3
        dVdq(3,1,k) = dVdq(3,1,k) - dVdr(3) * zOH3 / rOH3
c
c       The second atom (H1)
c

        dVdq(1,2,k) = dVdq(1,2,k) + dVdr(1) * xOH1 / rOH1
        dVdq(2,2,k) = dVdq(2,2,k) + dVdr(1) * yOH1 / rOH1
        dVdq(3,2,k) = dVdq(3,2,k) + dVdr(1) * zOH1 / rOH1

        dVdq(1,2,k) = dVdq(1,2,k) - dVdr(4) * xH1H2 / rH1H2
        dVdq(2,2,k) = dVdq(2,2,k) - dVdr(4) * yH1H2 / rH1H2
        dVdq(3,2,k) = dVdq(3,2,k) - dVdr(4) * zH1H2 / rH1H2

        dVdq(1,2,k) = dVdq(1,2,k) - dVdr(5) * xH1H3 / rH1H3
        dVdq(2,2,k) = dVdq(2,2,k) - dVdr(5) * yH1H3 / rH1H3
        dVdq(3,2,k) = dVdq(3,2,k) - dVdr(5) * zH1H3 / rH1H3
c
c       The third atom (H2)
c

        dVdq(1,3,k) = dVdq(1,3,k) + dVdr(2) * xOH2 / rOH2
        dVdq(2,3,k) = dVdq(2,3,k) + dVdr(2) * yOH2 / rOH2
        dVdq(3,3,k) = dVdq(3,3,k) + dVdr(2) * zOH2 / rOH2

        dVdq(1,3,k) = dVdq(1,3,k) + dVdr(4) * xH1H2 / rH1H2
        dVdq(2,3,k) = dVdq(2,3,k) + dVdr(4) * yH1H2 / rH1H2
        dVdq(3,3,k) = dVdq(3,3,k) + dVdr(4) * zH1H2 / rH1H2

        dVdq(1,3,k) = dVdq(1,3,k) - dVdr(6) * xH2H3 / rH2H3
        dVdq(2,3,k) = dVdq(2,3,k) - dVdr(6) * yH2H3 / rH2H3
        dVdq(3,3,k) = dVdq(3,3,k) - dVdr(6) * zH2H3 / rH2H3
c
c       The fourth atom (H3)
c

        dVdq(1,4,k) = dVdq(1,4,k) + dVdr(3) * xOH3 / rOH3
        dVdq(2,4,k) = dVdq(2,4,k) + dVdr(3) * yOH3 / rOH3
        dVdq(3,4,k) = dVdq(3,4,k) + dVdr(3) * zOH3 / rOH3

        dVdq(1,4,k) = dVdq(1,4,k) + dVdr(5) * xH1H3 / rH1H3
        dVdq(2,4,k) = dVdq(2,4,k) + dVdr(5) * yH1H3 / rH1H3
        dVdq(3,4,k) = dVdq(3,4,k) + dVdr(5) * zH1H3 / rH1H3

        dVdq(1,4,k) = dVdq(1,4,k) + dVdr(6) * xH2H3 / rH2H3
        dVdq(2,4,k) = dVdq(2,4,k) + dVdr(6) * yH2H3 / rH2H3
        dVdq(3,4,k) = dVdq(3,4,k) + dVdr(6) * zH2H3 / rH2H3

      end do
c      stop "in pooott"
      end subroutine egrad_oh3



      SUBROUTINE initialize_oh3
C
C   System:          OH3
C   Functional form: 
C   Common name:     
C   Number of derivatives: 1
C   Number of bodies: 4
C   Number of electronic surfaces: 1
C   Interface: 4-XS
C   Reference: G. C. Schatz and H. Elgersma
C              Chem. Phys. Lett. 73, 21 (1980).
C
      IMPLICIT REAL*8 (A-H,O-Z)
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /POT2CM_oh3/ R(4), ENERGY, DEDR(4)
      COMMON /CONCOM_oh3/ DE(3),BETA(3),RE(3),SATO,GAM(3),REOH,REHH,
     *   CON(7),ALP(4),CLAM(4),ACON(2)
C
      DO I=1,5                                                                 
          REF(I) = ' '                                                          
      END DO                                                                   
C                                                                               
       REF(1)='OH + H2 potential energy function'
       REF(2)='G. C. Schatz and H. Elgersma'
       REF(3)='Chem. Phys. Lett. 73, 21 (1980)'
C
       IRCTNT = 3
       INDEXES(1) = 8
       INDEXES(2) = 1
       INDEXES(3) = 1
       INDEXES(4) = 1
C
      CALL POTINFO_oh3
      CALL ANCVRT_oh3
C
      RETURN
      END
C
C*****
C
      SUBROUTINE pot_oh3(R,VTOT,DVDR)
C
C   This subprogram evaluates the OH + H2 potential energy and derivatives
C   of the potential with respect to the internal coordinates.
C   The subprogram PREPEF must be called once before any calls to this 
C   subprogram. 
C   All calculations in this subprogram are in hartree atomic units.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON/USROCM_oh3/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
      COMMON /PT1CM_oh3/  REX(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON /PT4CM_oh3/  ENGYES(ISURF),DEESDR(N3ATOM,ISURF)
      COMMON /PT5CM_oh3/  ENGYIJ(JSURF),DEIJDR(N3ATOM,JSURF)
      real(kind=8)::R(6),VTOT,DVDR(6)
C
c      COMMON /POTCM_oh3/ R(6), VTOT, DVDR(6)
      COMMON /POT2CM_oh3/ RSEND(4), ENERGY, DEDR(4)
      COMMON /CONCOM_oh3/ DE(3),BETA(3),RE(3),SATO,GAM(3),REOH,REHH,
     *   CON(7),ALP(4),CLAM(4),ACON(2)
C
C   Define the statement functions VMOR and DVMOR which evaluate the Morse
C   diatomic potential energy  and the derivatives of the Morse potential.
C
      VMOR(D,B,T,RR) = D*(1.0D0-EXP(-B*(RR-T)))**2
      DVMOR(D,B,T,RR) = 2.0D0*B*D*(1.0D0-EXP(-B*(RR-T)))*EXP(-B*(RR-T))
C
      CALL CARTOU_oh3
      CALL CARTTOR_oh3
C
c      DO I = 1,6
c         R(I) = REX(I)
c      ENDDO
C
C   Zero the array which will contain the derivatives of the energy 
C   with respect to the internal coordinates.
C
         IF (NDER .EQ. 1) THEN
             DO 10 I = 1, 6
                   DVDR(I) = 0.0D0
10           CONTINUE
         ENDIF
C
C   Calculate the Morse part of the potential energy.
      VTOT = VMOR(DE(1),BETA(1),RE(1),R(1))+VMOR(DE(2),BETA(2),RE(2),R(4
     *   ))+VMOR(DE(2),BETA(2),RE(2),R(5))
C   Calculate the derivatives of the Morse part of the potential energy.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DVMOR(DE(1),BETA(1),RE(1),R(1))
             DVDR(4) = DVDR(4)+DVMOR(DE(2),BETA(2),RE(2),R(4))
             DVDR(5) = DVDR(5)+DVMOR(DE(2),BETA(2),RE(2),R(5))
         ENDIF
C   Initialize the coordinates for the three-body LEPS part of the potential.
      RSEND(1) = R(2)
      RSEND(2) = R(3)
      RSEND(3) = R(6)
C
C   Calculate the three-body LEPS portion of the potential and update the 
C   energy term.
      CALL V3POT_oh3
      VTOT = VTOT+ENERGY
C   Update the array containing the derivatives with the LEPS derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(2) = DVDR(2)+DEDR(1)
             DVDR(3) = DVDR(3)+DEDR(2)
             DVDR(6) = DVDR(6)+DEDR(3)
         ENDIF
C   Initialize the coordinates for the H2O part of the potential for H1 and H2.
      RSEND(1) = R(1)
      RSEND(2) = R(2)
      RSEND(3) = R(4)
C   Calculate the H2O part of the potential and update the energy term.
      CALL VH2O_oh3
      VTOT = VTOT+ENERGY
C   Update the array containing the derivatives with the H2O derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DEDR(1)
             DVDR(2) = DVDR(2)+DEDR(2)
             DVDR(4) = DVDR(4)+DEDR(3)
         ENDIF
C   Initialize the coordinates for the H2O part of the potential for H1 and H3.
      RSEND(1) = R(1)
      RSEND(2) = R(3)
      RSEND(3) = R(5)
C   Calculate the H2O part of the potential and update the energy term.
      CALL VH2O_oh3
      VTOT = VTOT+ENERGY
C   Update the array containing the derivatives with the H2O derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(1) = DVDR(1)+DEDR(1)
             DVDR(3) = DVDR(3)+DEDR(2)
             DVDR(5) = DVDR(5)+DEDR(3)
         ENDIF
C   Initialize the coordinates for the four-body part of the potential.
      RSEND(1) = R(2)
      RSEND(2) = R(3)
      RSEND(3) = R(4)
      RSEND(4) = R(5)
C   Calculate the four-body part of the potential and update the energy term.
      CALL V4POT_oh3
      VTOT = VTOT+ENERGY
C   Update the array containing the derivatives.
         IF (NDER .EQ. 1) THEN
             DVDR(2) = DVDR(2)+DEDR(1)
             DVDR(3) = DVDR(3)+DEDR(2)
             DVDR(4) = DVDR(4)+DEDR(3)
             DVDR(5) = DVDR(5)+DEDR(4)
         ENDIF
C
C   Adjust the potential to the correct zero, which corresponds to
C   OH(R=RE) and H2(R=RE) at infinity.
C
      VTOT = VTOT-2.0D0*DE(2)+DE(3)
C
      ENGYGS = VTOT
      DO I = 1,6
         DEGSDR(I) = DVDR(I)
      ENDDO
C
      CALL EUNITZERO_oh3
      IF (NDER.NE.0) THEN
         CALL RTOCART_oh3
         IF (NFLAG(1)+NFLAG(2).NE.0) CALL DEDCOU_oh3
      ENDIF
C
      RETURN
      END
C
C*****
C
      SUBROUTINE V3POT_oh3
C
      IMPLICIT REAL*8 (A-H,O-Z)
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C
      COMMON /POT2CM_oh3/ R(4), ENERGY, DEDR(4)
      COMMON /CONCOM_oh3/ XDE(3), XBETA(3), XRE(3), SATO, GAM(3), REOH, 
     *                REHH, CON(7), ALP(4), CLAM(4), ACON(2)
      DIMENSION DE(3), BETA(3), RE(3), Z(3), ZPO(3), OP3Z(3), ZP3(3), 
     *          TZP3(3), TOP3Z(3), DO4Z(3), B(3), X(3), COUL(3), EXCH(3) 
C
C   Set up the values of DE, BETA, and RE for the three-body LEPS potential.
C
         DE(1)   = XDE(1)
         BETA(1) = XBETA(1)
         RE(1)   = XRE(1)
         DE(2)   = XDE(1)
         BETA(2) = XBETA(1)
         RE(2)   = XRE(1)
         DE(3)   = XDE(3)
         BETA(3) = XBETA(3)
         RE(3)   = XRE(3)
C
C   Compute useful constants.
C
      DO 10 I = 1, 3
         Z(I) = SATO
         ZPO(I)  = 1.0D0+Z(I)
         OP3Z(I)  = 1.0D0+3.0D0*Z(I)
         TOP3Z(I) = 2.0D0*OP3Z(I)
         ZP3(I)   = Z(I)+3.0D0
         TZP3(I)  = 2.0D0*ZP3(I)
         DO4Z(I)  = DE(I)/4.0D0/ZPO(I)
         B(I)     = BETA(I)*DO4Z(I)*2.0D0
   10 CONTINUE
C
C   Initialize the variable used for storing the energy.
C
      ENERGY = 0.D0
C
      DO 20 I = 1, 3
         X(I)    = EXP(-BETA(I)*(R(I)-RE(I)))
         COUL(I) = DO4Z(I)*(ZP3(I)*X(I)-TOP3Z(I))*X(I)
         EXCH(I) = DO4Z(I)*(OP3Z(I)*X(I)-TZP3(I))*X(I)
         ENERGY  = ENERGY + COUL(I)
   20 CONTINUE
      RAD = SQRT((EXCH(1)-EXCH(2))**2+(EXCH(2)-EXCH(3))**2+(EXCH(3)-EXCH
     *   (1))**2)
      ENERGY = ENERGY - RAD/SQRT(2.d0) 
C
C   Compute the derivatives of the energy with respect to the internal 
C   coordinates.
         IF (NDER .EQ. 1) THEN
             S = EXCH(1) + EXCH(2) + EXCH(3)
             DO 30 I = 1, 3
                   DEDR(I) = B(I)*X(I)*((3.0D0*EXCH(I)-S)/SQRT(2.d0)*
     *                       (OP3Z(I)*X(I)-ZP3(I))/RAD-
     *                       ZP3(I)*X(I)+OP3Z(I))
30           CONTINUE
         ENDIF
C
      RETURN
      END
C
C*****
C
      SUBROUTINE V4POT_oh3
      IMPLICIT REAL*8 (A-H,O-Z)
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /CONCOM_oh3/ DUM(22),A(4),C(4),COF(2)
      COMMON /POT2CM_oh3/ R(4),E,DEDR(4)
C
      T1 = EXP(-C(1)*(R(1)-A(1))**2-C(1)*(R(2)-A(1))**2-C(3)*(R(3)-A(3))
     *   **2-C(3)*(R(4)-A(3))**2)*COF(1)
      T2 = EXP(-C(2)*(R(1)-A(2))**2-C(2)*(R(2)-A(2))**2-C(4)*(R(3)-A(4))
     *   **2-C(4)*(R(4)-A(4))**2)*COF(2)
      E = T1+T2
         IF (NDER .EQ. 1) THEN
             DEDR(1) = -2.0D0*(T1*C(1)*(R(1)-A(1))+T2*C(2)*(R(1)-A(2)))
             DEDR(2) = -2.0D0*(T1*C(1)*(R(2)-A(1))+T2*C(2)*(R(2)-A(2)))
             DEDR(3) = -2.0D0*(T1*C(3)*(R(3)-A(3))+T2*C(4)*(R(3)-A(4)))
             DEDR(4) = -2.0D0*(T1*C(3)*(R(4)-A(3))+T2*C(4)*(R(4)-A(4)))
         ENDIF
C
      RETURN
      END
C*****
C
      SUBROUTINE VH2O_oh3
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION S(3),Q(3),DQ(3),X(3),DP(3)
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /CONCOM_oh3/ DUM1(10),GAM(3),REOH,REHH,C(7),DUM2(10)
      COMMON /POT2CM_oh3/ R(4),E,DEDR(4)
C
C     XMAX1 FOR WHEN TANH SET=1.0 ON VAX
C     XMAX2 FOR PREVENTING OVERFLOWS ON VAX
C
      DATA XMAX1,XMAX2 / 15.0D0,43.0D0 /
C
C     STATEMENT FUNCTION
C
      TANLG(XX) = 2.0D0/(1.0D0+EXP(2.0D0*XX))
      S(1) = R(1)-REOH
      S(3) = R(2)-REOH
      S(2) = R(3)-REHH
      DO 30 I = 1, 3
         X(I) = 0.5D0*GAM(I)*S(I)
         Q(I) = 1.0D0-TANH(X(I))
         IF (X(I).LT.XMAX1) GO TO 20
         IF (X(I).LT.XMAX2) GO TO 10
         Q(I) = 0.0D0
         DQ(I) = 0.0D0
         GO TO 30
   10    Q(I) = TANLG(X(I))
   20    DQ(I) = -0.5D0*GAM(I)/COSH(X(I))**2
   30 CONTINUE
      P = C(1)+C(2)*(S(1)+S(3))+C(3)*S(2)+0.5D0*C(4)*(S(1)*S(1)+S(3)*S(
     *    3))+0.5D0*C(5)*S(2)*S(2)+C(6)*S(2)*(S(1)+S(3))+C(7)*S(1)*S(3)
      E = Q(1)*Q(2)*Q(3)*P
      IF (NDER .EQ. 1) THEN 
          DP(1) = C(2)+C(4)*S(1)+C(6)*S(2)+C(7)*S(3)
          DP(2) = C(3)+C(5)*S(2)+C(6)*(S(1)+S(3))
          DP(3) = C(2)+C(4)*S(3)+C(6)*S(2)+C(7)*S(1)
          DO 40 I = 1, 3
                TRM1 = 0.0D0
                IF (Q(I).EQ.0.0D0) GO TO 40
                TRM1 = DQ(I)/Q(I)
                DEDR(I) = E*(TRM1+(DP(I)/P))
   40           CONTINUE
          TEMP = DEDR(2)
          DEDR(2) = DEDR(3)
          DEDR(3) = TEMP
      ENDIF
C
      RETURN
      END
C
C*****
C
      BLOCK DATA PTPACM_oh3
C                                                                               
      IMPLICIT REAL*8 (A-H,O-Z)                                       
C                                                                               
      CHARACTER*75 REF(5)                                                       
C                                                                               
      PARAMETER(N3ATOM = 75)                                                    
      PARAMETER (ISURF = 5)                                                     
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)                             
C                                                                               
      PARAMETER (PI = 3.141592653589793D0)                                      
      PARAMETER (NATOM = 25)                                                    
C                                                                               
      COMMON/PT3CM_oh3/ EZERO(ISURF+1)                                              
C                                                                               
      COMMON/INFOCM_oh3/ CARTNU(NATOM,3),INDEXES(NATOM),                            
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF                        
C                                                                               
      COMMON/USRICM_oh3/ CART(NATOM,3),ANUZERO,                                     
     +               NULBL(NATOM),NFLAG(20),                                    
     +               NASURF(ISURF+1,ISURF+1),NDER                               
C                                                                               
      COMMON /CONCOM_oh3/ DE(3),BETA(3),RE(3),SATO,GAM(3),REOH,REHH,
     *    CON(7), ALP(4),CLAM(4),ACON(2)
C
      DATA NASURF /1,35*0/                                                      
      DATA NDER /1/                                                             
      DATA NFLAG /1,1,15*0,6,0,0/                                              
C                                                                               
      DATA ANUZERO /0.0D0/                                                      
      DATA ICARTR,MSURF,MDER/2,0,1/                                             
      DATA NULBL /25*0/                                                         
      DATA NATOMS /4/                                                           
C
      DATA DE / 0.148201D0, 0.0275690D0, 0.151548D0 /  
      DATA BETA / 1.260580D0, 0.924180D0, 1.068620D0 /  
      DATA RE / 1.863300D0, 2.907700D0, 1.428600D0 /
      DATA SATO / 0.10D0 /
      DATA GAM / 2.399700D0, 1.058350D0, 2.399700D0 /
      DATA REOH / 1.808090D0 /
      DATA REHH / 2.861590D0 /
      DATA CON / -.0015920D0, 0.026963D0, 0.0014689D0, 0.080011D0,
     *           0.085816D0, -0.063179D0, 0.101380D0 /
      DATA ALP / 4.773D0, 7.14D0, 2.938D0, 5.28D0 /
      DATA CLAM / 0.10D0, 0.10D0, 0.20D0, 0.03D0 /
      DATA ACON / 0.10D0, 0.009D0 /
C
      END
C
C*****



