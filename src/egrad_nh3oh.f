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

      subroutine initialize_NH3OH()
        implicit double precision (a-h,o-z)
        DOUBLE PRECISION :: r, energy, dvdr
        DOUBLE PRECISION CART,ANUZERO
        INTEGER NULBL,NFLAG,NASURF,NDER
        INTEGER N3ATOM, NATOM, ISURF, JSURF
!
        PARAMETER (NATOM=25)
        PARAMETER (ISURF = 5)
        PARAMETER (JSURF = ISURF*(ISURF+1)/2)

        COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +                DIJCART(NATOM,3,JSURF)

        COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER


        INTEGER :: init
        data init /0/
        save init
        natoms = 6
        nder =  1
        nasurf =  0
        nasurf(1,1) = 1
        nflag(1) = 1
        nflag(2) = 1
        nflag(3) = 0
        nflag(4) = 0
        nflag(5) = 0

        init = 0
        call prepot_nh3oh
        init = 1

      end subroutine initialize_nh3oh

      subroutine egrad_nh3oh(q,Natoms,Nbeads,V,dVdq,info)
        implicit double precision (a-h,o-z)
        INTEGER :: Natoms, Nbeads
        DOUBLE PRECISION :: q(3,Natoms,Nbeads)
        DOUBLE PRECISION :: dVdq(3,Natoms,Nbeads), V(Nbeads)
        INTEGER :: k
        DOUBLE PRECISION :: r, energy, dvdr
        DOUBLE PRECISION CART,ANUZERO
        INTEGER NULBL,NFLAG,NASURF,NDER
        INTEGER N3ATOM, NATOM, ISURF, JSURF
        INTEGER :: INFO
!
        PARAMETER (NATOM=25)
        PARAMETER (ISURF = 5)
        PARAMETER (JSURF = ISURF*(ISURF+1)/2)

        COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +                DIJCART(NATOM,3,JSURF)

        COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER


        INTEGER :: init
        data init /0/
        save init
        natoms = 6
        nder =  1
        nasurf =  0
        nasurf(1,1) = 1
        info = 0
        nflag(1) = 1
        nflag(2) = 1
        nflag(3) = 0
        nflag(4) = 0
        nflag(5) = 0


        init = 1
        v = 0.d0
        dvdq(:,:,:)  = 0.d0
        dvdr = 0.d0
        cart(:,:) = 0.d0
        do k = 1, Nbeads
            call qtocart_nh3oh(q,Natoms,natom,Nbeads,k,cart)
            call pot_nh3oh

            V(k) = PENGYGS
        do i = 1,3
        do j = 1,natoms
        dvdq(i,j,k) = dgscart(j,i)
        end do
        end do
        END DO
      end subroutine egrad_nh3oh

      
      subroutine qtocart_nh3oh(q0,nat,natom,nb,inb,cart)
        implicit none
        DOUBLE PRECISION :: q0(3,nat,nb), cart(natom,3)
        INTEGER :: i,j,k,inb
        INTEGER :: nat, natom,nb
        do i= 1,3
        do j = 1,nat
        cart(j,i) = q0(i,j,inb)
        end do
        end do
    !
        return
      end subroutine qtocart_nh3oh

c**************************************************************************
C   System:           NH3OH
C   Number of electronic surfaces: 1
C   Number of derivatives: 0
C   Number of bodies: 6
C   Interface:        potlib2001
C   Common name:      NH3OH
C   Notes:
C
C   References: M. Monge-Palacios, C. Rangel and J. Espinosa-Garcia, 
C   J.Chem.Phys. 138, 084305 (2013) 


      SUBROUTINE POT_nh3oh
C
C        This potential is written such that:
C
C                       X(1)  - X(3)  : X, Y, Z for H1
C                       X(4)  - X(6)  : X, Y, Z for N
C                       X(7)  - X(9)  : X, Y, Z for H2
C                       X(10) - X(12) : X, Y, Z for H3
C                       X(13) - X(15) : X, Y, Z for O (Hb)
C                       X(16) - X(18) : X, Y, Z for H 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
c
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
      CALL CARTOU_nh3oh
      CALL CARTTOR_nh3oh
C
C Changing units from angstrom to bohr and initialize
C
      DO 10 I = 1, 18
         q(I) = R(I)*0.52918d0
         pdot(I)=0.D0
 10   CONTINUE
!      stop "Hphpi"
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0

       call coorden_nh3oh
c
c  calculate switching functions
c
       call switchf_nh3oh
c
c  calculate reference angles and their derivatives
c
       call refangles_nh3oh
       call oprefangles_nh3oh
c
c  calculate stretching potential
c
       call stretch_nh3oh(vstr)
c
c  calculate out of plane bending potential
c
c       call opbend(vop)
c
c  calculate in plane bending potential
c
       call ipbend_nh3oh(vip)
c
c  total potential energy is vstr+vop+vip
c
       en=vstr+vip
c      en=vstr+vop+vip
c
c 0.03812 conversion factor from 10(5) j/mol to hartrees
c
       en = en*0.03812D0
       ENGYGS = en
c
      CALL EUNITZERO_nh3oh
      IF(NDER.NE.0) THEN
c
c 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
c
         do i=1,18
            DEGSDR(i)=pdot(i)*0.0201723d0
         enddo
c         write(*,*) "degsdr",degsdr(1:18)
c         stop "Hphpi"
c JEG Calculo de derivadas numericas de la energia con
c     respecto a las coordenadas
      PASO= 1.0D-5
      do I=1,18
        q(I)=q(I) + PASO
        call coorden_nh3oh
        call switchf_nh3oh
        call refangles_nh3oh
        call stretch_nh3oh(vstr)
        call ipbend_nh3oh(vip)
        en=vstr+vip
        en=en*0.03812D0
        DEGSDR(I)=(en-ENGYGS)/PASO
        DEGSDR(I)=DEGSDR(I)*0.52918d0
        q(I)=q(I)-PASO
      enddo
         CALL RTOCART_nh3oh
         CALL DEDCOU_nh3oh
      ENDIF
C
       return
       end
c
c******************************************************
c
       subroutine coorden_nh3oh
c
c  calculates relative coordinates and bond lengths
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
c  calculate relative coordinates
c
       do ind=1,3
         tcb(ind)=q(nc(ind))-q(nhb(ind))
         do i=1,3
           tch(i,ind)=q(nc(ind))-q(nh(i,ind))
           tbh(i,ind)=q(nhb(ind))-q(nh(i,ind))
         enddo
       enddo
       do ind=1,3
         tno(ind)=q(no(ind))-q(nhb(ind)) 
       end do
c
c  calculate bond lengths
c
       rcb=sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
       rno=sqrt(tno(1)*tno(1)+tno(2)*tno(2)+tno(3)*tno(3))
       do i=1,3
         rch(i)=sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+
     *                tch(i,3)*tch(i,3))
         rbh(i)=sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+
     *                tbh(i,3)*tbh(i,3))
       enddo
c jcc-2010
       argmax = 19.d0
       P1 = 1.d0
       do i = 1, 3
          argp1 = (w1*(rch(i)-w2))
          if (argp1.lt.argmax) then
            t1tmp = 1 - tanh(argp1)
          else
            t1tmp = 0.d0
          endif
          P1 = P1 * t1tmp
       enddo
       r0ch = P1 * r0chr + (1-P1) * r0chp
c jcc-2010
       argmax = 19.d0
       P2 = 1.d0
          argp2 = (w3*(rno-w4))
          if (argp2.lt.argmax) then
            t2tmp = 1 - tanh(argp2)
          else
            t2tmp = 0.d0
          endif
          P2 = P2 * t2tmp
       r0hh = P2 * r0hhr + (1-P2) * r0hhp
       return
       end
c
c******************************************************
c
c
       subroutine refangles_nh3oh
c
c  subroutine calculates reference angles for the "in-plane" potential
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
C      tau=acos(-1.0d0/3.0d0)
C      pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
       ppito=(twopi-taunh2)/2.0d0
c
c  set diagonal elements to zero
c
       do i=1,3
         theta0(i,i)=0.0d0
         do k=1,4
           dtheta0(i,i,k)=0.0d0
         enddo
       enddo
c
c  calculate reference angles
c
       theta0(1,2)=tau+(tau-ppito)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-taunh2)*(stheta(3)-1.0d0)
       theta0(1,3)=tau+(tau-ppito)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-taunh2)*(stheta(2)-1.0d0)
       theta0(2,3)=tau+(tau-ppito)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-taunh2)*(stheta(1)-1.0d0)
c
c  calculate the derivatives of theta0(i,j) in terms of rch(k)
c  quantity calulated is dtheta0(i,j,k)
c
c  derivatives wrt rch(1)
c
       dtheta0(1,2,1)=(tau-ppito)*dsphi(1)*sphi(2)
       dtheta0(1,3,1)=(tau-ppito)*dsphi(1)*sphi(3)
       dtheta0(2,3,1)=(tau-taunh2)*dstheta(1)
c
c  derivatives wrt rch(2)
c
       dtheta0(1,2,2)=(tau-ppito)*sphi(1)*dsphi(2)
       dtheta0(1,3,2)=(tau-taunh2)*dstheta(2)
       dtheta0(2,3,2)=(tau-ppito)*dsphi(2)*sphi(3)
c
c  derivatives wrt rch(3)
c
       dtheta0(1,2,3)=(tau-taunh2)*dstheta(3)
       dtheta0(1,3,3)=(tau-ppito)*sphi(1)*dsphi(3)
       dtheta0(2,3,3)=(tau-ppito)*sphi(2)*dsphi(3)
c
c  fill in the other half of the matrix
c
        do i=1,2
          do j=i+1,3
            theta0(j,i)=theta0(i,j)
            do k=1,3
              dtheta0(j,i,k)=dtheta0(i,j,k)
            enddo
          enddo
        enddo
       return
       end
c
c******************************************************
c
c
       subroutine oprefangles_nh3oh
c
c  subroutine calculates reference angles for the "out-plane" potential
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
C       pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
       ppito=(twopi-taunh2)/2.0d0
c
c  set diagonal elements to zero
c
       do i=1,2
         optheta0(i,i)=0.0d0
         do k=1,3
           doptheta0(i,i,k)=0.0d0
         enddo
       enddo
c
c  calculate reference angles
c
       optheta0(1,1)=optau+(optau-halfpi)*(sphi(1)
     *          -1.0d0)+(optau-halfpi)*(stheta(3)-1.0d0)
     *          +(optau-halfpi)*(stheta(2)-1.0d0)
       optheta0(2,2)=optau+(optau-halfpi)*(sphi(2)
     *          -1.0d0)+(optau-halfpi)*(stheta(3)-1.0d0)
     *          +(optau-halfpi)*(stheta(1)-1.0d0)
       optheta0(3,3)=optau+(optau-halfpi)*(sphi(3)
     *          -1.0d0)+(optau-halfpi)*(stheta(1)-1.0d0)
     *          +(optau-halfpi)*(stheta(2)-1.0d0)
c
       optheta0(1,2)=optau+(optau-halfpi)*(sphi(1)*sphi(2)
     *          -1.0d0)+(optau-halfpi)*(stheta(3)-1.0d0)
       optheta0(1,3)=optau+(optau-halfpi)*(sphi(1)*sphi(3)
     *          -1.0d0)+(optau-halfpi)*(stheta(2)-1.0d0)
       optheta0(2,3)=optau+(optau-halfpi)*(sphi(2)*sphi(3)
     *          -1.0d0)+(optau-halfpi)*(stheta(1)-1.0d0)
c
c  calculate the derivatives of optheta0(i,j) in terms of rch(k)
c  quantity calulated is doptheta0(i,j,k)
c
c  derivatives wrt rch(1)
c
       doptheta0(1,1,1)=(optau-halfpi)*dsphi(1)
       doptheta0(2,2,1)=(optau-halfpi)*dstheta(1)
       doptheta0(3,3,1)=(optau-halfpi)*dstheta(1)
c
       doptheta0(1,2,1)=(optau-halfpi)*dsphi(1)*
     * sphi(2)
       doptheta0(1,3,1)=(optau-halfpi)*dsphi(1)*
     * sphi(3)
       doptheta0(2,3,1)=(optau-halfpi)*
     * dstheta(1)
c
c  derivatives wrt rch(2)
c
       doptheta0(1,1,2)=(optau-halfpi)*dstheta(2)
       doptheta0(2,2,2)=(optau-halfpi)*dsphi(2)
       doptheta0(3,3,2)=(optau-halfpi)*dstheta(2)
c
       doptheta0(1,2,2)=(optau-halfpi)*sphi(1)*
     * dsphi(2)
       doptheta0(1,3,2)=(optau-halfpi)*
     * dstheta(2)
       doptheta0(2,3,2)=(optau-halfpi)*dsphi(2)*
     * sphi(3)
c
c  derivatives wrt rch(3)
c
c
       doptheta0(1,1,3)=(optau-halfpi)*dstheta(3)
       doptheta0(2,2,3)=(optau-halfpi)*dstheta(3)
       doptheta0(3,3,3)=(optau-halfpi)*dsphi(3)
c
       doptheta0(1,2,3)=(optau-halfpi)*
     * dstheta(3)
       doptheta0(1,3,3)=(optau-halfpi)*sphi(1)*
     * dsphi(3)
       doptheta0(2,3,3)=(optau-halfpi)*sphi(2)*
     * dsphi(3)
c
c
c  fill in the other half of the matrix
c
        do i=1,2
          do j=i+1,3
            optheta0(j,i)=optheta0(i,j)
            do k=1,3
              doptheta0(j,i,k)=doptheta0(i,j,k)
            enddo
          enddo
        enddo
       return
       end
c
c******************************************************
c
c
       subroutine stretch_nh3oh(vstr)
c
c  subroutine to calculate leps-type stretching potential and its
c  derivatives
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension vqch(4),vjch(4),vqbh(4),vjbh(4),vq(4),vj(4),
     *           achdc(3),achdh(4,3)
c
c  calculate avergage bond length for the nh3 moiety
c
       rav=(rch(1)+rch(2)+rch(3))/3.0d0
c
c  initialise:
c
       vstr=0.0d0
c
       texp = exp(-(4.d0*(rav-rcbsp)/b3cb)**4.d0)
       d3cb = (d3cbi-a3cb) + a3cb * texp
c      dd3cb = -4.d0*a3cb*texp*(rav-rcbsp)**3.d0*(4.d0/b3cb)**4.d0
c  ach:
c
c  in double precision tanh(19.0d0)=1.0d0 and we put the if statement
c  in to avoid overflow/underflow errors
c
       arga=c1ch*(rav-r0ch)
       if(arga.lt.19.0d0)then
         ach=a1ch+b1ch*(tanh(arga)+1.0d0)*0.5d0
         dumach=b1ch*c1ch/(2.0d0*cosh(arga)**2)
       else
         ach=a1ch+b1ch
         dumach=0.0d0
       endif
c
c  calculate singlet: e1, triplet: e3 energies and vq and vj
c  terms for each bond
c
       e1=d1cb*(exp(-2.0d0*acb*(rcb-r0cb))-2.0d0*exp(-acb*(rcb-r0cb)))
       e3=d3cb*(exp(-2.0d0*acb*(rcb-r0cb))+2.0d0*exp(-acb*(rcb-r0cb)))
       vqcb=(e1+e3)*0.5d0
       vjcb=(e1-e3)*0.5d0
c O-H  new term (simple morse term)
c
       dt=(rno-r0hh)
       expterm=exp(-ahh*dt)
       vno=d1hh*(1.d0-expterm)**2.d0

       do i=1,3
         e1=d1ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              -2.0d0*exp(-ach*(rch(i)-r0ch)))
         e3=d3ch*(exp(-2.0d0*ach*(rch(i)-r0ch))
     *              +2.0d0*exp(-ach*(rch(i)-r0ch)))
         vqch(i)=(e1+e3)*0.5d0
         vjch(i)=(e1-e3)*0.5d0
         e1=d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         e3=d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         vqbh(i)=(e1+e3)*0.5d0
         vjbh(i)=(e1-e3)*0.5d0
c
c  calculate 3 body potential
c
         vq(i)=vqch(i)+vqcb+vqbh(i)
         vj(i)=-sqrt(((vjch(i)-vjcb)**2+(vjcb-vjbh(i))**2
     *                 +(vjbh(i)-vjch(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
       vstr=vstr+vno
c
c  partial derivatives
c  first we need the derivative of ach:
c
       do ind=1,3
         achdc(ind)=dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2)
     *            +tch(3,ind)/rch(3))/3.0d0
         do i=1,3
           achdh(i,ind)=-dumach*tch(i,ind)/rch(i)/3.0d0
         enddo
       enddo
       dumqcb=-acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))-
     *         (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
c
c  calculate cartesian derivatives:
c  looping over ch(i) and bh(i)
c
       do i=1,3
         dumqbh=-ahh*((d1hh+d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))-
     *           (d1hh-d3hh)*exp(-ahh*(rbh(i)-r0hh)))/rbh(i)
         factj=0.5d0/vj(i)
         dumjcb=-acb*((d1cb-d3cb)*exp(-2.0d0*acb*(rcb-r0cb))
     *            -(d1cb+d3cb)*exp(-acb*(rcb-r0cb)))*factj/rcb
         dumjbh=-ahh*((d1hh-d3hh)*exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *            -(d1hh+d3hh)*exp(-ahh*(rbh(i)-r0hh)))*factj/rbh(i)
         do ind=1,3
c
c  deriv wrt hb:
c
                  pdot(nhb(ind))=pdot(nhb(ind))
     *             -tcb(ind)*dumqcb+tbh(i,ind)*dumqbh
     *            +(vjch(i)-vjcb)*(dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*(-dumjcb*tcb(ind)-dumjbh*tbh(i,ind))
     *            +(vjbh(i)-vjch(i))*dumjbh*tbh(i,ind)
c
c  dvqch(i)/dc
c
           dumqch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
               pdot(nc(ind))=pdot(nc(ind))+dumqch+tcb(ind)*dumqcb
c
c  dvqch(i)/dh(i)
c
           dumqhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                 -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
              pdot(nh(i,ind))=pdot(nh(i,ind))+dumqhi-tbh(i,ind)*dumqbh
c
c  dvjch(i)/dc
c
           dumjch=-(ach*tch(i,ind)/rch(i)+achdc(ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnc(ind)
c
           pdot(nc(ind))=pdot(nc(ind))
     *            +(vjch(i)-vjcb)*(dumjch-dumjcb*tcb(ind))
     *            +(vjcb-vjbh(i))*dumjcb*tcb(ind)
     *            -(vjbh(i)-vjch(i))*dumjch
c
c  dvjch(i)/dh(i)
c
           dumjhi=(ach*tch(i,ind)/rch(i)-achdh(i,ind)*(rch(i)-r0ch))
     *              *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *               -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
c
c  dvj(i)/dnh(i,ind)
c
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *            +(vjch(i)-vjcb)*dumjhi
     *            +(vjcb-vjbh(i))*dumjbh*tbh(i,ind)
     *            +(vjbh(i)-vjch(i))*(-dumjbh*tbh(i,ind)-dumjhi)
c
c  dv(i)/dh(j)
c
           do k=1,2
             j=i+k
             if(j.gt.3)j=j-3
             dumqhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch+d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                    -(d1ch-d3ch)*exp(-ach*(rch(i)-r0ch)))
             dumjhj=-achdh(j,ind)*(rch(i)-r0ch)
     *                 *((d1ch-d3ch)*exp(-2.0d0*ach*(rch(i)-r0ch))
     *                  -(d1ch+d3ch)*exp(-ach*(rch(i)-r0ch)))*factj
             pdot(nh(j,ind))=pdot(nh(j,ind))+dumqhj
     *            +(vjch(i)-vjcb)*dumjhj
     *            -(vjbh(i)-vjch(i))*dumjhj
           enddo
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine opbend_nh3oh(vop)
c
c  subroutine calculates symmetrized vop potential and derivatives
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       double precision norma
       dimension sumd2(4),sumd4(4)
       dimension in(3),a(3),b(3),axb(3),c(4,3),argd(4)
c
c
       vop=0.0d0
c
c  calculate force constants and their derivatives
c
       call opforce_nh3oh
c
c  calculate out-of-plane angle and derivatives
c
       do i=1,3
         j=i+1
         if(j.gt.3)j=j-3
         k=j+1
         if(k.gt.3)k=k-3
c
c  modification to ensure that the set of methane CH bond vector
c  (rj,rk,rl) is a right-handed set
c
       in(1)=i
       in(2)=j
       in(3)=k
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(j,ind))-q(nh(i,ind))
         b(ind)=q(nh(k,ind))-q(nh(i,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
c
c  if argd < 0 we need to switch vectors k and l around
c
         if (argd(in(ii)).lt.0.d0) then
             itemp=j
             j=k
             k=itemp
         endif
       enddo
c
c  sum2 = sum delta**2
c  sum4 = sum delta**4
c
         call calcdelta_nh3oh(i,j,k,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
       enddo
c
       do i=1,3
         do j=1,3
c
c  overall derivatives of force constants i wrt the bond-length rch(j)
c
           ddr=dfdelta(i,j)*sumd2(i)+dhdelta(i,j)*sumd4(i)
c
c  calculate derivatives in terms of cartesian coordinates:
c
           do ind=1,3
             pdot(nh(j,ind))=pdot(nh(j,ind))-tch(j,ind)*ddr/rch(j)
             pdot(nc(ind))=pdot(nc(ind))+tch(j,ind)*ddr/rch(j)
           enddo
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipbend_nh3oh(vip)
c
c  subroutine calculates symmetrised in plane bend term
c  and its derivatives
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension costh(4,4),theta(4,4),dth(4,4)
c
c  initialise
c
       vip=0.0d0
c
c  calculate force constants: fk0(i,j), f1(i)
c  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
c
       call ipforce_nh3oh
c
c  calculate theta(i,j) and in plane bend potential
c
       do i=1,2
         do j=i+1,3
           costh(i,j)=tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2)
     *                       +tch(i,3)*tch(j,3)
           costh(i,j)=costh(i,j)/rch(i)/rch(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
c
c  calculate partial derivatives wrt cartesian coordinates
c
c  calculate pdots wrt theta:
c
           termth=-1.0d0/sqrt(1.0d0-costh(i,j)*costh(i,j))
           do ind=1,3
             dthi=-tch(j,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(i,ind)/rch(i)/rch(i)
             dthi=dthi*termth
             dthj=-tch(i,ind)/rch(i)/rch(j)
     *                  +costh(i,j)*tch(j,ind)/rch(j)/rch(j)
             dthj=dthj*termth
             dthc=-(dthi+dthj)
             pdot(nh(i,ind))=pdot(nh(i,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthi*dth(i,j)
             pdot(nh(j,ind))=pdot(nh(j,ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthj*dth(i,j)
             pdot(nc(ind))=pdot(nc(ind))
     *                       +fk0(i,j)*f1(i)*f1(j)*dthc*dth(i,j)
             do k=1,3
c
c  calculate pdots wrt force constants and wrt theta0
c
               dth0k=-dtheta0(i,j,k)*tch(k,ind)/rch(k)
               dth0c=-dth0k
               pdot(nh(k,ind))=pdot(nh(k,ind))
     *                -0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                -0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0k*dth(i,j)
               pdot(nc(ind))=pdot(nc(ind))
     *                +0.5d0*tch(k,ind)*dfdc(i,j,k)*dth(i,j)**2/rch(k)
     *                      -fk0(i,j)*f1(i)*f1(j)*dth0c*dth(i,j)
               pdot(nhb(ind))=pdot(nhb(ind))
     *                +0.5d0*tbh(k,ind)*dfdh(i,j,k)*dth(i,j)**2/rbh(k)
             enddo
           enddo
         enddo
       enddo
c  Now, calculate h2o bending energies.
c
c  First, calculate the angles:
c
cc       do i = 1,4
       do i = 1,3
         dot = 0.d0
         do j = 1,3
           dot = dot - tno(j)*tbh(i,j)
         enddo
         cosine = dot / (rno*rbh(i))
         cosine = min(1.d0, max(-1.d0, cosine))
         angh2o(i) = acos(cosine)
       enddo
c  Now, calculate each force constant fkh2o(i) as a function of the
c O-H(I)
c  distance, rbh(i)
c
CC
CC CRD 2009
       do i=1,3
          arga= alph2o * (rbh(i) - r0hh)
          if(arga.lt.19.0d0)then
            fkh2o(i) = fkh2oeq * (1 - tanh(arga))
          else
            fkh2o(i) = 0.0d0
          endif
       enddo
c  Now calculate the contribution to the energy for each H(I)-O-H(O)
c  harmonic bending term
c
CC
CC CRD 2009
      do i=1,3
           dang = (angh2o(i) - angh2oeq)
           vip=vip+0.5d0* fkh2o(i) * dang * dang
      enddo
       return
       end
c
c*************************************************************************
c
       subroutine calcdelta_nh3oh(i,j,k,sum2,sum4)
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension  delta(4),in(3),a(3),b(3),axb(3),c(4,3),argd(4),
     *            daxb(4,3,3),cdot(4,3,3),atemp2(3)
c
c  set i,j,k indices
c
       in(1)=i
       in(2)=j
       in(3)=k
c
c  vector a is rj-ri, vector b is rk-ri
c
       do ind=1,3
         a(ind)=q(nh(j,ind))-q(nh(i,ind))
         b(ind)=q(nh(k,ind))-q(nh(i,ind))
       enddo
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-optheta0(i,in(ii))
c
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
c  derivatives of axb wrt hi:
c
       daxb(i,1,1)=0.0d0
       daxb(i,1,2)=b(3)-a(3)
       daxb(i,1,3)=-b(2)+a(2)
       daxb(i,2,1)=-b(3)+a(3)
       daxb(i,2,2)=0.0d0
       daxb(i,2,3)=b(1)-a(1)
       daxb(i,3,1)=b(2)-a(2)
       daxb(i,3,2)=-b(1)+a(1)
       daxb(i,3,3)=0.0d0
c
c  derivatives of axb wrt hj:
c
       daxb(j,1,1)=0.0d0
       daxb(j,1,2)=-b(3)
       daxb(j,1,3)=b(2)
       daxb(j,2,1)=b(3)
       daxb(j,2,2)=0.0d0
       daxb(j,2,3)=-b(1)
       daxb(j,3,1)=-b(2)
       daxb(j,3,2)=b(1)
       daxb(j,3,3)=0.0d0
c
c  derivatives of axb wrt hl:
c
       daxb(k,1,1)=0.0d0
       daxb(k,1,2)=a(3)
       daxb(k,1,3)=-a(2)
       daxb(k,2,1)=-a(3)
       daxb(k,2,2)=0.0d0
       daxb(k,2,3)=a(1)
       daxb(k,3,1)=a(2)
       daxb(k,3,2)=-a(1)
       daxb(k,3,3)=0.0d0
c
c   loop over cdot(in(ii),ind,jind) where we consider deriv of c(in(ii))
c   wrt h(in(ii),jind) with components jind
c
       do ii=1,3
c
c  deriv of cdot(in(ii),x) wrt x, y, z
c
         cdot(in(ii),1,1)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),1)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,2)=tch(in(ii),2)*c(in(ii),1)/rch(in(ii))**2
         cdot(in(ii),1,3)=tch(in(ii),3)*c(in(ii),1)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),y) wrt x, y, z
c
         cdot(in(ii),2,1)=tch(in(ii),1)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,2)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),2)*c(in(ii),2)/rch(in(ii))**2
         cdot(in(ii),2,3)=tch(in(ii),3)*c(in(ii),2)/rch(in(ii))**2
c
c  deriv of cdot(in(ii),z) wrt x, y, z
c
         cdot(in(ii),3,1)=tch(in(ii),1)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,2)=tch(in(ii),2)*c(in(ii),3)/rch(in(ii))**2
         cdot(in(ii),3,3)=1.0d0/rch(in(ii))
     *                   +tch(in(ii),3)*c(in(ii),3)/rch(in(ii))**2
       enddo
c
       do ii=1,3
         do ind=1,3
            deldot=-doptheta0(i,ii,i)
c
c  derivative wrt h(i,ind)
c  for  rch(i) only terms are from the derivatives of theta0
c
            deldot=-deldot*tch(i,ind)/rch(i)
c
c  derivative wrt c(ind)
c
            deldot=-deldot
c
           do jj=1,3
c
c  partial derivatives wrt h(in(jj),ind), loop over delta(i,in(ii))
c
c   atemp1 is axb dot daxb wrt h(in(jj))
c
            atemp1=axb(1)*daxb(in(jj),ind,1)
     *            +axb(2)*daxb(in(jj),ind,2)
     *            +axb(3)*daxb(in(jj),ind,3)
            atemp1=atemp1/(norma**3)
c
c  atemp2 is deriv of normalised axb
c
            atemp2(1)=daxb(in(jj),ind,1)/norma-atemp1*axb(1)
            atemp2(2)=daxb(in(jj),ind,2)/norma-atemp1*axb(2)
            atemp2(3)=daxb(in(jj),ind,3)/norma-atemp1*axb(3)
c
c  atemp3 is daxb dot c(in(ii))

            atemp3=atemp2(1)*c(in(ii),1)+atemp2(2)*c(in(ii),2)
     *                             +atemp2(3)*c(in(ii),3)
c
c  atemp4 is axb dot cdot
c
            atemp4=0.0d0
            if(ii.eq.jj)then
c
c  ie deriv of c(in(ii)) wrt h(in(jj)) is non zero only for ii = jj
c
              atemp4=axb(1)*cdot(in(ii),1,ind)
     *                     +axb(2)*cdot(in(ii),2,ind)
     *                     +axb(3)*cdot(in(ii),3,ind)
              atemp4=atemp4/norma
            endif
c
c  atemp5 is deriv of theta0(i,in(ii)) wrt to nh(in(jj),ind)
c
            atemp5=-doptheta0(i,in(ii),in(jj))
c
c  deriv wrt h(in(jj)),ind):
c
            atemp5=-atemp5*tch(in(jj),ind)/rch(in(jj))
            deldot=atemp3+atemp4
            deldot=-1.0d0/sqrt(1.0d0-argd(in(ii))**2)*deldot
            deldot=deldot+atemp5
            pdot(nh(in(jj),ind))=pdot(nh(in(jj),ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  for carbon the only contributions are from axb dot cdot term and
c  from theta0 and derivative cdot wrt carbon=-cdot wrt hydrogen
c
            deldot=1.0d0/sqrt(1.0d0-argd(in(ii))**2)*atemp4
            deldot=deldot-atemp5
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
          enddo
        enddo
       enddo
       return
       end

c******************************************************
c
       subroutine opforce_nh3oh
c
c  calculates the out-of-plane bending force constants
c  and their derivatives
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension switch(4),dswitch(4,4)
c
c  calculate switching functions:
c
       switch(1)=s3(1)*s3(2)*s3(3)
       switch(2)=s3(2)*s3(3)*s3(1)
       switch(3)=s3(3)*s3(1)*s3(2)
c
c  calculate derivatives:
c  derivative of switch(1) wrt the 3 rch bond lengths
c
       dswitch(1,1)=ds3(1)*s3(2)*s3(3)
       dswitch(1,2)=s3(1)*ds3(2)*s3(3)
       dswitch(1,3)=s3(1)*s3(2)*ds3(3)
c
c  derivative of switch(2) wrt the 3 rch bond lengths
c
       dswitch(2,1)=s3(2)*s3(3)*ds3(1)
       dswitch(2,2)=ds3(2)*s3(3)*s3(1)
       dswitch(2,3)=s3(2)*ds3(3)*s3(1)
c
c  derivative of switch(3) wrt the 3 rch bond lengths
c
       dswitch(3,1)=s3(3)*ds3(1)*s3(2)
       dswitch(3,2)=s3(3)*s3(1)*ds3(2)
       dswitch(3,3)=ds3(3)*s3(1)*s3(2)
c
c  calculate the force constants and their derivatives
c
       do i=1,3
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
         do j=1,3
           dfdelta(i,j)=dswitch(i,j)*fch3
           dhdelta(i,j)=dswitch(i,j)*hch3
         enddo
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine ipforce_nh3oh
c
c  calculates the symmetrised in plane bend force constants and
c  all partial derivatives involving them
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension dfko(4),df1dc(4),df1dh(4)
c
c  set force constant at asymptotes
c
       f0=fkinf+ak
       f2=fkinf
c
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)-1.0d0)
c
       do i=1,3
c
c  calc derivatives of fko wrt each of the rch(i) bonds
c
          dfko(i)=-2.0d0*ak*bk*(rch(i)-r0ch)*exp(-argfk0)
c
c  calculate the terms f1(i)
c
         arga1=aa1*rbh(i)*rbh(i)
         arga2=aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rch(i)-r0ch)**2)
c
c  and calculate the derivatives wrt rch(i) and rbh(i)
c
         duma1=2.0d0*aa1*rbh(i)*exp(-arga1)
         duma2=-2.0d0*aa3*aa4*(rbh(i)-r0hh)*exp(-arga2)
         df1dc(i)=-2.0d0*(rch(i)-r0ch)*a1*a2*exp(-a2*(rch(i)-r0ch)**2)
         df1dh(i)=duma1*exp(-a2*(rch(i)-r0ch)**2)
     *             -duma2*(rch(i)-r0ch)**2*a1*exp(-a2*(rch(i)-r0ch)**2)
       enddo
c
c  derivative of total force constant f(i,j) wrt bond length rch(k)
c  is given by dfdc(i,j,k)
c
      dfdc(1,2,1)=dfko(1)*f1(1)*f1(2)+fko*df1dc(1)*f1(2)
      dfdc(1,2,2)=dfko(2)*f1(1)*f1(2)+fko*f1(1)*df1dc(2)
      dfdc(1,2,3)=dfko(3)*f1(1)*f1(2)
c
      dfdc(1,3,1)=dfko(1)*f1(1)*f1(3)+fko*df1dc(1)*f1(3)
      dfdc(1,3,2)=dfko(2)*f1(1)*f1(3)
      dfdc(1,3,3)=dfko(3)*f1(1)*f1(3)+fko*f1(1)*df1dc(3)
c
      dfdc(2,3,1)=dfko(1)*f1(2)*f1(3)
      dfdc(2,3,2)=dfko(2)*f1(2)*f1(3)+fko*df1dc(2)*f1(3)
      dfdc(2,3,3)=dfko(3)*f1(2)*f1(3)+fko*f1(2)*df1dc(3)
c
c  derivative of total force constant f(i,j) wrt bond length rbh(k)
c  is given by dfdh(i,j,k)
c
c  only non-zero derivatives are those from rbh(i) and rbh(j)
c
       dfdh(1,2,1)=fko*df1dh(1)*f1(2)
       dfdh(1,2,2)=fko*f1(1)*df1dh(2)
       dfdh(1,2,3)=0.0d0
c
       dfdh(1,3,1)=fko*df1dh(1)*f1(3)
       dfdh(1,3,2)=0.0d0
       dfdh(1,3,3)=fko*f1(1)*df1dh(3)
c
       dfdh(2,3,1)=0.0d0
       dfdh(2,3,2)=fko*df1dh(2)*f1(3)
       dfdh(2,3,3)=fko*f1(2)*df1dh(3)
c
       return
       end
c
c******************************************************
c
c
       subroutine switchf_nh3oh
c
c  calculates switching functions: s3,sphi,stheta
c  and their derivatives ds3,dsphi,dstheta
c
       implicit double precision (a-h,o-z)
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
      COMMON/PT1CM_nh3oh/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
      COMMON/PT4CM_nh3oh/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM_nh3oh/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM_nh3oh/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       a1s=1.5313681d-9
       b1s=-1.6696246d0
       a2s=1.0147402d-9
       b2s=-1.363798d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
       do i=1,3
         args1=a1s*(rch(i)-r0ch)*(rch(i)-b1s)**8
         if(args1.lt.argmax)then
           s1(i)=1.0d0-tanh(args1)
           ds1(i)=a1s*((rch(i)-b1s)**8
     *                 +8.0d0*(rch(i)-r0ch)*(rch(i)-b1s)**7)
           ds1(i)=-ds1(i)/cosh(args1)**2
         else
           s1(i)=0.0d0
           ds1(i)=0.0d0
         endif
c
         args2=a2s*(rch(i)-r0ch)*(rch(i)-b2s)**6
         if(args2.lt.argmax)then
           s2(i)=1.0d0-tanh(args2)
           ds2(i)=a2s*((rch(i)-b2s)**6
     *                 +6.0d0*(rch(i)-r0ch)*(rch(i)-b2s)**5)
           ds2(i)=-ds2(i)/cosh(args2)**2
         else
           s2(i)=0.0d0
           ds2(i)=0.0d0
         endif
c
c  calculate s3 and ds3
c
         args3=a3s*(rch(i)-r0ch)*(rch(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i)=1.0d0-tanh(args3)
           ds3(i)=a3s*(3.0d0*rch(i)**2-2.0d0*rch(i)*(r0ch+2.0d0*b3s)
     *          +b3s*(b3s+2.0d0*r0ch))
           ds3(i)=-ds3(i)/cosh(args3)**2
         else
           s3(i)=0.0d0
           ds3(i)=0.0d0
         endif
c
c  calculate sphi and dsphi
c
c  condition here is on the bondlength rch(i)
c  st argsphi is lt approx 19.0d0
c
         if(rch(i).lt.3.8d0)then
           argsphi=aphi*(rch(i)-r0ch)*exp(bphi*(rch(i)-cphi)**3)
           sphi(i)=1.0d0-tanh(argsphi)
           dsphi(i)=aphi*(1.0d0+3.0d0*bphi*(rch(i)-r0ch)
     *                      *(rch(i)-cphi)**2)
           dsphi(i)=dsphi(i)*exp(bphi*(rch(i)-cphi)**3)
           dsphi(i)=-dsphi(i)/cosh(argsphi)**2
         else
           sphi(i)=0.0d0
           dsphi(i)=0.0d0
         endif
c
c  calculate stheta and dstheta
c
         if(rch(i).lt.3.8d0)then
           argstheta=atheta*(rch(i)-r0ch)*exp(btheta*(rch(i)-ctheta)**3)
           stheta(i)=1.0d0-tanh(argstheta)
           dstheta(i)=atheta*(1.0d0+3.0d0*btheta*(rch(i)-r0ch)
     *           *(rch(i)-ctheta)**2)
           dstheta(i)=dstheta(i)*exp(btheta*(rch(i)-ctheta)**3)
           dstheta(i)=-dstheta(i)/cosh(argstheta)**2
         else
           stheta(i)=0.0d0
           dstheta(i)=0.0d0
         endif
       enddo
       return
       end
C
      SUBROUTINE PREPOT_nh3oh
C
C   N3TMMN = 3 * NATOMS
C   NATOMS = the number of atoms represented by this potential function
C
C   The variable N3TMMN is the minimum value of N3TM allowed to be
C   passed by the calling routine for the number of cartesian
C   coordinates needed to represent the full system represented by this
C   potential energy surface routine.
C   N3TM must be greater than or equal to N3TMMN.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*75 REF(5)
C
      PARAMETER(N3ATOM = 75)
      PARAMETER (ISURF = 5)
      PARAMETER (JSURF = ISURF*(ISURF+1)/2)
C
      PARAMETER (PI = 3.141592653589793D0)
      PARAMETER (NATOM = 25)
      PARAMETER (N3TMMN = 18)
C
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
C
C
C  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
C
      IF(NATOMS.GT.25) THEN
         WRITE(NFLAG(18),1111)
 1111    FORMAT(2X,'STOP. NUMBER OF ATOMS EXCEEDS ARRAY DIMENSIONS')
         STOP
      END IF
C
       DO I=1,5
          REF(I) = ' '
       END DO
C
       REF(1)='M.Monge-Palacios and J. Espinosa-Garcia'
       REF(2)='unpublished results, 2012'
C
      INDEXES(1) = 1
      INDEXES(2) = 7
      INDEXES(3) = 1
      INDEXES(4) = 1
      INDEXES(5) = 8
      INDEXES(6) = 1 
C
C
C
      IRCTNT=5
C
      CALL POTINFO_nh3oh
C
      CALL ANCVRT_nh3oh
c
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc(ind)=3*nnc+icount
         nhb(ind)=3*nnb+icount
         no(ind)=3*nno+icount
         do i=1,3
           nh(i,ind)=3*nnh(i)+icount
         enddo
       enddo
c
c  convert to appropriate units:
c
c  energy   in 1.0d+05 j/mol
c  time     in 1.0d-14 s
c  distance in 1.0d-10 m
c  angles   in radians
c  mass     in amu
c
C
       fact1=0.041840d0
       fact2=6.022045d0
       fact3=2.d0*3.1415926d0/360.d0
c
       d1ch=d1ch*fact1
       d3ch=d3ch*fact1
       d1cb=d1cb*fact1
       d3cbi=d3cbi*fact1
       a3cb=a3cb*fact1
       d1hh=d1hh*fact1
       d3hh=d3hh*fact1
       fch3=fch3*fact2
       hch3=hch3*fact2
       fkinf=fkinf*fact2
       ak=ak*fact2
       fkh2oeq=fkh2oeq*fact2
       angh2oeq=angh2oeq*fact3
c
      RETURN
      END
C
      BLOCK DATA PTPACM_nh3oh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      COMMON/PT3CM_nh3oh/ EZERO(ISURF+1)
C
      COMMON/INFOCM_nh3oh/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM_nh3oh/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM_nh3oh/ nnc,nnb,nnh(4),nno,
     +               r0chr,r0chp,w1,w2,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               d1hh,d3hh,ahh,
     +               r0hhr,r0hhp,w3,w4,
     +               r0cb,d1cb,d3cbi,acb,
     +               a3s,b3s,aphi,bphi,cphi,rphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5,
     +               tau,taunh2,optau,
     +               fkh2oeq,alph2o,angh2oeq,
     +               a3cb,b3cb,rcbsp
C
C      DIMENSION COORD(N3TM),DX(N3TM)
C
       common /angles/  theta0(4,4),dtheta0(4,4,4),angh2o(4)
       common /bonds/   rcb,rch(4),rbh(4),rno,r0ch,r0hh
       common /coords/  tcb(3),tno(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4),fkh2o(3)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3),no(3)
       common /opb1/    optheta0(4,4), doptheta0(4,4,4),ppito
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
      DATA NASURF /1,35*0/
      DATA NDER /1/
       DATA NFLAG /1,1,15*0,6,0,0/
C
      DATA ANUZERO /0.0D0/
      DATA ICARTR,MSURF,MDER/1,0,1/
      DATA NULBL /25*0/
      DATA NATOMS /6/
C
       DATA nnc     /2/
       DATA nnb     /5/
       DATA nnh     /1,3,4,0/
       DATA nno     /6/
       DATA r0chr   /   1.01417d0/
       DATA r0chp   /   1.02777d0/
       DATA w1      /   3.00000d0/
       DATA w2      /   1.01417d0/
       DATA d1ch    / 125.250d0/
       DATA d3ch    /  24.300d0/
       DATA a1ch    /   2.050000d0/
       DATA b1ch    /  -0.200000d0/
       DATA c1ch    /   200.4000d0/
       DATA d1hh    / 135.250d0/
       DATA d3hh    /  32.800d0/
       DATA ahh     /   2.0500d0/
       DATA r0hhr   /   0.9710d0/
       DATA r0hhp   /   0.9595d0/
       DATA w3      /   1.00d0/
       DATA w4      /   0.973d0/
       DATA r0cb    /   1.83800d0/
       DATA d1cb    /  80.900d0/
       DATA d3cbi   /  26.700d0/
       DATA acb     /   1.4800000d0/
       DATA a3s     /   0.2419100d0/
       DATA b3s     /  -0.4068400d0/
       DATA aphi    /   2.2287900d0/
       DATA bphi    /   0.0206600d0/
       DATA cphi    /   1.5209900d0/
       DATA rphi    /   1.01397d0/
       DATA atheta  /   1.1578700d0/
       DATA btheta  /   0.0358900d0/
       DATA ctheta  /   0.7115500d0/
       DATA fch3    /   0.0840000d0/
       DATA hch3    /   0.1915000d0/
       DATA fkinf   /   0.7100000d0/
       DATA ak      /  -0.0900000d0/
       DATA bk      /   2.7132000d0/
       DATA aa1     /   0.800000d0/
       DATA aa2     /   2.509960d0/
       DATA aa3     /   3.506600d0/
       DATA aa4     /   1.500000d0/
       DATA aa5     /   0.15263d0/
       DATA tau     /   1.9022600d0/
       DATA taunh2  /   1.8046700d0/
       DATA optau   /   1.1847700d0/
       DATA fkh2oeq /   0.7100000d0/
       DATA alph2o  /   0.7250d0/
       DATA angh2oeq /  103.5968d0/
       DATA a3cb     /  1.60d0/
       DATA b3cb     /  0.011d0/
       DATA rcbsp    /  1.63349d0/
C
       END


