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

       subroutine egrad_c2h7(q,Natoms,Nbeads,V,dVdq,info) 
       implicit none
        INTEGER :: Natoms, Nbeads, info
        DOUBLE PRECISION :: q(3,Natoms,Nbeads)
        DOUBLE PRECISION :: dVdq(3,Natoms,Nbeads), V(Nbeads)
        double precision :: X(Natoms,3), energy
	double precision :: shift, dvdqr(3,Natoms), dvdql(3,Natoms)
        integer :: i, j, k, i1, j1, k1
	DOUBLE PRECISION :: rH5H6,rH5H7,rH5H8, rmin1
	shift = 1.d-4
	info = 0 
	do k = 1, Nbeads
	X(:,:) = 0.d0
	do i=1,3
	do j=1,Natoms
	X(j,i)=q(i,j,k)
	end do
	end do 
        rH5H6 = 0.d0
        rH5H7 = 0.d0 
        rH5H8 = 0.d0 
	do i1 = 1,3
	rH5H6 = rH5H6 + (X(5,i1)-X(6,i1))**2.d0
	end do
	do i1 = 1,3
	rH5H7 = rH5H7 + (X(5,i1)-X(7,i1))**2.d0
	end do
	do i1 = 1,3
	rH5H8 = rH5H8 + (X(5,i1)-X(8,i1))**2.d0
	end do
	rmin1 = max(rH5H6,rH5H7,rH5H8)
	rmin1 = dsqrt(rmin1)
!	if (rmin1.gt.200.d0) then 
!	info = 0 !DEBUG 
!	v = 0.d0
!	dvdq = 0.d0
!	return
!	end if 

	call c2h7pot_au(X,energy)
	V(k) = energy 
	do i1 = 1,3
	do j1 = 1,Natoms
	X(j1,i1) = X(j1,i1) - shift		 
	call c2h7pot_au(X,energy)
	dvdql(i1,j1) = energy 
	x(j1,i1) = x(j1,i1) + 2.d0*shift
	call c2h7pot_au(X,energy)
	dvdqr(i1,j1) = energy 
	dvdq(i1,j1,k) = (dvdqr(i1,j1) - dvdql(i1,j1))/(2.d0*shift)
	x(j1,i1) = x(j1,i1) - shift		 
	end do	
	end do
	END DO
      return
      end 
C

C**************************************************************
      SUBROUTINE c2h7pot_au(Xau,eng_au)
C**************************************************************
       implicit none 
        double precision :: X(9,3), Xau(9,3)
      integer :: i,j
      double precision :: eng_kcal, eng_au
C
	X(:,:) = 0.d0
	eng_kcal = 0.d0
	eng_au = 0.d0 	
C converting to Angstrom
	do i = 1,3
	do j = 1,9
	X(j,i) = Xau(j,i)*0.5291772083d0
	end do 
	end do
!	end do 
        call c2h7pot(X,eng_kcal)
C converting energy to a.u.
	eng_au = eng_kcal/627.5095d0
C
      return
      end 
C
C   System:          C2H7
C   Functional form: 
C   Common name:     
C   Number of derivatives: 0
C   Number of bodies: 9
C   Number of electronic surfaces: 1
C   Interface: Section-2
C   Data file:
C
C   References:: Arindam Chakraborty, Yan Zhao, Hai Lin, and Donald G. Truhlar,
C      J. Chem. Phys., 124, 044315 (2006)
C
C   Notes:    PES for C2H6 + H --> C2H5 + H2
C
C
C     H2 H6
C     |  |
C  H4-C1-C5-H8.....H9
C     |  |
C     H3 H7
C
C input coordinate matrix: X(9,3)
C Coordinates in Ang
C Energy in kcal
C**************************************************************
      SUBROUTINE c2h7pot(X,eng_kcal)
C**************************************************************
C
      implicit none
      integer npmm
      integer npjg
      parameter(npmm=25)
      parameter(npjg=34)
      double precision X(9,3)
      double precision C1(3)
      double precision H2(3)
      double precision H3(3)
      double precision H4(3)
      double precision C5(3)
      double precision H6(3)
      double precision H7(3)
      double precision H8(3)
      double precision H9(3)
      double precision eng_kcal
      double precision Vmm3
      double precision Vjg
      double precision etotalmm
      double precision paramm(npmm)
      double precision parajg(npjg)

C average bond length
      double precision avgch
      double precision avgbh
      double precision rxncrd
      double precision r1,r2,r3
      double precision b1,b2,b3

C scale=C-H/C-C 
      double precision scale
C%%%  parameter(scale=0.7153679930d0)
      parameter(scale=0.7171830610d0)
      double precision H5prime(3)

      integer i
      logical LMM3
      parameter(LMM3=.FALSE.)
C--------------------------------------------
C Convert geometry from:
C
C     H2 H6
C     |  |
C  H4-C1-C5-H8.....H9
C     |  |
C     H3 H7
C
C to:
C
C     H6 H3
C     |  |
C  H8-C5-C1-H2.....H9
C     |  |
C     H7 H4
C
      do i=1,3
        C1(i) = X(5,i)
        H2(i) = X(8,i)
        H3(i) = X(6,i)
        H4(i) = X(7,i)
        C5(i) = X(1,i)
        H6(i) = X(2,i)
        H7(i) = X(3,i)
        H8(i) = X(4,i)
        H9(i) = X(9,i)
        H5prime(i) = (scale*C5(i)) + ((1.0d0-scale)*C1(i))
      end do


C Calculate avg C-H bond length for 3 QM H
      r1=0.0d0
      r2=0.0d0
      r3=0.0d0
      b1=0.0d0
      b2=0.0d0
      b3=0.0d0
      do i=1,3
        r1=r1+( (C1(i)-H2(i))**2 )
        r2=r2+( (C1(i)-H3(i))**2 )
        r3=r3+( (C1(i)-H4(i))**2 )
        b1=b1+( (H9(i)-H2(i))**2 )
        b2=b2+( (H9(i)-H3(i))**2 )
        b3=b3+( (H9(i)-H4(i))**2 )
      end do
      r1=dsqrt(r1)
      r2=dsqrt(r2)
      r3=dsqrt(r3)
      b1=dsqrt(b1)
      b2=dsqrt(b2)
      b3=dsqrt(b3)
      avgch=(r1+r2+r3)/3.0d0
      avgbh=(1.0d0/b1)+(1.0d0/b2)+(1.0d0/b3)
      rxncrd=avgch*avgbh
C
C Get energy
C
      if(LMM3) then
        call pmm3(npmm,paramm) 
        call ethane(npmm,paramm,rxncrd,
     1       C1,H2,H3,H4,C5,H6,H7,H8,Vmm3,etotalmm)
        eng_kcal=etotalmm
      else
        call pmm3_mod(npmm,paramm,avgch,rxncrd) 
        call pjg_mod(npjg,parajg,avgch) 
        call ethane(npmm,paramm,rxncrd,
     1       C1,H2,H3,H4,C5,H6,H7,H8,Vmm3,etotalmm)
        call xch5pot(npjg,parajg,rxncrd,
     1       C1,H2,H3,H4,H5prime,H9,Vjg)
        eng_kcal=Vjg+Vmm3
      end if



      return
      END
C end of routine c2h7pot

C**************************************************
      SUBROUTINE pmm3_mod(np,p,ravg,rxncrd) 
C**************************************************
      implicit none
C Output Variables
      integer np
      double precision p(np)
      double precision ravg
      double precision rxncrd
C Local Variables
      double precision r0ch,r0cc
      double precision ksch,kscc
      double precision t0hch,t0cch
      double precision t0hcRh,t0ccRh
      double precision kthch,ktcch
      double precision ksbhch,ksbcch
      double precision kbbhch,kbbcch
      double precision v1hcch,v2hcch,v3hcch
      double precision ksthcch
      double precision esph,espc,rvh,rvc,scalech
      double precision r0cRh,r0ccR
C     parameter(r0ch=1.1120d0,r0cc=1.52470d0)
C##   parameter(r0ch=1.09120d0,r0cc=1.51950d0)

      parameter(ksch=4.740d0,kscc=4.490d0)
C%%%  parameter(t0hch=107.80d0,t0cch=110.70d0)
C%%%  parameter(t0hcRh=107.60d0,t0ccRh=109.310d0)
      parameter(kthch=0.550d0,ktcch=0.590d0)
      parameter(ksbhch=0.0d0,ksbcch=0.080d0)
      parameter(kbbhch=0.0d0,kbbcch=0.30d0**2)
      parameter(v1hcch=0.0d0,v2hcch=0.0d0,v3hcch=0.2380d0)
      parameter(ksthcch=0.0590d0)
      parameter(espc=0.0270d0,rvc=2.040d0)
      parameter(esph=0.020d0 ,rvh=1.620d0)
      parameter(scalech=0.9230d0)
C%%%  parameter(r0cRh=1.1120d0,r0ccR=1.52470d0)
C----------------------------------------
C Switching 
      double precision ravg0
      double precision cr,cp,cs
C--------------------------------------
C Modified parameters
C--------------------------------------
      integer i,j
      integer ipt
      double precision x
      double precision t1,t2,t3
      double precision t4,t5,t6
      double precision z1,z2,z3
      double precision z4,z5,z6
      double precision Ar0ch,Ar0cc
      double precision Asch,Ascc
      double precision Av3hcch
      double precision Aksthcch
      double precision Akthch
      double precision Aktcch
      double precision Aesph
      double precision Arvh
      double precision At0hch
      double precision Ar0cRh
      double precision Ar0ccR
      double precision At0ccRh
C%%%  parameter(r0ch=1.1120d0,r0cc=1.52470d0)
C     parameter(r0ch=1.092510d0,r0cc=1.52720d0)


C%%% for ethane

      IPT = 4

      if(ipt .eq. 1) then
c       ETHANE
        r0ch=1.0850d0
        r0cc=1.510d0
        t0hch=107.52790d0
        t0cch=111.34990d0
        r0cRh=r0ch
        r0ccR=1.47110d0
        t0hcRh=120.0d0
        t0ccRh=120.0d0
      elseif(ipt .eq. 2) then
c       TS
        r0ch=1.0850d0
        r0cc=1.510d0
        t0hch=107.52790d0
        t0cch=111.34990d0
        r0cRh=r0ch
        r0ccR=1.47110d0
        t0hcRh=120.0d0
        t0ccRh=120.0d0
      elseif(ipt .eq. 3) then
c       ETHYL
        r0ch=1.0850d0
        r0cc=1.510d0
        t0hch=107.52790d0
        t0cch=111.34990d0
        r0cRh=r0ch
        r0ccR=1.47110d0
        t0hcRh=120.0d0
        t0ccRh=120.0d0
      elseif(ipt .eq. 4) then
C--------------------------------------
C geometry
        r0ch=1.0850d0
        r0cc=1.50d0
        t0hch=107.52790d0
        t0cch=111.34990d0
        r0cRh=r0ch
        r0ccR=1.47110d0
        t0hcRh=120.0d0
        t0ccRh=120.0d0
C%%%
        t0ccRh=117.2989954203791430d0
        t0ccRh=121.0d0
        t0cch=t0ccRh-7.10d0

        t1=6.0d0
        t2=0.90d0
        t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
        t4=1.0d0-t3
        t0cch=(t0cch*t3)+(111.34990d0*t4)
        t0ccRh=(t0ccRh*t3)+(111.34990d0*t4)

C       t1=6.0d0
C       t2=20.0d0
C       t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
C       t4=1.0d0-t3
C       t0cch=(t0cch*t4)+(120.0d0*t3)
C       t0ccRh=(t0ccRh*t4)+(120.0d0*t3)
C--------------------------------------
C freq and energy
        cs=10.0d0
        ravg0=1.17320d0
      cp=0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr=1.0d0-cp
 
        Ascc=0.70d0
        Aktcch=0.90d0
        Akthch=0.90d0
C
C       z1=10.0d0
C       z2=2.144106410736830d0
C       z3=dexp(-z1*((rxncrd-z2)**2))
C       Ascc=1.0d0-(0.30d0*z3)
C       Aktcch=1.0d0-(0.10d0*z3)
C       Akthch=1.0d0-(0.10d0*z3)

        t1=6.0d0
        t2=0.950d0
        t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
        t4=(1.0d0*(1.0d0-t3))+(0.70d0*t3)
        t5=(1.0d0*(1.0d0-t3))+(0.90d0*t3)
        t6=(1.0d0*(1.0d0-t3))+(0.90d0*t3)

        z1=6.0d0
        z2=10.0d0
        z3=0.50d0*(1.0d0+dtanh(z1*(rxncrd-z2)))

        z4=(t4*(1.0d0-z3))+(1.0d0*z3)
        z5=(t5*(1.0d0-z3))+(1.0d0*z3)
        z6=(t6*(1.0d0-z3))+(1.0d0*z3)

        Ascc   = z4
        Aktcch = z5
        Akthch = z6

      end if
C###
C     Ar0ch=0.980d0
C     Ar0cc=0.9950d0
C     Av3hcch=1.00d0
C     Aksthcch=1.0d0
C     Aesph=1.0d0
C     Arvh=1.0d0
C     At0hch=1.00d0
C     Akthch=1.00d0
C     Aktcch=1.00d0
C     Ar0cRh=1.0d0
C     Ar0ccR=0.980d0

C     At0ccRh=1.100d0

C     Asch=1.10d0

C%%%

      p(1)=r0ch 
      p(2)=r0cc 
      p(3)=ksch 
      p(4)=kscc * Ascc 
      p(5)=t0hch 
      p(6)=t0cch
      p(7)=t0hcRh
      p(8)=t0ccRh 
      p(9)=kthch * Akthch
      p(10)=ktcch * Aktcch
      p(11)=ksbhch
      p(12)=ksbcch
      p(13)=kbbhch
      p(14)=kbbcch
      p(15)=v1hcch
      p(16)=v2hcch
      p(17)=v3hcch 
      p(18)=ksthcch
      p(19)=espc
      p(20)=rvc
      p(21)=esph 
      p(22)=rvh 
      p(23)=scalech 
      p(24)=r0cRh 
      p(25)=r0ccR 


      return
      END   
C end of routine pmm3_mod
C**************************************************
      SUBROUTINE pjg_mod(np,para,ravg)
C**************************************************
      implicit none
C I/O variables
      integer np
      double precision para(np)
      double precision ravg
C Parameter list
      double precision a1ch
      double precision b1ch
      double precision c1ch
      double precision alpcb
      double precision d1cb
      double precision d3cb
      double precision r0cb
      double precision alphh
      double precision d1hh
      double precision d3hh
      double precision r0hh
      double precision alpch
      double precision d1ch
      double precision d3ch
      double precision r0ch
C
      double precision tau
      double precision aphi
      double precision bphi
      double precision cphi
      double precision atheta
      double precision btheta
      double precision ctheta
      double precision a3s
      double precision b3s
      double precision fd0
      double precision hd0
C
      double precision a1s
      double precision a2s
      double precision b1s
      double precision b2s
      double precision ak
      double precision kch3
      double precision aa1
      double precision aa2
      double precision aa3
      double precision aa4
C
      double precision fact1
      double precision fact2
      parameter(fact1=0.041840d0)
      parameter(fact2=6.022045d0)
C streching parameters
      parameter (a1ch=1.7130d0)
      parameter (b1ch=0.1350d0)
      parameter (c1ch=6.6140d0)
      parameter (alpcb=1.8530285d0)
      parameter (d1cb=26.4090d0*fact1)
      parameter (d3cb=20.0630d0*fact1)
C     parameter (r0cb=1.093970d0)
C%%%  parameter (r0cb=1.092510d0)
C%%%  parameter (r0ch=1.092510d0)

C##   parameter (r0cb=1.09030d0)
C##   parameter (r0ch=1.09070d0)

      parameter (alphh=1.94570d0)
      parameter (d1hh=109.4580d0*fact1)
      parameter (d3hh=39.6640d0*fact1)
C%%%  parameter (r0hh=0.741910d0)
C     parameter (alphch=CALC_AT_RUNTIME)
      parameter (d1ch=112.230d0*fact1)
      parameter (d3ch=38.8340d0*fact1)
C     parameter (r0ch=1.093970d0)
C##   parameter (r0ch=1.092510d0)

C out-of-plane bending parameters
C     parameter (tau=CALC_AT_RUNTIME)
      parameter (aphi=0.52879030d0)
      parameter (bphi=0.40066380d0)
      parameter (cphi=1.92099370d0)
      parameter (atheta=0.90787140d0)
      parameter (btheta=0.35488590d0)
      parameter (ctheta=1.89154970d0)
      parameter (a3s=0.14191470d0)
      parameter (b3s=-0.3068450d0)
      parameter (fd0=0.0957500d0*fact2)
      parameter (hd0=0.1915000d0*fact2)
C in-plane bending parameters
      parameter (a1s=1.53136810d-7)
      parameter (a2s=1.01474020d-7)
      parameter (b1s=-4.6696246d0)
      parameter (b2s=-12.363798d0)
      parameter (ak=0.1260d0*fact2)
      parameter (kch3=0.40770d0*fact2)
      parameter (aa1=3.2139520d0)
      parameter (aa2=1.5999630d0)
      parameter (aa3=2.1659530d0)
      parameter (aa4=11.5699770d0)
C----------------------------------------
C Switching 
      double precision ravg0
      double precision cr,cp,cs
C----------------------------------------
C###
      integer i
      integer ipt
      logical LABINITIO
      double precision Ad1hh
      double precision Ad1ch
      double precision Aalpcb
C
C%%% QM parameters

      IPT = 4

      Ad1hh  = 1.0d0
      Ad1ch  = 1.0d0
      Aalpcb = 1.0d0

      if(ipt .eq. 1) then
c       ETHANE
        r0ch=1.0850d0
        r0cb=1.400d0
        r0hh=0.680d0
      elseif(ipt .eq. 2) then
c       TS
        r0ch=1.080d0
        r0cb=1.400d0
        r0hh=0.660d0
      elseif(ipt .eq. 3) then
c       ETHYL
        r0ch=1.0750d0
        r0cb=1.400d0
        r0hh=0.73590d0
      elseif(ipt .eq. 4) then
C------------------------------------------------------
C for geometry
        cs=10.0d0
        ravg0=1.17320d0
      cp=0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr=1.0d0-cp
        r0ch=(1.0850d0*cr)+(1.0750d0*cp)
        r0cb=1.400d0
C       r0hh=0.73590d0
C------------------------------------------------------
C for rxn energy and freq
        cs=10.0d0
        ravg0=1.7320d0
      cp=0.50d0*(1.0d0+(dtanh(cs*(ravg-ravg0))))
        cr=1.0d0-cp
        Ad1hh=(1.0d0*cr)+(1.027420d0*cp)
        Ad1hh=(1.0d0*cr)+(1.027500d0*cp)
        r0hh=(0.720d0*cr)+(0.73590d0*cp)
        r0hh=(0.70d0*cr)+(0.73590d0*cp)
      end if
C======================================================
C Str
      para(1)=a1ch 
      para(2)=b1ch 
      para(3)=c1ch
      para(4)=alpcb 
      para(5)=d1cb
      para(6)=d3cb
      para(7)=r0cb 
      para(8)=alphh 
      para(9)=d1hh * Ad1hh 
      para(10)=d3hh 
      para(11)=r0hh 
      para(12)=d1ch 
      para(13)=d3ch 
      para(14)=r0ch
C out-of-plane
      para(15)=aphi
      para(16)=bphi
      para(17)=cphi
      para(18)=atheta
      para(19)=btheta
      para(20)=ctheta
      para(21)=a3s
      para(22)=b3s
      para(23)=fd0 
      para(24)=hd0
C in-plane bending parameters
      para(25)=a1s
      para(26)=a2s
      para(27)=b1s
      para(28)=b2s
      para(29)=ak
      para(30)=kch3 
      para(31)=aa1
      para(32)=aa2
      para(33)=aa3
      para(34)=aa4

      return
      END
C end of routine pjp_mod
C**************************************************
      SUBROUTINE pmm3(np,p) 
C**************************************************
      implicit none
C Output Variables
      integer np
      double precision p(np)
C Local Variables
      double precision r0ch,r0cc
      double precision ksch,kscc
      double precision t0hch,t0cch
      double precision t0hcRh,t0ccRh
      double precision kthch,ktcch
      double precision ksbhch,ksbcch
      double precision kbbhch,kbbcch
      double precision v1hcch,v2hcch,v3hcch
      double precision ksthcch
      double precision esph,espc,rvh,rvc,scalech
      double precision r0cRh,r0ccR

      parameter(r0ch=1.1120d0,r0cc=1.52470d0)
      parameter(ksch=4.740d0,kscc=4.490d0)
      parameter(t0hch=107.80d0,t0cch=110.70d0)
      parameter(t0hcRh=107.60d0,t0ccRh=109.310d0)
      parameter(kthch=0.550d0,ktcch=0.590d0)
      parameter(ksbhch=0.0d0,ksbcch=0.080d0)
      parameter(kbbhch=0.0d0,kbbcch=0.30d0**2)
      parameter(v1hcch=0.0d0,v2hcch=0.0d0,v3hcch=0.2380d0)
      parameter(ksthcch=0.0590d0)
      parameter(espc=0.0270d0,rvc=2.040d0)
      parameter(esph=0.020d0 ,rvh=1.620d0)
      parameter(scalech=0.9230d0)
      parameter(r0cRh=1.1120d0,r0ccR=1.52470d0)
      p(1)=r0ch
      p(2)=r0cc
      p(3)=ksch
      p(4)=kscc
      p(5)=t0hch
      p(6)=t0cch
      p(7)=t0hcRh
      p(8)=t0ccRh
      p(9)=kthch
      p(10)=ktcch
      p(11)=ksbhch
      p(12)=ksbcch
      p(13)=kbbhch
      p(14)=kbbcch
      p(15)=v1hcch
      p(16)=v2hcch
      p(17)=v3hcch
      p(18)=ksthcch
      p(19)=espc
      p(20)=rvc
      p(21)=esph
      p(22)=rvh
      p(23)=scalech
      p(24)=r0cRh
      p(25)=r0ccR
      return
      END   
C end of routine pmm3
C**************************************************
      SUBROUTINE pjg(np,para)
C**************************************************
      implicit none
C I/O variables
      integer np
      double precision para(np)
C Parameter list
      double precision a1ch
      double precision b1ch
      double precision c1ch
      double precision alpcb
      double precision d1cb
      double precision d3cb
      double precision r0cb
      double precision alphh
      double precision d1hh
      double precision d3hh
      double precision r0hh
      double precision alpch
      double precision d1ch
      double precision d3ch
      double precision r0ch
C
      double precision tau
      double precision aphi
      double precision bphi
      double precision cphi
      double precision atheta
      double precision btheta
      double precision ctheta
      double precision a3s
      double precision b3s
      double precision fd0
      double precision hd0
C
      double precision a1s
      double precision a2s
      double precision b1s
      double precision b2s
      double precision ak
      double precision kch3
      double precision aa1
      double precision aa2
      double precision aa3
      double precision aa4
C
      double precision fact1
      double precision fact2
      parameter(fact1=0.041840d0)
      parameter(fact2=6.022045d0)
C streching parameters
      parameter (a1ch=1.7130d0)
      parameter (b1ch=0.1350d0)
      parameter (c1ch=6.6140d0)
      parameter (alpcb=1.8530285d0)
      parameter (d1cb=26.4090d0*fact1)
      parameter (d3cb=20.0630d0*fact1)
      parameter (r0cb=1.093970d0)
      parameter (alphh=1.94570d0)
      parameter (d1hh=109.4580d0*fact1)
      parameter (d3hh=39.6640d0*fact1)
      parameter (r0hh=0.741910d0)
C     parameter (alphch=CALC_AT_RUNTIME)
      parameter (d1ch=112.230d0*fact1)
      parameter (d3ch=38.8340d0*fact1)
      parameter (r0ch=1.093970d0)
C out-of-plane bending parameters
C     parameter (tau=CALC_AT_RUNTIME)
      parameter (aphi=0.52879030d0)
      parameter (bphi=0.40066380d0)
      parameter (cphi=1.92099370d0)
      parameter (atheta=0.90787140d0)
      parameter (btheta=0.35488590d0)
      parameter (ctheta=1.89154970d0)
      parameter (a3s=0.14191470d0)
      parameter (b3s=-0.3068450d0)
      parameter (fd0=0.0957500d0*fact2)
      parameter (hd0=0.1915000d0*fact2)
C in-plane bending parameters
      parameter (a1s=1.53136810d-7)
      parameter (a2s=1.01474020d-7)
      parameter (b1s=-4.6696246d0)
      parameter (b2s=-12.363798d0)
      parameter (ak=0.1260d0*fact2)
      parameter (kch3=0.40770d0*fact2)
      parameter (aa1=3.2139520d0)
      parameter (aa2=1.5999630d0)
      parameter (aa3=2.1659530d0)
      parameter (aa4=11.5699770d0)

C Str
      para(1)=a1ch
      para(2)=b1ch
      para(3)=c1ch
      para(4)=alpcb
      para(5)=d1cb
      para(6)=d3cb
      para(7)=r0cb
      para(8)=alphh
      para(9)=d1hh
      para(10)=d3hh
      para(11)=r0hh
      para(12)=d1ch
      para(13)=d3ch
      para(14)=r0ch
C out-of-plane
      para(15)=aphi
      para(16)=bphi
      para(17)=cphi
      para(18)=atheta
      para(19)=btheta
      para(20)=ctheta
      para(21)=a3s
      para(22)=b3s
      para(23)=fd0
      para(24)=hd0
C in-plane bending parameters
      para(25)=a1s
      para(26)=a2s
      para(27)=b1s
      para(28)=b2s
      para(29)=ak
      para(30)=kch3
      para(31)=aa1
      para(32)=aa2
      para(33)=aa3
      para(34)=aa4
      return
      END
C end of routine pjp

C**************************************************
      SUBROUTINE ethane(np,para,rxncrd,
     1           C1,H2,H3,H4,C5,H6,H7,H8,Vmm3,etotalmm)
C**************************************************
      implicit none
C Input variables............
      integer np
      double precision para(np)
      double precision rxncrd
      double precision C1(3)
      double precision H2(3)
      double precision H3(3)
      double precision H4(3)
      double precision C5(3)
      double precision H6(3)
      double precision H7(3)
      double precision H8(3)
C Output variables............
      double precision Vmm3
      double precision estr
      double precision ebend
      double precision estrbnd
      double precision ebndbnd
      double precision etor
      double precision estrtor
      double precision evdwaal
      double precision etotalmm

C parameters:
C Angles in Deg.
C Dist in Ang.
      double precision r0ch,r0cc
      double precision ksch,kscc
      double precision t0hch,t0cch
      double precision t0hcRh,t0ccRh
      double precision kthch,ktcch
      double precision ksbhch,ksbcch
      double precision kbbhch,kbbcch
      double precision v1hcch,v2hcch,v3hcch
      double precision ksthcch
      double precision esph,espc,rvh,rvc,scalech


C Func called:
      double precision xdist
      double precision xangle
      double precision xtorsion

C Local parameters:
      integer IOUT
      logical LREACT
      logical LRADICAL
      logical LSWITCH
      logical LJG
      logical LDEBUGMM
      double precision PI
      parameter(IOUT=6)
      parameter(LREACT=.FALSE.)
      parameter(LRADICAL=.FALSE.)
      parameter(LSWITCH=.TRUE.)
      parameter(LJG=.TRUE.)
      parameter(LDEBUGMM=.FALSE.)
      parameter(PI=3.1415926540d0)

C Stretch varaibles
      double precision rch(6),rcc
      double precision esch,escc

C Bend varaibles
      double precision thch(6),tcch(6)
      double precision ebhch,ebcch

C Stretch-Bend varaibles
      double precision esbcch

C Bend-Bend varaibles
      double precision ebbhch,ebbcch

C Torsional variables
      double precision phcch(9)
      double precision ethcch

C Stretch-xtorsion variables
      double precision esthcch

C van der Waal varaibles
      double precision sH2(3)
      double precision sH3(3)
      double precision sH4(3)
      double precision sH5(3)
      double precision sH6(3)
      double precision sH7(3)
      double precision sH8(3)
      double precision rvdwhh(9)
      double precision evdwhh

C energy compponents
      double precision astr(7) 
      double precision abhch(6)
      double precision abcch(6)
      double precision asbcch(6)
      double precision abbhch(6)
      double precision abbcch(6)
      double precision athcch(9)
      double precision asthcch(9)
      double precision avdwhh(9) 
      double precision tstr(7) 
      double precision tbhch(6)
      double precision tbcch(6)
      double precision tsbcch(6)
      double precision tbbhch(6)
      double precision tbbcch(6)
      double precision tthcch(9)
      double precision tsthcch(9)
      double precision tvdwhh(9) 

C Switching func
      double precision t1
      double precision t2
      double precision scale
      double precision sw1(3)
      double precision sw2
      double precision swt1(3)
      double precision swt2
      double precision r0cRh
      double precision r0ccR
      double precision r0chSW
      double precision r0ccSW
      double precision t0hchSW
      double precision t0cchSW
      parameter(t1=22.0d0)
      parameter(t2=1.354030d0)

C Misc varaibles
      integer i,j,k
      double precision xrv,xesp
      double precision xsumRX
      double precision xsumJG
      double precision xstr
      double precision xbend
      double precision xstrbnd
      double precision xbndbnd
      double precision xtor
      double precision xstrtor
      double precision xvdw
      double precision z1
      double precision z2
      double precision z3
      double precision z4
      double precision y1
      double precision y2
      double precision y3
C-------------------------------------------
C Assign parameters
      r0ch=para(1)
      r0cc=para(2)
      ksch=para(3)
      kscc=para(4)
      t0hch=para(5)
      t0cch=para(6)
      t0hcRh=para(7)
      t0ccRh=para(8)
      kthch=para(9)
      ktcch=para(10)
      ksbhch=para(11)
      ksbcch=para(12)
      kbbhch=para(13)
      kbbcch=para(14)
      v1hcch=para(15)
      v2hcch=para(16)
      v3hcch=para(17)
      ksthcch=para(18)
      espc=para(19)
      rvc=para(20)
      esph=para(21)
      rvh=para(22)
      scalech=para(23)
      r0cRh=para(24)
      r0ccR=para(25)
C-------------------------------------------
C     write(IOUT,*)
C     write(IOUT,*)'LREACT   :',LREACT
C     write(IOUT,*)'LRADICAL :',LRADICAL
C-------------------------------------------
C Evaluate bond xdist
      rch(1) = xdist(C1,H2)
      rch(2) = xdist(C1,H3)
      rch(3) = xdist(C1,H4)
      rch(4) = xdist(C5,H6)
      rch(5) = xdist(C5,H7)
      rch(6) = xdist(C5,H8)
      rcc    = xdist(C1,C5)
C-------------------------------------------
C###
C calculate the swithing func
      swt2=1.0d0
      do i=1,3
        swt1(i)=(1.0d0-dtanh(t1*(rch(i)-t2)))*0.50d0
        swt2=swt2*swt1(i)
      end do
C switch C-C xdist
      r0cc=(swt2*r0cc)+((1.0d0-swt2)*r0ccR)
C%%%
C switch force constants
      z1=10.0d0
      z2=2.144106410736830d0
      z3=dexp(-z1*((rxncrd-z2)**2))
      scale=1.0d0-(0.050d0*z3)
C     scale=0.950d0
      
      ksch=ksch*scale
      kthch=kthch*scale
      ktcch=ktcch*scale
C-------------------------------------------
C calculate str energy
      call Es(r0cc,kscc,rcc,escc)   
      estr = escc
      do i=1,6
            call Es(r0ch,ksch,rch(i),esch)
            estr = estr + esch
            astr(i) = esch
      end do
      astr(7) = escc
C-------------------------------------------
C calculate bend energy
C there are 12 bending terms
      thch(1) = xangle(H2,C1,H3)*180.0d0/PI
      thch(2) = xangle(H2,C1,H4)*180.0d0/PI
      thch(3) = xangle(H3,C1,H4)*180.0d0/PI
      thch(4) = xangle(H6,C5,H7)*180.0d0/PI
      thch(5) = xangle(H6,C5,H8)*180.0d0/PI
      thch(6) = xangle(H7,C5,H8)*180.0d0/PI

      tcch(1) = xangle(H2,C1,C5)*180.0d0/PI
      tcch(2) = xangle(H3,C1,C5)*180.0d0/PI
      tcch(3) = xangle(H4,C1,C5)*180.0d0/PI
      tcch(4) = xangle(H6,C5,C1)*180.0d0/PI
      tcch(5) = xangle(H7,C5,C1)*180.0d0/PI
      tcch(6) = xangle(H8,C5,C1)*180.0d0/PI

      ebend = 0.0d0
      do i=1,6
        call Eb(t0hch,kthch,thch(i),ebhch)
        call Eb(t0cch,ktcch,tcch(i),ebcch)
        ebend = ebend + ebhch + ebcch
        abhch(i)= ebhch
        abcch(i)= ebcch
      end do
C-------------------------------------------
C Calculate Stretch-bend terms
C there are 6 stretch-bend cross terms
C H2-C1-C5
C H3-C1-C5
C H4-C1-C5
C H6-C5-C1
C H7-C5-C1
C H8-C5-C1
C
      estrbnd=0.0d0
      do i=1,6
        call Esb(r0ch,r0cc,t0cch,ksbcch,rch(i),rcc,tcch(i),esbcch)
        estrbnd = estrbnd + esbcch
        asbcch(i) = esbcch
      end do
C-------------------------------------------
C Bend-Bend terms
C there are 12 bend-bend terms
C H2-C1-H3 & H2-C1-H4
C H2-C1-H3 & H3-C1-H4
C H2-C1-H4 & H3-C1-H4
C H6-C5-H7 & H6-C5-H8
C H6-C5-H7 & H7-C5-H8
C H6-C5-H8 & H7-C5-H8
C
C H2-C1-C5 & H3-C1-C5
C H2-C1-C5 & H4-C1-C5
C H3-C1-C5 & H4-C1-C5
C H6-C5-C1 & H6-C5-C1
C H6-C5-C1 & H7-C5-C1
C H7-C5-C1 & H8-C5-C1
      
      ebndbnd=0.0d0
      k = 0
      do i=1,3
      do j=i+1,3
        k = k + 1
        call Ebb(t0hch,t0hch,kbbhch,thch(i),thch(j),ebbhch)
        call Ebb(t0cch,t0cch,kbbcch,tcch(i),tcch(j),ebbcch)
        ebndbnd = ebndbnd + ebbhch + ebbcch
        abbhch(k) = ebbhch
        abbcch(k) = ebbcch
      end do
      end do
      do i=4,6
      do j=i+1,6
        k = k + 1
        call Ebb(t0hch,t0hch,kbbhch,thch(i),thch(j),ebbhch)
        call Ebb(t0cch,t0cch,kbbcch,tcch(i),tcch(j),ebbcch)
        ebndbnd = ebndbnd + ebbhch + ebbcch
        abbhch(k) = ebbhch
        abbcch(k) = ebbcch
      end do
      end do
C-------------------------------------------
C there are 9 xtorsional terms
C H2-C1-C5-H6
C H2-C1-C5-H7
C H2-C1-C5-H8
C H3-C1-C5-H6
C H3-C1-C5-H7
C H3-C1-C5-H8
C H4-C1-C5-H6
C H4-C1-C5-H7
C H4-C1-C5-H8
C
C Scale the C-H bonds and evaluate the 
C coordinates of the scaled H
      phcch(1)=xtorsion(H2,C1,C5,H6)
      phcch(2)=xtorsion(H2,C1,C5,H7)
      phcch(3)=xtorsion(H2,C1,C5,H8)
      phcch(4)=xtorsion(H3,C1,C5,H6)
      phcch(5)=xtorsion(H3,C1,C5,H7)
      phcch(6)=xtorsion(H3,C1,C5,H8)
      phcch(7)=xtorsion(H4,C1,C5,H6)
      phcch(8)=xtorsion(H4,C1,C5,H7)
      phcch(9)=xtorsion(H4,C1,C5,H8)

      etor=0.0d0
      do i=1,9
        call Et(v1hcch,v2hcch,v3hcch,phcch(i),ethcch)
        etor = etor + ethcch
        athcch(i) = ethcch
        if(LDEBUGMM) then
          write(*,*) 'Torsion Angle:',phcch(i)
        end if
      end do
C-------------------------------------------
C stretch-xtorsion terms
C interaction of the C-C str with xtorsion
C there are 9 such terms
C H2-C1-C5-H6
C H2-C1-C5-H7
C H2-C1-C5-H8
C H3-C1-C5-H6
C H3-C1-C5-H7
C H3-C1-C5-H8
C H4-C1-C5-H6
C H4-C1-C5-H7
C H4-C1-C5-H8
C
      estrtor=0.0d0
      do i=1,9
        call Est(r0cc,ksthcch,rcc,phcch(i),esthcch)
        estrtor = estrtor + esthcch
        asthcch(i) = esthcch
      end do
C-------------------------------------------
C Van der Waals terms
C all interaction which are greater than 1-3
C are van der Waal intereactions
C there are 9 vdw terms
C H2-H6
C H2-H7
C H2-H8
C H3-H6
C H3-H7
C H3-H8
C H4-H6
C H4-H7
C H4-H8
C
      do i=1,3
        sH2(i) = (scalech*H2(i)) + ((1.0d0-scalech)*C1(i))
        sH3(i) = (scalech*H3(i)) + ((1.0d0-scalech)*C1(i))
        sH4(i) = (scalech*H4(i)) + ((1.0d0-scalech)*C1(i))
        sH6(i) = (scalech*H6(i)) + ((1.0d0-scalech)*C5(i))
        sH7(i) = (scalech*H7(i)) + ((1.0d0-scalech)*C5(i))
        sH8(i) = (scalech*H8(i)) + ((1.0d0-scalech)*C5(i))
      end do
      
      rvdwhh(1) = xdist(sH2,sH6) 
      rvdwhh(2) = xdist(sH2,sH7) 
      rvdwhh(3) = xdist(sH2,sH8) 
      rvdwhh(4) = xdist(sH3,sH6) 
      rvdwhh(5) = xdist(sH3,sH7) 
      rvdwhh(6) = xdist(sH3,sH8) 
      rvdwhh(7) = xdist(sH4,sH6) 
      rvdwhh(8) = xdist(sH4,sH7) 
      rvdwhh(9) = xdist(sH4,sH8) 

      evdwaal=0.0d0
      do i=1,9
        xrv = rvh + rvh
        xesp = dsqrt(esph*esph)
        call Evdw(xrv,xesp,rvdwhh(i),evdwhh)
        evdwaal = evdwaal + evdwhh
        avdwhh(i) = evdwhh
      end do

      etotalmm = estr+ebend+estrbnd+ebndbnd+
     1         etor+estrtor+evdwaal
C-------------------------------------------------------
C Switch between the parameters of ethane and ethyl radical
C swt2=1 for ethane and 0 for ethyl
C###
      r0ccSW=r0cc
      r0chSW=(swt2*r0ch)+((1.0d0-swt2)*r0cRh)
      t0hchSW=(swt2*t0hch)+((1.0d0-swt2)*t0hcRh)
      t0cchSW=(swt2*t0cch)+((1.0d0-swt2)*t0ccRh)
C-----------------------------------------------------------
C%%%
      z2=2.141145010007599490d0
      z2=1.50d0
      z1=0.50d0*(1.0d0+dtanh(6.0d0*(rxncrd-z2)))
      z4=((1.0d0-z1)*t0cch)+(z1*t0ccRh)
      z1=117.2989954203791430d0

      sw2=1.0d0
      do i=1,3
        z1=10.0d0
        sw1(i)=(1.0d0-dtanh(z1*(rch(i)-t2)))*0.50d0
        sw2=sw2*sw1(i)
      end do
      t0cchSW=(sw2*t0cch)+((1.0d0-sw2)*t0ccRh)
C-----------------------------------------------------------

      if(LSWITCH) then
        do i=1,3
          call Eb(t0hchSW,kthch,thch(i),ebhch)
          call Eb(t0cchSW,ktcch,tcch(i),ebcch)
          call Esb(r0chSW,r0ccSW,t0cchSW,ksbcch,
     6             rch(i),rcc,tcch(i),esbcch)
          abhch(i) = ebhch
          abcch(i) = ebcch
          asbcch(i)= esbcch
        end do
        k = 0
        do i=1,3
        do j=i+1,3
          k = k + 1
          call Ebb(t0hchSW,t0hchSW,kbbhch,thch(i),thch(j),ebbhch)
          call Ebb(t0cchSW,t0cchSW,kbbcch,tcch(i),tcch(j),ebbcch)
          abbhch(k) = ebbhch
          abbcch(k) = ebbcch
        end do
        end do
C compute engergies
        do i=1,7
          tstr(i)=astr(i)
        end do
        do i=1,6
          tbhch(i) =abhch(i)
          tbcch(i) =abcch(i)
          tsbcch(i)=asbcch(i)
          tbbhch(i)=abbhch(i)
          tbbcch(i)=abbcch(i)
        end do
        do i=1,9
          tthcch(i)=athcch(i)
          tsthcch(i)=asthcch(i)
          tvdwhh(i)=avdwhh(i)
        end do

C       switch C-H str term
        tstr(1) = tstr(1)*swt1(1)
        tstr(2) = tstr(2)*swt1(2)
        tstr(3) = tstr(3)*swt1(3)

C       switch the following H-C-H bend terms
C       H2-C1-H3
C       H2-C1-H4
C       H3-C1-H4
        tbhch(1) = tbhch(1)*swt1(1)*swt1(2)
        tbhch(2) = tbhch(2)*swt1(1)*swt1(3)
        tbhch(3) = tbhch(3)*swt1(2)*swt1(3)

C       switch the following C-C-H bend terms
C       C5-C1-H2
C       C5-C1-H3
C       C5-C1-H4
        tbcch(1) = tbcch(1)*swt1(1)
        tbcch(2) = tbcch(2)*swt1(2)
        tbcch(3) = tbcch(3)*swt1(3)

C       switch the following C-C-H strbnd term
C       C5-C1-H2
C       C5-C1-H3
C       C5-C1-H4
        tsbcch(1) = tsbcch(1)*swt1(1)
        tsbcch(2) = tsbcch(2)*swt1(2)
        tsbcch(3) = tsbcch(3)*swt1(3)

C       remove three H-C-H2 bnd-bnd terms
        tbbhch(1) = tbbhch(1)*swt1(1)*swt1(2)
        tbbhch(2) = tbbhch(2)*swt1(1)*swt1(3)
        tbbhch(3) = tbbhch(3)*swt1(2)*swt1(3)
        
C       remove two C-C-H2 bnd-bnd terms
        tbbcch(1) = tbbcch(1)*swt1(1)*swt1(2)
        tbbcch(2) = tbbcch(2)*swt1(1)*swt1(3)
        tbbcch(3) = tbbcch(3)*swt1(2)*swt1(3)

C       switch the following  H-C-C-H xtorsional terms
C       H2-C1-C5-H6
C       H2-C1-C5-H7
C       H2-C1-C5-H8
C       H3-C1-C5-H6
C       H3-C1-C5-H7
C       H3-C1-C5-H8
C       H4-C1-C5-H6
C       H4-C1-C5-H7
C       H4-C1-C5-H8
        tthcch(1) = tthcch(1)*swt1(1)
        tthcch(2) = tthcch(2)*swt1(1)
        tthcch(3) = tthcch(3)*swt1(1)
        tthcch(4) = tthcch(4)*swt1(2)
        tthcch(5) = tthcch(5)*swt1(2)
        tthcch(6) = tthcch(6)*swt1(2)
        tthcch(7) = tthcch(7)*swt1(3)
        tthcch(8) = tthcch(8)*swt1(3)
        tthcch(9) = tthcch(9)*swt1(3)
C
C       switch the following H-C-C-H str-xtorsional terms
C       H2-C1-C5-H6
C       H2-C1-C5-H7
C       H2-C1-C5-H8
C       H3-C1-C5-H6
C       H3-C1-C5-H7
C       H3-C1-C5-H8
C       H4-C1-C5-H6
C       H4-C1-C5-H7
C       H4-C1-C5-H8
        tsthcch(1) = tsthcch(1)*swt1(1)
        tsthcch(2) = tsthcch(2)*swt1(1)
        tsthcch(3) = tsthcch(3)*swt1(1)
        tsthcch(4) = tsthcch(4)*swt1(2)
        tsthcch(5) = tsthcch(5)*swt1(2)
        tsthcch(6) = tsthcch(6)*swt1(2)
        tsthcch(7) = tsthcch(7)*swt1(3)
        tsthcch(8) = tsthcch(8)*swt1(3)
        tsthcch(9) = tsthcch(9)*swt1(3)

C       switch the following vdw terms
C       H2-H6
C       H2-H7
C       H2-H8
C       H3-H6
C       H3-H7
C       H3-H8
C       H4-H6
C       H4-H7
C       H4-H8
        tvdwhh(1) = tvdwhh(1)*swt1(1)
        tvdwhh(2) = tvdwhh(2)*swt1(1)
        tvdwhh(3) = tvdwhh(3)*swt1(1)
        tvdwhh(4) = tvdwhh(4)*swt1(2)
        tvdwhh(5) = tvdwhh(5)*swt1(2)
        tvdwhh(6) = tvdwhh(6)*swt1(2)
        tvdwhh(7) = tvdwhh(7)*swt1(3)
        tvdwhh(8) = tvdwhh(8)*swt1(3)
        tvdwhh(9) = tvdwhh(9)*swt1(3)

        xsumRX = 0.0d0
        xstr=0.0d0
        do i=1,7
          xstr=xstr+tstr(i)
        end do
        xbend =0.0d0
        xstrbnd=0.0d0
        xbndbnd=0.0d0
        do i=1,6
          xbend=xbend+tbhch(i)+tbcch(i)
          xstrbnd=xstrbnd+tsbcch(i)
          xbndbnd=xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor=0.0d0
        xstrtor=0.0d0
        xvdw=0.0d0
        do i=1,9
          xtor=xtor+tthcch(i)
          xstrtor=xstrtor+tsthcch(i)
          xvdw=xvdw+tvdwhh(i)
        end do
        xsumRX=xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw

C       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
C    1             xtor,xstrtor,xvdw,xsumRX)

      etotalmm=xsumRX
      end if
C-------------------------------------------------------
      if(LRADICAL) then
C       modify the h-cR-h and h-cR-C bends and c-cR-h str-bnd
C       for C-sp3 bonded to only two atoms
        do i=1,3
          call Eb(t0hcRh,kthch,thch(i),ebhch)
          call Eb(t0ccRh,ktcch,tcch(i),ebcch)
          call Esb(r0ch,r0cc,t0ccRh,ksbcch,rch(i),rcc,tcch(i),esbcch)
          abhch(i) = ebhch
          abcch(i) = ebcch
          asbcch(i)= esbcch
        end do
        k = 0
        do i=1,3
        do j=i+1,3
          k = k + 1
          call Ebb(t0hcRh,t0hcRh,kbbhch,thch(i),thch(j),ebbhch)
          call Ebb(t0ccRh,t0ccRh,kbbcch,tcch(i),tcch(j),ebbcch)
          abbhch(k) = ebbhch
          abbcch(k) = ebbcch
        end do
        end do
      end if

      if(LREACT) then
        do i=1,7
          tstr(i)=astr(i)
        end do
        do i=1,6
          tbhch(i) =abhch(i)
          tbcch(i) =abcch(i)
          tsbcch(i)=asbcch(i)
          tbbhch(i)=abbhch(i)
          tbbcch(i)=abbcch(i)
        end do
        do i=1,9
          tthcch(i)=athcch(i)
          tsthcch(i)=asthcch(i)
          tvdwhh(i)=avdwhh(i)
        end do

C       remove one C-H2 str term
        tstr(1) = 0.0d0

C       remove two H-C-H2 bend terms
        tbhch(1) = 0.0d0
        tbhch(2) = 0.0d0

C       remove one C-C-H2 bend term
        tbcch(1) = 0.0d0

C       remove one C-C-H2 strbnd term
        tsbcch(1) = 0.0d0

C       remove three H-C-H2 bnd-bnd terms
        tbbhch(1) = 0.0d0
        tbbhch(2) = 0.0d0
        tbbhch(3) = 0.0d0
        
C       remove two C-C-H2 bnd-bnd terms
        tbbcch(1) = 0.0d0
        tbbcch(2) = 0.0d0

C       remove three H-C-C-H2 xtorsional terms
        tthcch(1) = 0.0d0
        tthcch(2) = 0.0d0
        tthcch(3) = 0.0d0

C       remove three H-C-C-H2 str-xtorsional terms
        tsthcch(1) = 0.0d0
        tsthcch(2) = 0.0d0
        tsthcch(3) = 0.0d0

C       remove three vdw terms
        tvdwhh(1) = 0.0d0
        tvdwhh(2) = 0.0d0
        tvdwhh(3) = 0.0d0

        xsumRX = 0.0d0
        xstr=0.0d0
        do i=1,7
          xstr=xstr+tstr(i)
        end do
        xbend =0.0d0
        xstrbnd=0.0d0
        xbndbnd=0.0d0
        do i=1,6
          xbend=xbend+tbhch(i)+tbcch(i)
          xstrbnd=xstrbnd+tsbcch(i)
          xbndbnd=xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor=0.0d0
        xstrtor=0.0d0
        xvdw=0.0d0
        do i=1,9
          xtor=xtor+tthcch(i)
          xstrtor=xstrtor+tsthcch(i)
          xvdw=xvdw+tvdwhh(i)
        end do
        xsumRX=xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw

C       call xdisp(xstr,xbend,xstrbnd,xbndbnd,
C    1             xtor,xstrtor,xvdw,xsumRX)

      etotalmm=xsumRX
      end if 
            
        
C start of if-block label: jg_case
C....................................
      if(LJG) then
        do i=1,7
          tstr(i)=astr(i)
        end do
        do i=1,6
          tbhch(i) =abhch(i)
          tbcch(i) =abcch(i)
          tsbcch(i)=asbcch(i)
          tbbhch(i)=abbhch(i)
          tbbcch(i)=abbcch(i)
        end do
        do i=1,9
          tthcch(i)=athcch(i)
          tsthcch(i)=asthcch(i)
          tvdwhh(i)=avdwhh(i)
        end do
 
C       remove the following C-H str terms
C       C-H2
C       C-H3
C       C-H4
        tstr(1) = 0.0d0
        tstr(2) = 0.0d0
        tstr(3) = 0.0d0

C       remove the following H-C-H bend terms
C       H2-C1-H3
C       H2-C1-H4
C       H3-C1-H4
        tbhch(1) = 0.0d0
        tbhch(2) = 0.0d0
        tbhch(3) = 0.0d0

C       switch the following C-C-H bend terms
C       C5-C1-H2
C       C5-C1-H3
C       C5-C1-H4
        tbcch(1) = tbcch(1)*swt1(1)
        tbcch(2) = tbcch(2)*swt1(2)
        tbcch(3) = tbcch(3)*swt1(3)

C       switch the following C-C-H strbnd term
C       C5-C1-H2
C       C5-C1-H3
C       C5-C1-H4
        tsbcch(1) = tsbcch(1)*swt1(1)
        tsbcch(2) = tsbcch(2)*swt1(2)
        tsbcch(3) = tsbcch(3)*swt1(3)

C       remove three H-C-H2 bnd-bnd terms
        tbbhch(1) = 0.0d0
        tbbhch(2) = 0.0d0
        tbbhch(3) = 0.0d0

C       switch C-C-H bnd-bnd terms involving H2, H3 and H4
        tbbcch(1) = tbbcch(1)*swt1(1)
        tbbcch(2) = tbbcch(2)*swt1(2)
        tbbcch(3) = tbbcch(3)*swt1(3)

C       switch the following  H-C-C-H xtorsional terms
C       H2-C1-C5-H6
C       H2-C1-C5-H7
C       H2-C1-C5-H8
C       H3-C1-C5-H6
C       H3-C1-C5-H7
C       H3-C1-C5-H8
C       H4-C1-C5-H6
C       H4-C1-C5-H7
C       H4-C1-C5-H8
        tthcch(1) = tthcch(1)*swt1(1)
        tthcch(2) = tthcch(2)*swt1(1)
        tthcch(3) = tthcch(3)*swt1(1)
        tthcch(4) = tthcch(4)*swt1(2)
        tthcch(5) = tthcch(5)*swt1(2)
        tthcch(6) = tthcch(6)*swt1(2)
        tthcch(7) = tthcch(7)*swt1(3)
        tthcch(8) = tthcch(8)*swt1(3)
        tthcch(9) = tthcch(9)*swt1(3)
C
C       switch the following H-C-C-H str-xtorsional terms
C       H2-C1-C5-H6
C       H2-C1-C5-H7
C       H2-C1-C5-H8
C       H3-C1-C5-H6
C       H3-C1-C5-H7
C       H3-C1-C5-H8
C       H4-C1-C5-H6
C       H4-C1-C5-H7
C       H4-C1-C5-H8
        tsthcch(1) = tsthcch(1)*swt1(1)
        tsthcch(2) = tsthcch(2)*swt1(1)
        tsthcch(3) = tsthcch(3)*swt1(1)
        tsthcch(4) = tsthcch(4)*swt1(2)
        tsthcch(5) = tsthcch(5)*swt1(2)
        tsthcch(6) = tsthcch(6)*swt1(2)
        tsthcch(7) = tsthcch(7)*swt1(3)
        tsthcch(8) = tsthcch(8)*swt1(3)
        tsthcch(9) = tsthcch(9)*swt1(3)
 
C       switch the following vdw terms
C       H2-H6
C       H2-H7
C       H2-H8
C       H3-H6
C       H3-H7
C       H3-H8
C       H4-H6
C       H4-H7
C       H4-H8
        tvdwhh(1) = tvdwhh(1)*swt1(1)
        tvdwhh(2) = tvdwhh(2)*swt1(1)
        tvdwhh(3) = tvdwhh(3)*swt1(1)
        tvdwhh(4) = tvdwhh(4)*swt1(2)
        tvdwhh(5) = tvdwhh(5)*swt1(2)
        tvdwhh(6) = tvdwhh(6)*swt1(2)
        tvdwhh(7) = tvdwhh(7)*swt1(3)
        tvdwhh(8) = tvdwhh(8)*swt1(3)
        tvdwhh(9) = tvdwhh(9)*swt1(3)

        xsumJG = 0.0d0
        xstr=0.0d0
        do i=1,7
          xstr=xstr+tstr(i)
        end do
        xbend =0.0d0
        xstrbnd=0.0d0
        xbndbnd=0.0d0
        do i=1,6
          xbend=xbend+tbhch(i)+tbcch(i)
          xstrbnd=xstrbnd+tsbcch(i)
          xbndbnd=xbndbnd+tbbhch(i)+tbbcch(i)
        end do
        xtor=0.0d0
        xstrtor=0.0d0
        xvdw=0.0d0
        do i=1,9
          xtor=xtor+tthcch(i)
          xstrtor=xstrtor+tsthcch(i)
          xvdw=xvdw+tvdwhh(i)
        end do
        xsumJG=xstr+xbend+xstrbnd+xbndbnd+xtor+xstrtor+xvdw
        if(LDEBUGMM) then
          write(*,1000)
          write(*,1000) '----------------------------'
          write(*,1000) 'Energy from MM3'
          write(*,1000) '----------------------------'
          call xdisp(xstr,xbend,xstrbnd,xbndbnd,
     1             xtor,xstrtor,xvdw,xsumJG)

        end if

      end if 
C end of if-block label: jg_case
C....................................

      Vmm3=xsumJG

1000  format(A,F14.8)
      return
      END 
C end sub ethane
C**************************************************
      SUBROUTINE xdisp(xstr,xbend,xstrbnd,xbndbnd,
     1 xtor,xstrtor,xvdw,xsum)
C**************************************************
      implicit none
      double precision xstr
      double precision xbend
      double precision xstrbnd
      double precision xbndbnd
      double precision xtor
      double precision xstrtor
      double precision xvdw
      double precision xsum
      write(*,1000)
      write(*,1000) 'Energy in kcal/mol'
      write(*,1000)'Stretch   :',xstr
      write(*,1000)'Bend      :',xbend
      write(*,1000)'Str-Bnd   :',xstrbnd
      write(*,1000)'Bnd-Bnd   :',xbndbnd
      write(*,1000)'Torsional :',xtor
      write(*,1000)'Str-Tor   :',xstrtor
      write(*,1000)'vdW       :',xvdw
      write(*,1000)
      write(*,1000)'Total MM eng.  :',xsum
1000  format(A,F14.4)
      return
      END
C end sub xdisp
C**************************************************
      FUNCTION xdist(x,y)
C**************************************************
      implicit none
      double precision x(3)
      double precision y(3)
      double precision xdist
      integer i
      double precision z
      z = 0.0d0
      do i=1,3
        z = z + ( (y(i)-x(i))**2 )
      end do
      xdist = dsqrt(z)

      return
      END
C end of func xdist
C**************************************************
      FUNCTION xangle(a,b,c)
C**************************************************
      implicit none
      double precision a(3)
      double precision b(3)
      double precision c(3)
      double precision xangle
      double precision xdist

      integer i
      double precision sum,rab(3),rcb(3)

      do i=1,3
        rab(i) = a(i)-b(i)
        rcb(i) = c(i)-b(i)
      end do
      sum = 0.0d0
      do i=1,3
        sum = sum + (rab(i)*rcb(i))
      end do
      sum = sum/(xdist(a,b)*xdist(c,b))
      xangle = dabs(dacos(sum))
            
      return
      END
C end of xangle

C***********************************************************
      FUNCTION xtorsion(a,b,c,d) 
C***********************************************************
      implicit none
      double precision a(3)
      double precision b(3)
      double precision c(3)
      double precision d(3)
      double precision xtorsion
C
C This func calculates the dihedral xangle
C between a-b and c-d
C a
C |
C b---c
C     |
C     d
C we will use the cross product definition of the 
C dihedral xangle
C
      integer i,j
      double precision x,y,z
      double precision s
      double precision phi
      double precision dpuv
      double precision rab(3)
      double precision rcb(3)
      double precision rcd(3)
      double precision u(3)
      double precision v(3)
C
C define bond-vectors
      do i=1,3
        rab(i) = a(i)-b(i)
        rcb(i) = c(i)-b(i)
        rcd(i) = c(i)-d(i)
      end do
C
C define the normal vectors 
C u = rab X rcb
C v = rcb X rcd
      call crossprod(rab,rcb,u)
      call crossprod(rcb,rcd,v)
C
C evaluate the norm and dot product 
C u and v 
C dpuv = u . v
C s = rab . v and is used for sign of
C the dihedral xangle
C
      x=0.0d0
      y=0.0d0
      s = 0.0d0
      dpuv = 0.0d0
      do i=1,3
        x = x+(u(i)*u(i))
        y = y+(v(i)*v(i))
        dpuv = dpuv + (u(i)*v(i))
        s = s + (rab(i)*v(i))
      end do
      x = sqrt(x)
      y = sqrt(y)
C
C z is cos(phi) and should be in the
C range [-1,1]
C
      if(x*y .ne. 0.0d0) then
        z = dpuv/(x*y)
      else
        z = dpuv
      end if
      if( z .gt.  1.0d0) z=1.0d0
      if( z .lt. -1.0d0) z=-1.0d0
C get dihedral xangle
      phi = dacos(z)
C
C get the sign
      if(s /= 0.0d0) then
        s = s/dabs(s)
      else
        s = 1.0d0
      end if

      xtorsion = s * phi
      
      return
      END
C end of sub xtorsion

C**************************************************
      SUBROUTINE crossprod(a,b,c)
C**************************************************
      implicit none
      double precision a(3)
      double precision b(3)
      double precision c(3)
      c(1) = (a(2)*b(3))-(a(3)*b(2))
      c(2) = (a(3)*b(1))-(a(1)*b(3))
      c(3) = (a(1)*b(2))-(a(2)*b(1))
      return
      END 

C23456789
C*****************************************
      SUBROUTINE Es(r0,ks,r,eng)
C*****************************************
      implicit none
      double precision r0
      double precision ks
      double precision r
      double precision eng

      double precision dr
      double precision c1,c2,c3
      parameter (c1=71.940d0)
      parameter (c2=2.550d0)
      parameter (c3=7.0d0/12.0d0)

      eng = 0.0d0
      dr  = (r-r0)
      eng = (c1*ks*dr*dr) *
     1 (1.0d0-(c2*dr)+(c3*c2*c2*dr*dr))
C    1 (1.0d0-(c2*dr)+(c3*c2*dr*dr))
      
      return
      END 
C end of sub Es
C*****************************************
      SUBROUTINE Eb(t0,kt,t,eng)
C*****************************************
      implicit none
      double precision t0
      double precision kt
      double precision t
      double precision eng

      double precision dt
      double precision c1,c2,c3,c4,c5
C Alligers (JACS, vol 111, pg 8551) parameters
C     parameter(c1=0.0219140d0)
C     parameter(c2=0.0140d0)
C     parameter(c3=5.60d-5)
C     parameter(c4=7.0d-7)
C     parameter(c5=9.0d-10)
C Tinker 4.1 MM3 parameter 
      parameter(c1=0.021914180d0)
      parameter(c2=0.0140d0)
      parameter(c3=0.0000560d0)
      parameter(c4=0.00000070d0)
      parameter(c5=0.0000000220d0)

      dt = t-t0
      eng = c1*kt*dt*dt*
     1 ( 1.0d0-(c2*dt)+(c3*dt*dt)-(c4*dt*dt*dt)
     2 +(c5*dt*dt*dt*dt) )     
      
      return
      END
C end of sub Eb
C*****************************************
      SUBROUTINE Et(v1,v2,v3,w,eng)
C*****************************************
      implicit none
      double precision v1
      double precision v2
      double precision v3
      double precision w
      double precision eng

      eng = 0.50d0*(
     1 ( v1*(1.0d0+dcos(w)) ) +
     2 ( v2*(1.0d0-dcos(2.0d0*w)) ) +
     3 ( v3*(1.0d0+dcos(3.0d0*w)) ) )


      return
      END 
C end of Et
C*****************************************
      SUBROUTINE Esb(ra0,rb0,t0,kst,ra,rb,t,eng)
C*****************************************
      implicit none
      double precision ra0
      double precision rb0
      double precision t0
      double precision kst
      double precision ra
      double precision rb
      double precision t
      double precision eng

      double precision c1
      parameter(c1=2.511180d0)

      eng = c1*kst*(t-t0)*((ra-ra0)+(rb-rb0))

      END 
C end of Esb
C*****************************************
      SUBROUTINE Est(r0,kws,r,w,eng)
C*****************************************
      implicit none
      double precision r0
      double precision kws
      double precision r
      double precision w
      double precision eng

      double precision c1
      parameter(c1=-11.9950d0)
C In Alligers paper(JACS, vol 111, pg 8551) 
C c1 is +ve but TINKER (version 4.1) defines
C c1/2 as -ve

      eng = 0.50d0*c1*kws*(r-r0)*(1.0d0+dcos(3.0d0*w))

      return
      END
C end of Est
C*****************************************
      SUBROUTINE Ebb(ta0,tb0,kbb,ta,tb,eng)
C*****************************************
      implicit none
      double precision ta0
      double precision tb0
      double precision kbb
      double precision ta
      double precision tb
      double precision eng

      double precision c1
      parameter(c1=-0.021914180d0)

      eng = c1*kbb*(ta-ta0)*(tb-tb0)
C     write(*,1000) c1,kbb,(ta-ta0)*(tb-tb0),eng
1000  format(5F12.5)
      return
      END 
C end of Ebb
C*****************************************
      SUBROUTINE Evdw(rv,esp,r,eng)
C*****************************************
      implicit none
      double precision rv
      double precision esp
      double precision r
      double precision eng

      double precision x
      double precision c1,c2,c3
      parameter(c1=2.250d0)
      parameter(c2=1.840d+5)
      parameter(c3=12.00d0)

      x = rv/r
      eng = esp*((-c1*x*x*x*x*x*x)+(c2*dexp(-c3/x)))

      return
      END
C end of vdw

C23456789
C References:
C [1]. Joseph et al.   : JCP, vol 87(12),  pg. 7036, 1987
C [2]. Jordan et al.   : JCP, vol 102(14), pg. 5669, 1995
C [3]. Steckler et al. : JCP, vol 87(15),  pg. 7024, 1987
C

C*****************************************************************
      SUBROUTINE xch5pot(np,para,rxncrd,C1,H2,H3,H4,H5,Hb,xpot_kcal)
C*****************************************************************
C
C Input : Coordintes are in Ang.
C Output: xpot_kcal is in kcal/mol
C
C All energy variables local to this
C routine are in 10(5)Joules/mol 
C
      implicit none
C Input variables
      integer np
      double precision para(np)
      double precision rxncrd
      double precision C1(3)
      double precision H2(3)
      double precision H3(3)
      double precision H4(3)
      double precision H5(3)
      double precision Hb(3)
C Output variables
      double precision xpot_kcal
C Parameter list
      double precision a1ch
      double precision b1ch
      double precision c1ch
      double precision alpcb
      double precision d1cb
      double precision d3cb
      double precision r0cb
      double precision alphh
      double precision d1hh
      double precision d3hh
      double precision r0hh
      double precision alpch
      double precision d1ch
      double precision d3ch
      double precision r0ch
C
      double precision tau
      double precision aphi
      double precision bphi
      double precision cphi
      double precision atheta
      double precision btheta
      double precision ctheta
      double precision a3s
      double precision b3s
      double precision fd0
      double precision hd0
C
      double precision a1s
      double precision a2s
      double precision b1s
      double precision b2s
      double precision ak
      double precision kch3
      double precision aa1
      double precision aa2
      double precision aa3
      double precision aa4
C
      double precision fact1
      double precision fact2
      parameter(fact1=0.041840d0)
      parameter(fact2=6.022045d0)
C Funct called
      double precision xangle
C Local variables
C     Local parameters
      integer IOUT
      logical LDEBUG
      double precision KJ_IN_KCAL
      parameter(IOUT=6)
      parameter(LDEBUG=.FALSE.)
      parameter(KJ_IN_KCAL=0.2390057)
C
C     bond lenght and bond vectors
C
      double precision rch2
      double precision rch3
      double precision rch4
      double precision rch5
      double precision rbh2
      double precision rbh3
      double precision rbh4
      double precision rbh5
      double precision rcb
      double precision rchav
      double precision bvch2(3)
      double precision bvch3(3)
      double precision bvch4(3)
      double precision bvch5(3)
      double precision bvcb(3)
C
C     Vstr varaibles
C
      double precision cbq
      double precision cbj
      double precision bh2q
      double precision bh2j
      double precision bh3q
      double precision bh3j
      double precision bh4q
      double precision bh4f
      double precision bh5q
      double precision bh5j
      double precision ch2q
      double precision ch2j
      double precision ch3q
      double precision ch3j
      double precision ch4q
      double precision ch4j
      double precision ch5q
      double precision ch5j
      double precision vqh2
      double precision vqh3
      double precision vqh4
      double precision vqh5
      double precision vjh2
      double precision vjh3
      double precision vjh4
      double precision vjh5
      double precision v3h2
      double precision v3h3
      double precision v3h4
      double precision v3h5
      double precision vstr
C
C     Bond xangles and out-of-plane-bendings
C
      double precision PI
      double precision theta023
      double precision theta024
      double precision theta025
      double precision theta032
      double precision theta034
      double precision theta035
      double precision theta042
      double precision theta043
      double precision theta045
      double precision theta052
      double precision theta053
      double precision theta054
      double precision delta23
      double precision delta24
      double precision delta25
      double precision delta32
      double precision delta34
      double precision delta35
      double precision delta42
      double precision delta43
      double precision delta45
      double precision delta52
      double precision delta53
      double precision delta54
      double precision d23
      double precision d24
      double precision d25
      double precision d32
      double precision d34
      double precision d35
      double precision d42
      double precision d43
      double precision d45
      double precision d52
      double precision d53
      double precision d54
      double precision fd2
      double precision fd3
      double precision fd4
      double precision fd5
      double precision hd2
      double precision hd3
      double precision hd4
      double precision hd5
      double precision voph2
      double precision voph3
      double precision voph4
      double precision voph5
      double precision vop
C
C     in-plane bending parameters     
C
      double precision xk023
      double precision xk024
      double precision xk025
      double precision xk034
      double precision xk035
      double precision xk045
      double precision xk2
      double precision xk3
      double precision xk4
      double precision xk5
      double precision theta23
      double precision theta24
      double precision theta25
      double precision theta34
      double precision theta35
      double precision theta45
      double precision vip23
      double precision vip24
      double precision vip25
      double precision vip34
      double precision vip35
      double precision vip45
      double precision vip
C
C     misc and temp variables
C
      integer i,j,k
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision t6
      double precision z1
      double precision z2
      double precision z3
      double precision z4
      double precision scat0
      double precision xpot
C
C     CH5 surf
C
      double precision vsurf
      double precision R(18)
      double precision dv(18)
C
C Assign parameters
C
      a1ch=para(1)
      b1ch=para(2)
      c1ch=para(3)
      alpcb=para(4)
      d1cb=para(5)
      d3cb=para(6)
      r0cb=para(7)
      alphh=para(8)
      d1hh=para(9)
      d3hh=para(10)
      r0hh=para(11)
      d1ch=para(12)
      d3ch=para(13)
      r0ch=para(14)
C out-of-plane
      aphi=para(15)
      bphi=para(16)
      cphi=para(17)
      atheta=para(18)
      btheta=para(19)
      ctheta=para(20)
      a3s=para(21)
      b3s=para(22)
      fd0=para(23)
      hd0=para(24)
C in-plane bending parameters
      a1s=para(25)
      a2s=para(26)
      b1s=para(27)
      b2s=para(28)
      ak=para(29)
      kch3=para(30)
      aa1=para(31)
      aa2=para(32)
      aa3=para(33)
      aa4=para(34)
C
C Define j,k,l indices wrt to each of the 4 H
C H2 : 3,4,5
C H3 : 5,4,2
C H4 : 5,2,3
C H5 : 2,4,3
C
C Calculate bond lenghts and bond vectors
C bv is the bond vectors and rch is the 
C corresponding bond lenghts
C rchav is the avg of the four ch xdist
C
      rch2=0.0d0
      rch3=0.0d0
      rch4=0.0d0
      rch5=0.0d0
      rbh2=0.0d0
      rbh3=0.0d0
      rbh4=0.0d0
      rbh5=0.0d0
      rcb =0.0d0
      do i=1,3
        bvch2(i)=H2(i)-C1(i)
        bvch3(i)=H3(i)-C1(i)
        bvch4(i)=H4(i)-C1(i)
        bvch5(i)=H5(i)-C1(i)
        bvcb(i) =Hb(i)-C1(i)
        rch2    =rch2+(bvch2(i)*bvch2(i))
        rch3    =rch3+(bvch3(i)*bvch3(i))
        rch4    =rch4+(bvch4(i)*bvch4(i))
        rch5    =rch5+(bvch5(i)*bvch5(i))
        rcb     =rcb +( bvcb(i)*bvcb(i) )
        rbh2    =rbh2+((Hb(i)-H2(i))*(Hb(i)-H2(i)))
        rbh3    =rbh3+((Hb(i)-H3(i))*(Hb(i)-H3(i)))
        rbh4    =rbh4+((Hb(i)-H4(i))*(Hb(i)-H4(i)))
        rbh5    =rbh5+((Hb(i)-H5(i))*(Hb(i)-H5(i)))
      end do
      rch2 =dsqrt(rch2)
      rch3 =dsqrt(rch3)
      rch4 =dsqrt(rch4)
      rch5 =dsqrt(rch5)
      rcb  =dsqrt(rcb)
      rbh2 =dsqrt(rbh2)
      rbh3 =dsqrt(rbh3)
      rbh4 =dsqrt(rbh4)
      rbh5 =dsqrt(rbh5)
      rchav=0.250d0*(rch2+rch3+rch4+rch5)
C
C use switch func to change r0ch from
C ethane to ethyl value
C
C####
      t1=1.09250d0
      t2=1.08850d0
      t3=c1ch*10.0d0
      t4=1.12750d0
      t5=0.50d0*(1.0d0+dtanh(t3*(rchav-t4)))
C     r0ch=(t5*t2)+((1.0d0-t5)*t1)

C
C calculate the Morse parameter alphaCH (alpch)
C used in the Vstr using Equation 2 of
C Ref. [2]
C
      t1=0.50d0*( 1.0d0+dtanh(c1ch*(rchav-r0ch)) )
      alpch=a1ch+(b1ch*t1)
C
C calculate Vstr
C first get the coulomb and resonance integrals
C involving Hb and all the other atoms
C cbq and cbj : coulomb and resonance integral between C1 and Hb
C
      call V3qj77(alpcb,d1cb,d3cb,r0cb,rcb,cbq,cbj)
      call V3qj77(alphh,d1hh,d3hh,r0hh,rbh2,bh2q,bh2j)
      call V3qj77(alphh,d1hh,d3hh,r0hh,rbh3,bh3q,bh3j)
      call V3qj77(alphh,d1hh,d3hh,r0hh,rbh4,bh4q,bh4f)
      call V3qj77(alphh,d1hh,d3hh,r0hh,rbh5,bh5q,bh5j)
C%%%--------------------------------------------------
C     t1 = 0.9950d0
C     bh2q=bh2q*t1
C     bh3q=bh3q*t1
C     bh4q=bh4q*t1
C     bh5q=bh5q*t1

C     bh2j=bh2j*t1
C     bh3j=bh3j*t1
C     bh4f=bh4f*t1
C     bh5j=bh5j*t1
C%%%--------------------------------------------------
C
C C1 and (H2...H5) terms
C
C####
      t1=0.50d0*(1.0d0+dtanh(10.0d0*(rbh2-1.50d0)))
      t2=1.0d0-t1
      t3=(1.0d0+(t2*0.080d0))*r0ch
C     call V3qj77(alpch,d1ch,d3ch,t3,rch2,ch2q,ch2j)

      call V3qj77(alpch,d1ch,d3ch,r0ch,rch2,ch2q,ch2j)
      call V3qj77(alpch,d1ch,d3ch,r0ch,rch3,ch3q,ch3j)
      call V3qj77(alpch,d1ch,d3ch,r0ch,rch4,ch4q,ch4j)
      call V3qj77(alpch,d1ch,d3ch,r0ch,rch5,ch5q,ch5j)
C
C coulomb and resonance integral for H2...H5
C Equation 2 of Ref. [1]
C
      vqh2=cbq+ch2q+bh2q
      vqh3=cbq+ch3q+bh3q
      vqh4=cbq+ch4q+bh4q
      vqh5=cbq+ch5q+bh5q


      t2=((ch2j-cbj)**2) +
     1   ((cbj-bh2j )**2) +
     2   ((bh2j-ch2j)**2)

      t3=((ch3j-cbj)**2) +
     1   ((cbj-bh3j )**2) +
     2   ((bh3j-ch3j)**2) 

      t4=((ch4j-cbj)**2) +
     1   ((cbj-bh4f )**2) +
     2   ((bh4f-ch4j)**2) 

      t5=((ch5j-cbj)**2) +
     1   ((cbj-bh5j )**2) +
     2   ((bh5j-ch5j)**2)

      vjh2=-dsqrt(0.50d0*t2)
      vjh3=-dsqrt(0.50d0*t3)
      vjh4=-dsqrt(0.50d0*t4)
      vjh5=-dsqrt(0.50d0*t5)

      v3h2=vqh2+vjh2
      v3h3=vqh3+vjh3
      v3h4=vqh4+vjh4
      v3h5=vqh5+vjh5
      vstr=v3h2+v3h3+v3h4+v3h5

C%%%--------------------------------------------------
C     t1=10.0d0
C     t2=2.144106410736830d0
C     t3=0.013850d0
C     t4=1.0d0+(t3*dexp(-t1*((rxncrd-t2)**2)))

      t1=6.0d0
      t2=0.950d0
      t3=0.50d0*(1.0d0+dtanh(t1*(rxncrd-t2)))
      t4=1.013850d0
      
      call xgaussfit(rxncrd,t4)
C     t4=t4*0.999930d0
C     t4=t4*1.00205d0
C     t4=t4*1.0159995120d0/1.016043066097427340d0
      t4=t4*1.016043066097427340d0/1.0159995120d0
  
      t5=(1.0d0*(1.0d0-t3))+(t3*t4)

      z1=6.0d0
      z2=10.0d0
      z3=0.50d0*(1.0d0+dtanh(z1*(rxncrd-z2)))
      z4=(t5*(1.0d0-z3))+(1.0d0*z3)

      vstr = vstr * z4

C     vstr = vstr * 1.013850d0
C%%%--------------------------------------------------


      if(LDEBUG) then
        write(*,*) 
        write(*,*)'H2 Vstr :',v3h2*1.0d2*KJ_IN_KCAL 
        write(*,*)'H3 Vstr :',v3h3*1.0d2*KJ_IN_KCAL 
        write(*,*)'H4 Vstr :',v3h4*1.0d2*KJ_IN_KCAL 
        write(*,*)'H5 Vstr :',v3h5*1.0d2*KJ_IN_KCAL 
        write(*,*) 'VSTR FROM XPOT:',vstr*1.0d2*KJ_IN_KCAL
      end if
C
C calculate the out-of-plane bend terms
C
C evaluate tau and PI
      PI =dacos(-1.0d0)
      tau=dacos(-1.0d0/3.0d0) 
C
C calculate theta0 
C
      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch2,rch3,rch4,rch5,theta023)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch2,rch4,rch5,rch3,theta024)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch2,rch5,rch3,rch4,theta025)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch3,rch2,rch5,rch4,theta032)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch3,rch4,rch2,rch5,theta034)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch3,rch5,rch4,rch2,theta035)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch4,rch2,rch3,rch5,theta042)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch4,rch3,rch5,rch2,theta043)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch4,rch5,rch2,rch3,theta045)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch5,rch2,rch4,rch3,theta052)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch5,rch3,rch2,rch4,theta053)

      call gjtheta0(PI,tau,aphi,bphi,cphi,
     1     atheta,btheta,ctheta,r0ch,
     2     rch5,rch4,rch3,rch2,theta054)

C     scat0=1.0d0
C     theta023=theta023*scat0
C     theta024=theta024*scat0
C     theta025=theta025*scat0
C     theta032=theta032*scat0
C     theta034=theta034*scat0
C     theta035=theta035*scat0
C     theta042=theta042*scat0
C     theta043=theta043*scat0
C     theta045=theta045*scat0
C     theta052=theta052*scat0
C     theta053=theta053*scat0
C     theta054=theta054*scat0


C
C calculate delta
C
      call steckdel(bvch3,bvch4,bvch5,bvch3,d23)
      call steckdel(bvch3,bvch4,bvch5,bvch4,d24)
      call steckdel(bvch3,bvch4,bvch5,bvch5,d25)

      call steckdel(bvch5,bvch4,bvch2,bvch2,d32)
      call steckdel(bvch5,bvch4,bvch2,bvch4,d34)
      call steckdel(bvch5,bvch4,bvch2,bvch5,d35)

      call steckdel(bvch5,bvch2,bvch3,bvch2,d42)
      call steckdel(bvch5,bvch2,bvch3,bvch3,d43)
      call steckdel(bvch5,bvch2,bvch3,bvch5,d45)

      call steckdel(bvch2,bvch4,bvch3,bvch2,d52)
      call steckdel(bvch2,bvch4,bvch3,bvch3,d53)
      call steckdel(bvch2,bvch4,bvch3,bvch4,d54)

      delta23=d23-theta023
      delta24=d24-theta024
      delta25=d25-theta025

      delta32=d32-theta032
      delta34=d34-theta034
      delta35=d35-theta035

      delta42=d42-theta042
      delta43=d43-theta043
      delta45=d45-theta045

      delta52=d52-theta052
      delta53=d53-theta053
      delta54=d54-theta054

      if(LDEBUG) then
        write(IOUT,*) 
        write(IOUT,*) 'theta0 and delta terms:'
        write(IOUT,*) 'H2-H3:',theta023,delta23
        write(IOUT,*) 'H2-H4:',theta024,delta24
        write(IOUT,*) 'H2-H5:',theta025,delta25
        write(IOUT,*) 'H3-H4:',theta034,delta34
        write(IOUT,*) 'H3-H5:',theta035,delta35
        write(IOUT,*) 'H4-H5:',theta045,delta45
      end if
C
C calculate fd and hd terms
      call calcfh(a3s,b3s,fd0,hd0,r0ch,rch2,rch3,rch4,rch5,fd2,hd2)
      call calcfh(a3s,b3s,fd0,hd0,r0ch,rch3,rch5,rch4,rch2,fd3,hd3)
      call calcfh(a3s,b3s,fd0,hd0,r0ch,rch4,rch5,rch2,rch3,fd4,hd4)
      call calcfh(a3s,b3s,fd0,hd0,r0ch,rch5,rch2,rch4,rch3,fd5,hd5)

      if(LDEBUG) then
        write(IOUT,*)
        write(IOUT,*)'FD and HD terms:'
        write(IOUT,*)'H2',fd2,hd2
        write(IOUT,*)'H3',fd3,hd3
        write(IOUT,*)'H4',fd4,hd4
        write(IOUT,*)'H5',fd5,hd5
      end if
C
C calculate vop contributions for each H
      t1=delta23*delta23
      t2=delta24*delta24
      t3=delta25*delta25 
      voph2=(fd2*(t1+t2+t3))+
     1      ( hd2*((t1*t1)+(t2*t2)+(t3*t3)) )

      t1=delta32*delta32
      t2=delta34*delta34
      t3=delta35*delta35 
      voph3=(fd3*(t1+t2+t3))+
     1      ( hd3*((t1*t1)+(t2*t2)+(t3*t3)) )

      t1=delta42*delta42
      t2=delta43*delta43
      t3=delta45*delta45 
      voph4=(fd4*(t1+t2+t3))+
     1      ( hd4*((t1*t1)+(t2*t2)+(t3*t3)) )

      t1=delta52*delta52
      t2=delta53*delta53
      t3=delta54*delta54 
      voph5=(fd5*(t1+t2+t3))+
     1      ( hd5*((t1*t1)+(t2*t2)+(t3*t3)) )


C total vop contribution
      vop=voph2+voph3+voph4+voph5

      if(LDEBUG) then
        write(*,*) 
        write(*,*)'H2 Vop : ',voph2*1.0d2*KJ_IN_KCAL 
        write(*,*)'H3 Vop : ',voph3*1.0d2*KJ_IN_KCAL 
        write(*,*)'H4 Vop : ',voph4*1.0d2*KJ_IN_KCAL 
        write(*,*)'H5 Vop : ',voph5*1.0d2*KJ_IN_KCAL 
        write(IOUT,*) 'VOP FROM XPOT :',vop*1.0d2*KJ_IN_KCAL
      end if
C
C Calculate in-plane terms

C get the k0ij terms
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch2,rch3,rch4,rch5,xk023)
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch2,rch4,rch5,rch3,xk024)
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch2,rch5,rch3,rch4,xk025)

C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
C    1               rch3,rch2,rch5,rch4,xk032)
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch3,rch4,rch2,rch5,xk034)
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch3,rch5,rch4,rch2,xk035)

C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
C    1               rch4,rch2,rch3,rch5,xk042)
C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
C    1               rch4,rch3,rch5,rch2,xk043)
      call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,
     1               rch4,rch5,rch2,rch3,xk045)
      
C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch2,rch4,rch3,xk052)
C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch3,rch2,rch4,xk053)
C     call calcvipk0(a1s,b1s,a2s,b2s,ak,kch3,r0ch,rch5,rch4,rch3,rch2,xk054)

C
C get ki terms
      call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch2,rbh2,xk2)
      call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch3,rbh3,xk3)
      call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch4,rbh4,xk4)
      call calcvipk(aa1,aa2,aa3,aa4,r0ch,r0hh,rch5,rbh5,xk5)
C
C get bond xangles
      theta23=xangle(H2,C1,H3)
      theta24=xangle(H2,C1,H4)
      theta25=xangle(H2,C1,H5)
      theta34=xangle(H3,C1,H4)
      theta35=xangle(H3,C1,H5)
      theta45=xangle(H4,C1,H5)

      vip23=xk023*xk2*xk3*(theta23-theta023)*(theta23-theta023)
      vip24=xk024*xk2*xk4*(theta24-theta024)*(theta24-theta024)
      vip25=xk025*xk2*xk5*(theta25-theta025)*(theta25-theta025)
      vip34=xk034*xk3*xk4*(theta34-theta034)*(theta34-theta034)
      vip35=xk035*xk3*xk5*(theta35-theta035)*(theta35-theta035)
      vip45=xk045*xk4*xk5*(theta45-theta045)*(theta45-theta045)

      vip=0.50d0*(vip23+vip24+vip25+vip34+vip35+vip45)

C
C remove H5 str term
C total potential in 10(5)Joules/mol and kcal/mol
      xpot=vstr+vop+vip-v3h5-voph5-((vip25+vip35+vip45)*0.50d0)
      xpot_kcal=xpot*100.0d0*KJ_IN_KCAL

      if(LDEBUG) then
        write(*,*) 
        write(*,*)'H2-H3 Vip :',0.50d0*vip23*1.0d2*KJ_IN_KCAL 
        write(*,*)'H2-H4 Vip :',0.50d0*vip24*1.0d2*KJ_IN_KCAL 
        write(*,*)'H2-H5 Vip :',0.50d0*vip25*1.0d2*KJ_IN_KCAL 
        write(*,*)'H3-H4 Vip :',0.50d0*vip34*1.0d2*KJ_IN_KCAL 
        write(*,*)'H3-H5 Vip :',0.50d0*vip35*1.0d2*KJ_IN_KCAL 
        write(*,*)'H4-H5 Vip :',0.50d0*vip45*1.0d2*KJ_IN_KCAL 
        write(*,*) 'VIP FROM XPOT:',vip*1.0d2*KJ_IN_KCAL
        write(*,*) 'VSTR+VIP+VOP :',(vstr+vip+vop)*1.0d2*KJ_IN_KCAL

        write(*,*)
        write(*,*)'Energy after removing H5 term'
        write(*,*)'Vstr:',(vstr-v3h5)*1.0d2*KJ_IN_KCAL
        write(*,*)'Vop :',(vop-voph5)*1.0d2*KJ_IN_KCAL
        write(*,*)'Vip :',vip-((vip25+vip35+vip45)*0.50d0)*
     1                    1.0d2*KJ_IN_KCAL
        write(*,*)'Total xpot:',xpot_kcal
      end if


      if(.FALSE.) then
        R(1)=H2(1)/0.52918d0
        R(2)=H2(2)/0.52918d0
        R(3)=H2(3)/0.52918d0
        R(4)=C1(1)/0.52918d0
        R(5)=C1(2)/0.52918d0
        R(6)=C1(3)/0.52918d0
        R(7)=H3(1)/0.52918d0
        R(8)=H3(2)/0.52918d0
        R(9)=H3(3)/0.52918d0
        R(10)=H4(1)/0.52918d0
        R(11)=H4(2)/0.52918d0
        R(12)=H4(3)/0.52918d0
        R(13)=H5(1)/0.52918d0
        R(14)=H5(2)/0.52918d0
        R(15)=H5(3)/0.52918d0
        R(16)=Hb(1)/0.52918d0
        R(17)=Hb(2)/0.52918d0
        R(18)=Hb(3)/0.52918d0
        call CSURF(vsurf,R,dv,18)
C       vsurf=(vsurf-v3h5)*100.0d0*KJ_IN_KCAL
        write(*,*) 'VSURF:',vsurf*627.510d0
      write(*,*) 'MYVSURF:',(vstr+vop+vip)*0.03812D0*627.510d0
      end if
      
      return
      END 
C23456789   
C References:
C [1]. Joseph et al. : JCP, vol 87(12),  pg 7036, 1987
C [2]. Jordan et al. : JCP, vol 102(14), pg 5669, 1995
C
C************************************************
      SUBROUTINE V3qj77(alpha,d1,d3,r0,r,vq,vj)
C************************************************
C this subroutine calculates the Coulomb (vq) and 
C resonance (vj) integrals from singlet (e1) and 
C triplet (E3) states
C Equations: 2, 3, 4, 5, and 8 of Ref. [1]
C
      implicit none
C Input variables
      double precision alpha
      double precision d1
      double precision d3
      double precision r0
      double precision r
C Output variables
      double precision vq
      double precision vj
C Local variables
      double precision z1,z2
      double precision e1,e3

      z1=dexp(-alpha*(r-r0))
      z2=dexp(-2.0d0*alpha*(r-r0))

      e1=d1*(z2-(2.0d0*z1)) 
      e3=d3*(z2+(2.0d0*z1))
      vq=0.50d0*(e1+e3)
      vj=0.50d0*(e1-e3)
      return
      END 
C end of V3qj



C23456789   
C References:
C [1]. Joseph et al. : JCP, vol 87(12),  pg. 7036, 1987
C [2]. Jordan et al. : JCP, vol 102(14), pg. 5669, 1995
C [3]. Steckler et al. : JCP, vol 87(15), pg. 7024, 1987 
C
C************************************************
      FUNCTION funcS3(alpha,beta,r0,r)
C************************************************
C This func. calculates the S3 switching func
C Equation 4 of Ref. [2]
C
      implicit none
      double precision alpha
      double precision beta
      double precision r0
      double precision r
      double precision funcS3
      funcS3=1.0d0-dtanh(alpha*(r-r0)*(r-beta)*(r-beta))
      return
      END 
C end of func funcS3
C***********************************************************************
      SUBROUTINE calcfh(a3,b3,fd0,hd0,r0,ri,rj,rk,rl,fdi,hdi)
C***********************************************************************
C This sub. calculates the fdelta_i and hdelta_i
C as defined by Equation 4 of Ref. [2]
C
      implicit none
C Input variables
      double precision a3
      double precision b3
      double precision fd0
      double precision hd0
      double precision r0
      double precision ri
      double precision rj
      double precision rk
      double precision rl
C Output variables
      double precision fdi
      double precision hdi
C Local variables
      double precision s3i
      double precision s3j
      double precision s3k
      double precision s3l

      s3i=1.0d0-tanh( a3*(ri-r0)*(ri-b3)*(ri-b3) )
      s3j=1.0d0-tanh( a3*(rj-r0)*(rj-b3)*(rj-b3) )
      s3k=1.0d0-tanh( a3*(rk-r0)*(rk-b3)*(rk-b3) )
      s3l=1.0d0-tanh( a3*(rl-r0)*(rl-b3)*(rl-b3) )
      fdi=(1.0d0-s3i)*s3j*s3k*s3l*fd0
      hdi=(1.0d0-s3i)*s3j*s3k*s3l*hd0

      return
      END
C end of sub calcfh

C***********************************************************************
C23456789
      SUBROUTINE gjtheta0(PI,tau,Aphi,Bphi,Cphi,Atheta,Btheta,Ctheta,
     1                    r0,ri,rj,rk,rl,theta0ij)
C***********************************************************************
C This sub. computes the theta0 for bending terms
C Equation 5 in Ref. [2]
C
      implicit none
C Input variables
      double precision PI
      double precision tau
      double precision Aphi
      double precision Bphi
      double precision Cphi
      double precision Atheta
      double precision Btheta
      double precision Ctheta
      double precision r0
      double precision ri
      double precision rj
      double precision rk
      double precision rl
C Output variables
      double precision theta0ij
C Local variables
      double precision si
      double precision sj
      double precision sk
      double precision sl
      double precision s1
      double precision s2
      double precision t1
      double precision t2
      double precision sphii
      double precision sphij
      double precision sthetak
      double precision sthetal

      si=Aphi*(ri-r0)*dexp(Bphi*(ri-Cphi)*(ri-Cphi)*(ri-Cphi))
      sj=Aphi*(rj-r0)*dexp(Bphi*(rj-Cphi)*(rj-Cphi)*(rj-Cphi))
      sk=Atheta*(rk-r0)*dexp(Btheta*(rk-Ctheta)*(rk-Ctheta)*(rk-Ctheta))
      sl=Atheta*(rl-r0)*dexp(Btheta*(rl-Ctheta)*(rl-Ctheta)*(rl-Ctheta))
      sphii   =1.0d0-dtanh(si) 
      sphij   =1.0d0-dtanh(sj) 
      sthetak =1.0d0-dtanh(sk) 
      sthetal =1.0d0-dtanh(sl) 
      s1      =(sphii*sphij)-1.0d0 
      s2      =(sthetak*sthetal)-1.0d0 
      t1      = tau-(0.50d0*PI)
      t2      = tau-(2.0d0*PI/3.0d0)
      theta0ij=tau+(t1*s1)+(t2*s2)

      return
      END 
C end of gjtheta0
C***********************************************************************
      SUBROUTINE steckdel(v1,v2,v3,vi,deltai)
C***********************************************************************
C This sub. computes the out-of-plane delta variable
C as defined by Steckler et al. in Ref. [3] equation 18.
C Note that is *NOT* symmetric with respect to all
C the three methyl hydrogen.
C Symmetrization (as given in Ref. [2] equation 4) 
C will be done in the calling suroutine. 
C
      implicit none
C Input variables
      double precision v1(3)
      double precision v2(3)
      double precision v3(3)
      double precision vi(3)
C Output variables
      double precision deltai
C Local variables
      integer i
      double precision ri
      double precision rn
      double precision r21
      double precision r31
      double precision doti
      double precision v21(3)
      double precision v31(3)
      double precision pn(3)
C
C r21, r31 and ri are the lenghts of vectors
C v21, v31 and vi resp.
C v21=v2-v1
C v31=v3-v1
C
      do i=1,3
        v21(i)=v2(i)-v1(i)
        v31(i)=v3(i)-v1(i)
      end do
C
C the normal to the plane containing 
C points 1, 2, and 3 is given by
C pn = r21 X r31
C rn is lenght of pn
C
      call crossprod(v21,v31,pn)
C
C The xangle between pn and vi is deltaI
C evalute the dot product:
C doti = pn . vi
C
      ri  =0.0d0
      rn  =0.0d0
      doti=0.0d0
      do i=1,3
        ri  =ri+(vi(i)*vi(i))
        rn  =rn+(pn(i)*pn(i))
        doti=doti+(pn(i)*vi(i))
      end do
C doti should be -ve by convension
C for right-handed set 
      doti=-dabs(doti)
      ri  =dsqrt(ri)
      rn  =dsqrt(rn)

      deltai=dacos(doti/(ri*rn))


      return
      END
C end of sub steckdel
C23456789   
C References:
C [1]. Joseph et al. : JCP, vol 87(12),  pg. 7036, 1987
C [2]. Jordan et al. : JCP, vol 102(14), pg. 5669, 1995
C [3]. Steckler et al. : JCP, vol 87(15), pg. 7024, 1987 
C
C***********************************************************************
      SUBROUTINE calcvipk(a1,a2,a3,a4,r0ch,r0bh,rch,rbh,vipk)
C***********************************************************************
C The sub. calculates the force constant (vipk) for the in-plane 
C bendinging as given by equation 6 of Ref. [2]
C
      implicit none
C Input varaibles
      double precision a1
      double precision a2
      double precision a3
      double precision a4
      double precision r0ch
      double precision r0bh
      double precision rch
      double precision rbh
C Output variables
      double precision vipk
C Local variables
      double precision uA1
      double precision uA2

      uA1 =1.0d0-dexp(-a1*rbh*rbh)
      uA2 =a2+(a3*dexp(-a4*(rbh-r0bh)*(rbh-r0bh)))
      vipk=uA1*dexp(-uA2*(rch-r0ch)*(rch-r0ch))
      return
      END
C end of sub calcvipk
C***********************************************************************
      SUBROUTINE calcvipk0(a1,b1,a2,b2,ak,kch3,r0ch,ri,rj,rk,rl,vipk0ij)
C***********************************************************************
C This sub. calculates the equilibrium force constant (k0ij) used 
C for in-plane-bending terms
C Equation 9 of Ref. [2]
      implicit none
C Input varaibles
      double precision a1
      double precision b1
      double precision a2
      double precision b2
      double precision kch3
      double precision ak
      double precision r0ch
      double precision ri
      double precision rj
      double precision rk
      double precision rl
C Output varaibles
      double precision vipk0ij
C Local varaibles
      double precision kch4
      double precision s1i
      double precision s1j
      double precision s2k
      double precision s2l
      
      kch4=kch3+ak
C calculate switching func
      s1i=1.0d0-dtanh( a1*(ri-r0ch)*((ri-b1)**8) )
      s1j=1.0d0-dtanh( a1*(rj-r0ch)*((rj-b1)**8) )
      s2k=1.0d0-dtanh( a2*(rk-r0ch)*((rk-b2)**6) )
      s2l=1.0d0-dtanh( a2*(rl-r0ch)*((rl-b2)**6) )

C calculate k0 using Eq. 9 of Ref. [2]
      vipk0ij=kch4+( kch4*((s1i*s1j)-1.0d0) ) 
     1       +( (kch4-kch3)*((s2k*s2l)-1.0d0) )
      
      return
      END
C end of sub. calcvipk0
C23456789
C************************************************
      SUBROUTINE xlimit(a,b,s,x0,x,y)
C************************************************
C This subroutine uses a tanh function
C to smoothly connect two points (a,b)
C The tanh function is centered at x0
C with a scale factor of s.
C
C The output value is given by y
C
      implicit none
      double precision a
      double precision b
      double precision s
      double precision x0
      double precision x
      double precision y

      double precision r
      double precision ca
      double precision cb
C---------------------
C When s is +ve:
C---------------------
C As x-->b
C cb--->1 & ca-->0
C
C As x-->a
C ca--->1 & cb-->0
C
      r=x-x0
      cb=0.50d0*(1.0d0+dtanh(s*r))
      ca=1.0d0-cb
      y=(ca*a)+(cb*b)

      return
      END
C end of xlimit 
C==========================================================
C==========================================================
C==========================================================

C
C System:          CH5
C Functional form:
C Common name:     JG (2nd coding)
C Reference  :      M. J. T. Jordan, R. G. Gilbert, J. Chem. Phys., 
C                   Vol. 102, p. 5669, 1995.
C               
C Number of bodies: 6
C Number of electronic states: 1
C Number of derivatives: 0
C Interface: Section-2
c Notes: This routine contains the same JG CH5 surface as the routine 'ch5jg.f'
c        but has a different interface.

      SUBROUTINE CSURF(ENGYGS,R,DEGSDR,N3TM)
C
C        This potential is written such that:
C
C                       X(1)  - X(3)  : X, Y, Z for H1
C                       X(4)  - X(6)  : X, Y, Z for C
C                       X(7)  - X(9)  : X, Y, Z for H3
C                       X(10) - X(12) : X, Y, Z for H4
C                       X(13) - X(15) : X, Y, Z for H2
C                       X(16) - X(18) : X, Y, Z for H5
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (PI = 3.141592653589793D0)
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
      DIMENSION R(N3TM),DEGSDR(N3TM)
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)

        CALL PREPOT
C Changing units from angstrom to bohr and initialize
C
      DO 10 I = 1, 18
         q(I) = R(I)*0.52918d0
         pdot(I)=0.D0
 10   CONTINUE
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0

       call coorden
c
c  calculate switching functions
c
       call switchf
c
c  calculate reference xangles and their derivatives
c
       call refxangles
c
c  calculate stretching potential
c
       call stretch(vstr)
c
c  calculate out of plane bending potential
c
       call opbend(vop)
c
c  calculate in plane bending potential
c
       call ipbend(vip)
c
c  total potential energy is vstr+vop+vip
c
       en=vstr+vop+vip
c
c 0.03812 conversion factor from 10(5) j/mol to hartrees
c
C##--------------------------------------------------
C       write(*,*) 
C       write(*,*) 'vstr 10(5)j/mol and au : ',vstr,vstr*0.03812D0
C       write(*,*) 'vop  10(5)j/mol and au : ',vop,vop*0.03812D0
C       write(*,*) 'vip  10(5)j/mol and au : ',vip,vip*0.03812D0
C       write(*,*) 'pot  10(5)j/mol and au : ',en,en*0.03812D0
C##--------------------------------------------------
       en = en*0.03812D0
       ENGYGS = en
c
c 0.0201723 conversion factor from 10(5)j/mol/A to hartrees/bohr
c
         do i=1,18
            DEGSDR(i)=pdot(i)*0.0201723d0
         enddo
C
       return
       end
c
c******************************************************
c
       subroutine coorden
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
c  calculate relative coordinates
c
       do ind=1,3
         tcb(ind)=q(nc(ind))-q(nhb(ind))
         do i=1,4
           tch(i,ind)=q(nc(ind))-q(nh(i,ind))
           tbh(i,ind)=q(nhb(ind))-q(nh(i,ind))
         enddo
       enddo
c
c  calculate bond lengths
c
       rcb=sqrt(tcb(1)*tcb(1)+tcb(2)*tcb(2)+tcb(3)*tcb(3))
       do i=1,4
         rch(i)=sqrt(tch(i,1)*tch(i,1)+tch(i,2)*tch(i,2)+
     *                tch(i,3)*tch(i,3))
         rbh(i)=sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+
     *                tbh(i,3)*tbh(i,3))
       enddo
       return
       end
c
c******************************************************
c
c
       subroutine refxangles
c
c  subroutine calculates reference xangles for the "in-plane" potential
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       tau=acos(-1.0d0/3.0d0)
C       pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
c
c  set diagonal elements to zero
c
       do i=1,4
         theta0(i,i)=0.0d0
         do k=1,4
           dtheta0(i,i,k)=0.0d0
         enddo
       enddo
c
c  calculate reference xangles
c
C##
       theta0(1,2)=tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
       theta0(1,3)=tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
       theta0(1,4)=tau+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
       theta0(2,3)=tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
       theta0(2,4)=tau+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
       theta0(3,4)=tau+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c
c  calculate the derivatives of theta0(i,j) in terms of rch(k)
c  quantity calulated is dtheta0(i,j,k)
c
c  derivatives wrt rch(1)
c
       dtheta0(1,2,1)=(tau-halfpi)*dsphi(1)*sphi(2)
       dtheta0(1,3,1)=(tau-halfpi)*dsphi(1)*sphi(3)
       dtheta0(1,4,1)=(tau-halfpi)*dsphi(1)*sphi(4)
       dtheta0(2,3,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(4)
       dtheta0(2,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(3)
       dtheta0(3,4,1)=(tau-twopi/3.0d0)*dstheta(1)*stheta(2)
c
c  derivatives wrt rch(2)
c
       dtheta0(1,2,2)=(tau-halfpi)*sphi(1)*dsphi(2)
       dtheta0(1,3,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(4)
       dtheta0(1,4,2)=(tau-twopi/3.0d0)*dstheta(2)*stheta(3)
       dtheta0(2,3,2)=(tau-halfpi)*dsphi(2)*sphi(3)
       dtheta0(2,4,2)=(tau-halfpi)*dsphi(2)*sphi(4)
       dtheta0(3,4,2)=(tau-twopi/3.0d0)*stheta(1)*dstheta(2)
c
c  derivatives wrt rch(3)
c
       dtheta0(1,2,3)=(tau-twopi/3.0d0)*dstheta(3)*stheta(4)
       dtheta0(1,3,3)=(tau-halfpi)*sphi(1)*dsphi(3)
       dtheta0(1,4,3)=(tau-twopi/3.0d0)*stheta(2)*dstheta(3)
       dtheta0(2,3,3)=(tau-halfpi)*sphi(2)*dsphi(3)
       dtheta0(2,4,3)=(tau-twopi/3.0d0)*stheta(1)*dstheta(3)
       dtheta0(3,4,3)=(tau-halfpi)*dsphi(3)*sphi(4)
c
c  derivatives wrt rch(4)
c
       dtheta0(1,2,4)=(tau-twopi/3.0d0)*stheta(3)*dstheta(4)
       dtheta0(1,3,4)=(tau-twopi/3.0d0)*stheta(2)*dstheta(4)
       dtheta0(1,4,4)=(tau-halfpi)*sphi(1)*dsphi(4)
       dtheta0(2,3,4)=(tau-twopi/3.0d0)*stheta(1)*dstheta(4)
       dtheta0(2,4,4)=(tau-halfpi)*sphi(2)*dsphi(4)
       dtheta0(3,4,4)=(tau-halfpi)*sphi(3)*dsphi(4)
c
c  fill in the other half of the matrix
c
        do i=1,3
          do j=i+1,4
            theta0(j,i)=theta0(i,j)
            do k=1,4
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
       subroutine stretch(vstr)
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension vqch(4),vjch(4),vqbh(4),vjbh(4),vq(4),vj(4),
     *           achdc(3),achdh(4,3)
c
c  calculate avergage bond length for the methane moiety
c
       rav=(rch(1)+rch(2)+rch(3)+rch(4))/4.0d0
c
c  initialise:
c
       vstr=0.0d0
c
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
       do i=1,4
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
c  partial derivatives
c  first we need the derivative of ach:
c
       do ind=1,3
         achdc(ind)=dumach*(tch(1,ind)/rch(1)+tch(2,ind)/rch(2)
     *            +tch(3,ind)/rch(3)+tch(4,ind)/rch(4))/4.0d0
         do i=1,4
           achdh(i,ind)=-dumach*tch(i,ind)/rch(i)/4.0d0
         enddo
       enddo
       dumqcb=-acb*((d1cb+d3cb)*exp(-2.0d0*acb*(rcb-r0cb))-
     *         (d1cb-d3cb)*exp(-acb*(rcb-r0cb)))/rcb
c
c  calculate cartesian derivatives:
c  looping over ch(i) and bh(i)
c
       do i=1,4
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
           do k=1,3
             j=i+k
             if(j.gt.4)j=j-4
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
       subroutine opbend(vop)
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension sumd2(4),sumd4(4)
c
c
       vop=0.0d0
c
c  calculate force constants and their derivatives
c
       call opforce
c
c  calculate out-of-plane xangle and derivatives
c
       do i=1,4
         j=i+1
         if(j.gt.4)j=j-4
         k=j+1
         if(k.gt.4)k=k-4
         l=k+1
         if(l.gt.4)l=l-4
c
c  if i is an even number then we have switched the vectors
c  from a right handed set to a left handed set
c
c  in this case we need to switch vectors k and l around
c
         if((i.eq.2).or.(i.eq.4))then
           itemp=k
           k=l
           l=itemp
         endif
c
c  subroutine performs sum over j, k, l
c  sum2 = sum delta**2
c  sum4 = sum delta**4
c
         call calcdelta(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
       enddo
       do i=1,4
         do j=1,4
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
       subroutine ipbend(vip)
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
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
       call ipforce
c
c  calculate theta(i,j) and in plane bend potential
c
       do i=1,3
         do j=i+1,4
           costh(i,j)=tch(i,1)*tch(j,1)+tch(i,2)*tch(j,2)
     *                       +tch(i,3)*tch(j,3)
           costh(i,j)=costh(i,j)/rch(i)/rch(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
C##
C                 write(*,*) i,j,f1(j)
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
             do k=1,4
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
       return
       end
c
c*************************************************************************
c
       subroutine calcdelta(i,j,k,l,sum2,sum4)
c
c  subroutine calculates out of plane xangle delta, loops
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension  delta(4),in(3),a(3),b(3),axb(3),c(4,3),argd(4),
     *            daxb(4,3,3),cdot(4,3,3),atemp2(3)
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  set j,k,l indices
c
       in(1)=j
       in(2)=k
       in(3)=l
c
c  vector a is rk-rj, vector b is rl-rj
c
       do ind=1,3
         a(ind)=q(nh(k,ind))-q(nh(j,ind))
         b(ind)=q(nh(l,ind))-q(nh(j,ind))
       enddo

C##
C        write(*,*) 
C        write(*,*) ' i,j,k,l : ',i,j,k,l 
C        write(*,1000) ' i vec :',q(nh(i,1:3))
C        write(*,1000) ' j vec :',q(nh(j,1:3))
C        write(*,1000) ' k vec :',q(nh(k,1:3))
C        write(*,1000) ' l vec :',q(nh(l,1:3))
C        write(*,*) 
C1000    format(A,3F8.2)
c
c  axb is vector a cross b
c
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
C##
C        write(*,*) 'axb1 : ',i,axb(1)
C        write(*,*) 'axb2 : ',i,axb(2)
C        write(*,*) 'axb3 : ',i,axb(3)
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
C##
       do ii=1,3
         do ind=1,3
           c(in(ii),ind)=-tch(in(ii),ind)/rch(in(ii))
         enddo
       enddo
C##
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma

        delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
C##        write(*,*) 'theta0 ',i,in(ii),theta0(i,in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
c  derivatives of axb wrt hj:
c
       daxb(j,1,1)=0.0d0
       daxb(j,1,2)=b(3)-a(3)
       daxb(j,1,3)=-b(2)+a(2)
       daxb(j,2,1)=-b(3)+a(3)
       daxb(j,2,2)=0.0d0
       daxb(j,2,3)=b(1)-a(1)
       daxb(j,3,1)=b(2)-a(2)
       daxb(j,3,2)=-b(1)+a(1)
       daxb(j,3,3)=0.0d0
c
c  derivatives of axb wrt hk:
c
       daxb(k,1,1)=0.0d0
       daxb(k,1,2)=-b(3)
       daxb(k,1,3)=b(2)
       daxb(k,2,1)=b(3)
       daxb(k,2,2)=0.0d0
       daxb(k,2,3)=-b(1)
       daxb(k,3,1)=-b(2)
       daxb(k,3,2)=b(1)
       daxb(k,3,3)=0.0d0
c
c  derivatives of axb wrt hl:
c
       daxb(l,1,1)=0.0d0
       daxb(l,1,2)=a(3)
       daxb(l,1,3)=-a(2)
       daxb(l,2,1)=-a(3)
       daxb(l,2,2)=0.0d0
       daxb(l,2,3)=a(1)
       daxb(l,3,1)=a(2)
       daxb(l,3,2)=-a(1)
       daxb(l,3,3)=0.0d0
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
            deldot=-dtheta0(i,in(ii),i)
c
c  derivative wrt h(i,ind)
c  for  rch(i) only terms are from the derivatives of theta0
c
            deldot=-deldot*tch(i,ind)/rch(i)
            pdot(nh(i,ind))=pdot(nh(i,ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
c
c  derivative wrt c(ind)
c
            deldot=-deldot
            pdot(nc(ind))=pdot(nc(ind))
     *                   +2.0d0*fdelta(i)*delta(in(ii))*deldot
     *                   +4.0d0*hdelta(i)*delta(in(ii))**3*deldot
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
            atemp5=-dtheta0(i,in(ii),in(jj))
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
       subroutine opforce
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension switch(4),dswitch(4,4)
c
c  calculate switching functions:
c
       switch(1)=(1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
       switch(2)=(1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
       switch(3)=(1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
       switch(4)=(1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
c
c  calculate derivatives:
c  derivative of switch(1) wrt the 4 rch bond lengths
c
       dswitch(1,1)=-ds3(1)*s3(2)*s3(3)*s3(4)
       dswitch(1,2)=(1.0d0-s3(1))*ds3(2)*s3(3)*s3(4)
       dswitch(1,3)=(1.0d0-s3(1))*s3(2)*ds3(3)*s3(4)
       dswitch(1,4)=(1.0d0-s3(1))*s3(2)*s3(3)*ds3(4)
c
c  derivative of switch(2) wrt the 4 rch bond lengths
c
       dswitch(2,1)=(1.0d0-s3(2))*s3(3)*s3(4)*ds3(1)
       dswitch(2,2)=-ds3(2)*s3(3)*s3(4)*s3(1)
       dswitch(2,3)=(1.0d0-s3(2))*ds3(3)*s3(4)*s3(1)
       dswitch(2,4)=(1.0d0-s3(2))*s3(3)*ds3(4)*s3(1)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(3,1)=(1.0d0-s3(3))*s3(4)*ds3(1)*s3(2)
       dswitch(3,2)=(1.0d0-s3(3))*s3(4)*s3(1)*ds3(2)
       dswitch(3,3)=-ds3(3)*s3(4)*s3(1)*s3(2)
       dswitch(3,4)=(1.0d0-s3(3))*ds3(4)*s3(1)*s3(2)
c
c  derivative of switch(3) wrt the 4 rch bond lengths
c
       dswitch(4,1)=(1.0d0-s3(4))*ds3(1)*s3(2)*s3(3)
       dswitch(4,2)=(1.0d0-s3(4))*s3(1)*ds3(2)*s3(3)
       dswitch(4,3)=(1.0d0-s3(4))*s3(1)*s3(2)*ds3(3)
       dswitch(4,4)=-ds3(4)*s3(1)*s3(2)*s3(3)
c
c  calculate the force constants and their derivatives
c
       do i=1,4
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
         do j=1,4
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
       subroutine ipforce
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
       dimension dfk0(4,4,4),df1dc(4),df1dh(4)
c
c  set force constant at asymptotes
c
       f0=fkinf+ak
       f2=fkinf
c
C##
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
       fk0(1,4)=f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
       fk0(2,4)=f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
       fk0(3,4)=f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)

C##
C           write(*,*)
C           write(*,*) 'f0    : ',f0
C           write(*,*) 'f2    : ',f2
C           write(*,*) 's1(1) : ',s1(1)
C           write(*,*) 's1(2) : ',s1(2)
C           write(*,*) 's2(3) : ',s2(3)
C           write(*,*) 's2(4) : ',s2(4)
C           write(*,*)

c
c  derivative of fk0
c
       dfk0(1,2,1)=f0*ds1(1)*s1(2)
       dfk0(1,2,2)=f0*s1(1)*ds1(2)
       dfk0(1,2,3)=(f0-f2)*ds2(3)*s2(4)
       dfk0(1,2,4)=(f0-f2)*s2(3)*ds2(4)
c
       dfk0(1,3,1)=f0*ds1(1)*s1(3)
       dfk0(1,3,2)=(f0-f2)*ds2(2)*s2(4)
       dfk0(1,3,3)=f0*s1(1)*ds1(3)
       dfk0(1,3,4)=(f0-f2)*s2(2)*ds2(4)
c
       dfk0(1,4,1)=f0*ds1(1)*s1(4)
       dfk0(1,4,2)=(f0-f2)*ds2(2)*s2(3)
       dfk0(1,4,3)=(f0-f2)*s2(2)*ds2(3)
       dfk0(1,4,4)=f0*s1(1)*ds1(4)
c
       dfk0(2,3,1)=(f0-f2)*ds2(1)*s2(4)
       dfk0(2,3,2)=f0*ds1(2)*s1(3)
       dfk0(2,3,3)=f0*s1(2)*ds1(3)
       dfk0(2,3,4)=(f0-f2)*s2(1)*ds2(4)
c
       dfk0(2,4,1)=(f0-f2)*ds2(1)*s2(3)
       dfk0(2,4,2)=f0*ds1(2)*s1(4)
       dfk0(2,4,3)=(f0-f2)*s2(1)*ds2(3)
       dfk0(2,4,4)=f0*s1(2)*ds1(4)
c
       dfk0(3,4,1)=(f0-f2)*ds2(1)*s2(2)
       dfk0(3,4,2)=(f0-f2)*s2(1)*ds2(2)
       dfk0(3,4,3)=f0*ds1(3)*s1(4)
       dfk0(3,4,4)=f0*s1(3)*ds1(4)
c
       do i=1,4
c
c  calculate the terms f1(i)
c
C##
         arga1=aa1*rbh(i)*rbh(i)
         arga2=aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rch(i)-r0ch)**2)

C           write(*,*) 'A1 :',a1
C           write(*,*) 'A2 :',a2
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
      dfdc(1,2,1)=dfk0(1,2,1)*f1(1)*f1(2)+fk0(1,2)*df1dc(1)*f1(2)
      dfdc(1,2,2)=dfk0(1,2,2)*f1(1)*f1(2)+fk0(1,2)*f1(1)*df1dc(2)
      dfdc(1,2,3)=dfk0(1,2,3)*f1(1)*f1(2)
      dfdc(1,2,4)=dfk0(1,2,4)*f1(1)*f1(2)
c
      dfdc(1,3,1)=dfk0(1,3,1)*f1(1)*f1(3)+fk0(1,3)*df1dc(1)*f1(3)
      dfdc(1,3,2)=dfk0(1,3,2)*f1(1)*f1(3)
      dfdc(1,3,3)=dfk0(1,3,3)*f1(1)*f1(3)+fk0(1,3)*f1(1)*df1dc(3)
      dfdc(1,3,4)=dfk0(1,3,4)*f1(1)*f1(3)
c
      dfdc(1,4,1)=dfk0(1,4,1)*f1(1)*f1(4)+fk0(1,4)*df1dc(1)*f1(4)
      dfdc(1,4,2)=dfk0(1,4,2)*f1(1)*f1(4)
      dfdc(1,4,3)=dfk0(1,4,3)*f1(1)*f1(4)
      dfdc(1,4,4)=dfk0(1,4,4)*f1(1)*f1(4)+fk0(1,4)*f1(1)*df1dc(4)
c
      dfdc(2,3,1)=dfk0(2,3,1)*f1(2)*f1(3)
      dfdc(2,3,2)=dfk0(2,3,2)*f1(2)*f1(3)+fk0(2,3)*df1dc(2)*f1(3)
      dfdc(2,3,3)=dfk0(2,3,3)*f1(2)*f1(3)+fk0(2,3)*f1(2)*df1dc(3)
      dfdc(2,3,4)=dfk0(2,3,4)*f1(2)*f1(3)
c
      dfdc(2,4,1)=dfk0(2,4,1)*f1(2)*f1(4)
      dfdc(2,4,2)=dfk0(2,4,2)*f1(2)*f1(4)+fk0(2,4)*df1dc(2)*f1(4)
      dfdc(2,4,3)=dfk0(2,4,3)*f1(2)*f1(4)
      dfdc(2,4,4)=dfk0(2,4,4)*f1(2)*f1(4)+fk0(2,4)*f1(2)*df1dc(4)
c
      dfdc(3,4,1)=dfk0(3,4,1)*f1(3)*f1(4)
      dfdc(3,4,2)=dfk0(3,4,2)*f1(3)*f1(4)
      dfdc(3,4,3)=dfk0(3,4,3)*f1(3)*f1(4)+fk0(3,4)*df1dc(3)*f1(4)
      dfdc(3,4,4)=dfk0(3,4,4)*f1(3)*f1(4)+fk0(3,4)*f1(3)*df1dc(4)
c
c  derivative of total force constant f(i,j) wrt bond length rbh(k)
c  is given by dfdh(i,j,k)
c
c  only non-zero derivatives are those from rbh(i) and rbh(j)
c
       dfdh(1,2,1)=fk0(1,2)*df1dh(1)*f1(2)
       dfdh(1,2,2)=fk0(1,2)*f1(1)*df1dh(2)
       dfdh(1,2,3)=0.0d0
       dfdh(1,2,4)=0.0d0
c
       dfdh(1,3,1)=fk0(1,3)*df1dh(1)*f1(3)
       dfdh(1,3,2)=0.0d0
       dfdh(1,3,3)=fk0(1,3)*f1(1)*df1dh(3)
       dfdh(1,3,4)=0.0d0
c
       dfdh(1,4,1)=fk0(1,4)*df1dh(1)*f1(4)
       dfdh(1,4,2)=0.0d0
       dfdh(1,4,3)=0.0d0
       dfdh(1,4,4)=fk0(1,4)*f1(1)*df1dh(4)
c
       dfdh(2,3,1)=0.0d0
       dfdh(2,3,2)=fk0(2,3)*df1dh(2)*f1(3)
       dfdh(2,3,3)=fk0(2,3)*f1(2)*df1dh(3)
       dfdh(2,3,4)=0.0d0
c
       dfdh(2,4,1)=0.0d0
       dfdh(2,4,2)=fk0(2,4)*df1dh(2)*f1(4)
       dfdh(2,4,3)=0.0d0
       dfdh(2,4,4)=fk0(2,4)*f1(2)*df1dh(4)
c
       dfdh(3,4,1)=0.0d0
       dfdh(3,4,2)=0.0d0
       dfdh(3,4,3)=fk0(3,4)*df1dh(3)*f1(4)
       dfdh(3,4,4)=fk0(3,4)*f1(3)*df1dh(4)
c
       return
       end
c
c******************************************************
c
c
       subroutine switchf
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
      COMMON/PT1CM/ R(N3ATOM), ENGYGS, DEGSDR(N3ATOM)
      COMMON/PT3CM/ EZERO(ISURF+1)
      COMMON/PT4CM/ ENGYES(ISURF), DEESDR(N3ATOM,ISURF)
      COMMON/PT5CM/ ENGYIJ(JSURF), DEIJDR(N3ATOM,JSURF)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
      COMMON/USROCM/ PENGYGS,PENGYES(ISURF),
     +               PENGYIJ(JSURF),
     +               DGSCART(NATOM,3),DESCART(NATOM,3,ISURF),
     +               DIJCART(NATOM,3,JSURF)
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /xangles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
       common /op1/     s3(4),ds3(4)
       common /qpdot_pl/   q(150),pdot(150)
       common /switch1/ sphi(4),dsphi(4),stheta(4),dstheta(4)
C
C##
       a1s=1.5313681d-7
       b1s=-4.6696246d0
       a2s=1.0147402d-7
       b2s=-12.363798d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
       do i=1,4
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
C##          write(*,*) i,s1(i),rch(i)

      
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
      SUBROUTINE PREPOT
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
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
       common /ndx/     nc(3),nhb(3),nh(4,3)
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
       REF(1)='M. J. T. Jordan and R. G. Gilbert'
       REF(2)='J. Chem. Phys. Vol. 102, p. 5669, 1995'
C
      INDEXES(1) = 1
      INDEXES(2) = 6
      INDEXES(3) = 1
      INDEXES(4) = 1
      INDEXES(5) = 1
      INDEXES(6) = 1
C
C
C
      IRCTNT=6
C
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc(ind)=3*nnc+icount
         nhb(ind)=3*nnb+icount
         do i=1,4
           nh(i,ind)=3*nnh(i)+icount
         enddo
       enddo
c
c  convert to appropriate units:
c
c  energy   in 1.0d+05 j/mol
c  time     in 1.0d-14 s
c  xdistance in 1.0d-10 m
c  xangles   in radians
c  mass     in amu
c
       fact1=0.041840d0
       fact2=6.022045d0
c
C##----------------------------------
       d1ch=d1ch*fact1
       d3ch=d3ch*fact1
       d1cb=d1cb*fact1
       d3cb=d3cb*fact1
       d1hh=d1hh*fact1
       d3hh=d3hh*fact1
       fch3=fch3*fact2
       hch3=hch3*fact2
       fkinf=fkinf*fact2
       ak=ak*fact2
C##----------------------------------
C
      RETURN
      END
        SUBROUTINE CSETUP(N3TM)
C       dummy subroutine to make the code compatable with chk.f

        END 
C
      BLOCK DATA PTPACM
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
      COMMON/PT3CM/ EZERO(ISURF+1)
C
      COMMON/INFOCM/ CARTNU(NATOM,3),INDEXES(NATOM),
     +               IRCTNT,NATOMS,ICARTR,MDER,MSURF,REF
C
C
      COMMON/USRICM/ CART(NATOM,3),ANUZERO,
     +               NULBL(NATOM),NFLAG(20),
     +               NASURF(ISURF+1,ISURF+1),NDER
C
      COMMON /POTCM/ nnc,nnb,nnh(4),
     +               r0ch,d1ch,d3ch,
     +               a1ch,b1ch,c1ch,
     +               r0hh,d1hh,d3hh,ahh,
     +               r0cb,d1cb,d3cb,acb,
     +               a3s,b3s,aphi,bphi,cphi,
     +               atheta,btheta,ctheta,
     +               fch3,hch3,
     +               fkinf,ak,bk,aa1,aa2,aa3,aa4
C
       common /angles/  theta0(4,4),dtheta0(4,4,4)
       common /bonds/   rcb,rch(4),rbh(4)
       common /coords/  tcb(3),tch(4,3),tbh(4,3)
       common /delta1/  fdelta(4),hdelta(4)
       common /delta2/  dfdelta(4,4),dhdelta(4,4)
       common /force1/  fk0(4,4),f1(4),dfdc(4,4,4),dfdh(4,4,4)
       common /fsw1/    a1s,b1s,a2s,b2s
       common /ip1/     s1(4),ds1(4),s2(4),ds2(4)
       common /ndx/     nc(3),nhb(3),nh(4,3)
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
       DATA nnc     /2/
       DATA nnb     /6/
       DATA nnh     /3,4,5,1/
       DATA r0ch    /   1.09397d0/
       DATA d1ch    / 112.230d0/
       DATA d3ch    /  38.834d0/
       DATA a1ch    /   1.71300d0/
       DATA b1ch    /   0.13500d0/
       DATA c1ch    /   6.61404d0/
       DATA r0hh    /   0.74191d0/
       DATA d1hh    / 109.458d0/
       DATA d3hh    /  39.664d0/
       DATA ahh     /   1.9457d0/
       DATA r0cb    /   1.09397d0/
       DATA d1cb    /  26.409d0/
       DATA d3cb    /  20.063d0/
       DATA acb     /   1.8530285d0/
       DATA a3s     /   0.1419147d0/
       DATA b3s     /  -0.3068450d0/
       DATA aphi    /   0.5287903d0/
       DATA bphi    /   0.4006638d0/
       DATA cphi    /   1.9209937d0/
       DATA atheta  /   0.9078714d0/
       DATA btheta  /   0.3548859d0/
       DATA ctheta  /   1.8915497d0/
       DATA fch3    /   0.0957500d0/
       DATA hch3    /   0.1915000d0/
       DATA fkinf   /   0.4077000d0/
       DATA ak      /   0.1260000d0/
       DATA bk      /  10.7132d0/
       DATA aa1     /   3.213952d0/
       DATA aa2     /   1.599963d0/
       DATA aa3     /   2.165953d0/
       DATA aa4     /  11.569977d0/
C
       END

C23456789
      SUBROUTINE xgaussfit(r,fit)
      implicit none
      double precision r
      double precision fit

      integer NG
      parameter (NG=6)

      integer i,j,k
      double precision d
      double precision ans
      double precision r0(NG)
      double precision alpha(NG)
      double precision cf(NG)
      double precision fx
    

      r0(1) = 1.949230d0
      r0(2) = 1.013860d0
      r0(3) = 2.790270d0
      r0(4) = (r0(1)+r0(2))*0.50d0
      r0(5) = (r0(1)+r0(3))*0.50d0
      r0(6) = (r0(2)+r0(3))*0.50d0
      do i=1,NG
        alpha(i)=0.250d0
      end do

      cf(1)=-1178.890067483491520d0
      cf(2)=32.80653355230102620d0
      cf(3)=10.31188564665787060d0
      cf(4)=-159.8465369524321600d0
      cf(5)=17.33838608018073660d0
      cf(6)=1280.367802966489080d0

      ans=0.0d0
      do i=1,NG
        d=r-r0(i) 
        fx=dexp(-alpha(i)*d*d)
        ans=ans+(cf(i)*fx)
      end do

      fit=ans

      return
      END 

