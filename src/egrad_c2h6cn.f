c:*************************************************************************
      SUBROUTINE egrad_c2h6cn(V, COORD, DX, N3TM)
C
C   System:    C2H6 + CN --> C2H5 + HCN.
C              Joaquin Espinosa-Garcia and Cipriano Rangel, Junio 2023
C
C
C   All the information passed to and from the potential energy surface
C   routine is in hartree atomic units.
C
C        This potential is written such that:
C                       X(1)  - X(3)  : X, Y, Z for C1
C                       X(4)  - X(6)  : X, Y, Z for C2
C                       X(7)  - X(9)  : X, Y, Z for H1
C                       X(10) - X(12) : X, Y, Z for H2
C                       X(13) - X(15) : X, Y, Z for H3
C                       X(16) - X(18) : X, Y, Z for H4
C                       X(19) - X(21) : X, Y, Z for H5
C                       X(22) - X(24) : X, Y, Z for H6
C                       X(25) - X(27) : X, Y, Z for C
C                       X(28) - X(30) : X, Y, Z for N
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C CRD 2019 aumentamos de 27 a 30 las matrices
C
      DIMENSION COORD(30),DX(30)
      DIMENSION dxvstr(30),dxvtor(30),dxvop(30),dxvip(30)
C
      include 'c2h6cn.inc'
C
C     PUT COORDINATES IN PROPER ARRAYS
C
C Inicialization de variables
C
c CRD 2019 inicializamos a 30
C
      DO I=1,30
          q(I)=0.D0
          pdot(I)=0.D0
      ENDDO
c
C Changing units from angstrom to bohr
C
c CRD 2019 aumentada la matriz a 30
c
      DO 10 I = 1, 30
         q(I) = COORD(I)*0.52918d0
 10   CONTINUE
c
c  Read the parameters from the input file CONST
c
c  FORMATTED read statements
c
c       open(unit=7,file='CONST',status='old')
c
c  read in indices for the carbon and hydrogen atoms
c
c       read(7,80)comlin
c       read(7,80)comlin
c       read(7,80)comlin
c
c      read(7,910)nnc1,nnc2,nnb,(nnh(i),i=1,6)
       nnc1=1
       nnc2=2
       nnb=9
       nnh(1)=3
       nnh(2)=4
       nnh(3)=5
       nnh(4)=6
       nnh(5)=7
       nnh(6)=8
       nno=10
C CRD 2019 añadimos el H del OH como nno
C
c       read(7,910)nnc1,nnc2,nnb,(nnh(i),i=1,6),nno
c
c80     format(a79)
c910    format(9i4)
c910    format(10i4)
c
c  calculate indexes for coordinates
c
       do ind=1,3
         icount=ind-3
         nc1(ind)=3*nnc1+icount
         nc2(ind)=3*nnc2+icount
         nhb(ind)=3*nnb+icount
c
c CRD 2019 añadimos el caluculo de las coordenadas de H
c
         no(ind)=3*nno+icount
c
c fin CRD 2019
c
         do i=1,6
           nh(i,ind)=3*nnh(i)+icount
         enddo
       enddo
c
c  read in parameters for the stretching term
c
c morse parameters for ch bonds: r0c1hr,r0c1hp,w3,w4, d1c1h, d3c1h
       r0c1hr=1.08897
       r0c1hp=1.07577
       w3=1.00000
       w4=1.08897
       d1c1h=109.23000
       d3c1h=23.00400
c attenuation parameters for ach: a1c1h, b1c1h, c1c1h
       a1c1h=1.68000
       b1c1h=0.14500
       c1c1h=8.21404
c morse parameters for ch bonds: r0c2hr,r0c2hp,w5,w6, d1c2h, d3c2h
       r0c2hr=1.08897
       r0c2hp=1.07577
       w5=1.00000
       w6=1.08897
       d1c2h=109.23000
       d3c2h=23.00400
c attenuation parameters for ach: a1c2h, b1c2h, c1c2h
       a1c2h=1.68000
       b1c2h=0.14500
       c1c2h=8.21404
c morse parameters for cx bonds: r0cc, d1cc, acc
       r0cc=1.52007
       d1cc=117.949
       acc=1.934
c parametros para el fk0: a1c, a2c
       a1c=0.90007
       a2c=0.949
c morse parameters for hh bonds: r0hh, d1hh, d3hh and ahh
       r0hh=1.07101
       d1hh=131.858
       d3hh=35.564
       ahh=1.9057
c morse parameters for cb bonds: r0cb, d1cb 112, d3cb 38 and acb  1.77
       r0cb=1.73397
       d1cb=56.009
       d3cb=18.470
       acb=1.1910285
c out-of-plane bend: a3s, b3s, aphi, bphi, cphi, rphi
       a3s=1.8419147
       b3s=1.088930
       aphi=1.0087903
       bphi=1.2006638
       cphi=1.0889937
       rphi=1.08897
c out-of-plane bend: atheta, btheta, ctheta
       atheta=1.2078714
       btheta=0.3548859
       ctheta=1.8915497
c out-of-plane bend: fch, hch3
       fch=0.1357500
       hch3=0.2615000
c in-plane bend: fkinf, ak, bk, aa1, aa2, aa3, aa4, aa5 (no se usa)
       fkinf=0.2977000
       ak=0.306000
       bk=0.2032000
       aa1=0.450952
       aa2=0.499963
       aa3=6.465953
       aa4=1.569977
       aa5=0.1526290
c h2o in-plane bend: fkh2oeq, alph2o, angh2oeq
       fkh2oeq=0.26000
       alph2o=3.20800
       angh2oeq=180.00000
c vtorsion V3, w1, w2
       V3=13.7300000
       w1=9.3548859
       w2=0.91397
c a1s b1s a2s b2s 1.5313681 -4.669625 1.0147402 -12.36380
       a1s=2.5313681
       b1s=1.088625
       a2s=2.0147402
       b2s=1.088980
c a1scc b1scc a2scc b2scc a3scc b3scc
       a1scc=2.5313681
       b1scc=1.550025
       a2scc=2.0147402
       b2scc=1.550080
       a3scc=0.531368
       b3scc=1.550080
c vdwe vdwr vdwc1 vdwc2, vdt1, vdt2, vdtr anulado en el fuente
       vdwe=0.0049500
       vdwr=2.183100
       vdwc1=-6.700000
       vdwc2=6.500000
       vdt1=1.000000
       vdt2=8.000000
       vdtr=1.088390


c       read(7,80)comlin
c      read(7,930)r0ch,d1ch,d3ch
c       read(7,930)r0c1hr,r0c1hp,w3,w4,d1c1h,d3c1h
c       read(7,80)comlin
c       read(7,930)a1c1h,b1c1h,c1c1h
c       read(7,80)comlin
c       read(7,930)r0c2hr,r0c2hp,w5,w6,d1c2h,d3c2h
c       read(7,80)comlin
c       read(7,930)a1c2h,b1c2h,c1c2h
c
c  read in parameters for the stretching term C-C
c
c       read(7,80)comlin
c       read(7,930)r0cc,d1cc,acc
c       read(7,80)comlin
c       read(7,930)a1c,a2c
c
c       read(7,80)comlin
c       read(7,930)r0hh,d1hh,d3hh,ahh
c
c       read(7,80)comlin
c       read(7,930)r0cb,d1cb,d3cb,acb
c
c930    format(8f10.5)
c
c  read in parameters for the out of plane bending term
c
c       read(7,80)comlin
c       read(7,930)a3s,b3s,aphi,bphi,cphi,rphi

c       read(7,80)comlin
c       read(7,930)atheta,btheta,ctheta

c       read(7,80)comlin
c       read(7,930)fch3,hch3
c
c  read in parameters for the in plane bending term
c
c       read(7,80)comlin
c       read(7,930)fkinf,ak,bk,aa1,aa2,aa3,aa4,aa5
c
c CRD 2019 añadimos los terminos del H2O
c
c
c  read in parameters for the h2o bending term
c
c       read(7,80) comlin
c       read(7,930) fkh2oeq, alph2o, angh2oeq
c fin CRD 2019
c
c CRD 2011 terminos del CONST de vtorsion
c
c       read(7,80)comlin
c       read(7,930)V3,w1,w2
c jcc
c       read(7,80)comlinc
c       read(7,930)a1s,b1s,a2s,b2s
c       read(7,80)comlin
c       read(7,930)a1scc,b1scc,a2scc,b2scc,a3scc,b3scc
c
c      read(7,80)comlin
c      read(7,930)vdwe, vdwr, vdwc1, vdwc2, vdt1, vdt2, vdtr
c
c jcc
c
c
c  convert to appropriate units:
c
c  From:
c
c  kcal/mol
c  mdyn A-1        -> 1.0d+05 j/mol...
c  mdyn A rad-1
c
c  To
c
c  energy   in 1.0d+05 j/mol
c
c distance in 1.0d-10 m
c  angles   in radians
c  mass     in amu
c
       fact1=0.041840d0
       fact2=6.022045d0
c CRD 2019 fact3
       fact3=2.d0*3.1415926d0/360.d0
c fin CRD 2019
       d1ch=d1ch*fact1
       d3ch=d3ch*fact1
       d1c1h=d1c1h*fact1
       d3c1h=d3c1h*fact1
       d1c2h=d1c2h*fact1
       d3c2h=d3c2h*fact1
       d1cc=d1cc*fact1
c      d3cc=d3cc*fact1
       d1bc=d1bc*fact1
       d3bc=d3bc*fact1
       d1cb=d1cb*fact1
       d3cb=d3cb*fact1
       d1hh=d1hh*fact1
       d3hh=d3hh*fact1
       fch3=fch3*fact2
       hch3=hch3*fact2
       fkinf=fkinf*fact2
       ak=ak*fact2
c CRD 2019
       fkh2oeq=fkh2oeq*fact2
       angh2oeq=angh2oeq*fact3
c fin CRD 2019
c
       a1s=a1s*1.D-8
       a2s=a2s*1.D-8
       a1scc=a1scc*1.D-8
       a2scc=a2scc*1.D-8
       vdwc2=vdwc2*1.D+4
c
       close(unit=7)
c
c
c  calculate relative coordinates and bond lengths
c
       en=0.0d0
       call coorden_c2h6cn
c
c CRD 2011 incluimos los terminos de potencial de torsion
c
       call torsion_c2h6cn(vtor)
c
c  calculate switching functions
c
       call switchf_c2h6cn
c
c  calculate reference angles and their derivatives
c
       call refangles_c2h6cn
c
c  calculate stretching potential
c
       call stretch_c2h6cn(vstr)
c
c  calculate out of plane bending potential
c
       call opbend_c2h6cn(vop)
c
c  calculate in plane bending potential
c
       call ipbend_c2h6cn(vip)
c
c
c jcc calculate vdw term
c
c      call vander(evdw)
c
C calculo de la energia
C
       en=vstr+vop+vip+vtor
C
c  convert from 10(5) j/mol to au
c
       en = en*0.03812D0
       V = en
c
c  Initialize derivatives
c
       do n=1,30
           DX(n)=0.0d0
       enddo
c
c  copy the pdots to the
c  appropriate elements in dx and
c  convert from 10(5) j/mol/A to au/bohr
c
CC       do ind=1,30
CC           DX(ind)=pdot(ind)*0.0201723d0
CC       enddo
c
c JEG Calculo de derivadas numericas de la energia con
c     respecto a las coordenadas
c
      PASO= 1.0D-7
      do I=1,30
        q(I)=q(I) + PASO
        call coorden_c2h6cn
        call switchf_c2h6cn
        call refangles_c2h6cn
        call torsion_c2h6cn(vtor)
        call stretch_c2h6cn(vstr)
        call opbend_c2h6cn(vop)
        call ipbend_c2h6cn(vip)
        en=vstr+vop+vip+vtor
C
c  convert from 10(5) j/mol to au
c
        en=en*0.03812D0
c
        DX(I)=(en-V)/PASO
        DX(I)=DX(I)*0.529177d0
        q(I)=q(I)-PASO
      enddo
c
       return
       end
c******************************************************
c
c------------------------------------------------------
       subroutine coorden_c2h6cn
c------------------------------------------------------
c
c  calculates relative coordinates and bond lengths
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
c  calculate relative coordinates
c
       do ind=1,3
         tc1b(ind)=q(nc1(ind))-q(nhb(ind))
         tc2b(ind)=q(nc2(ind))-q(nhb(ind))
         tcc(ind)=q(nc1(ind))-q(nc2(ind))
c CRD 2019
         tno(ind)=q(no(ind))-q(nhb(ind))
c
       enddo
c
       do ind=1,3
         do i=1,3
           tc1h(i,ind)=q(nc1(ind))-q(nh(i,ind))
           tc2h(i,ind)=q(nc2(ind))-q(nh(i+3,ind))
           tbh(i,ind)=q(nhb(ind))-q(nh(i,ind))
           tbh(i+3,ind)=q(nhb(ind))-q(nh(i+3,ind))
         enddo
       enddo
c
c
c coordenadas entre los hidrogenos de los dos carbonos
c
c      do i=1,3
c         do j=4,6
c            do ind=1,3
c          thh(i,j,ind)=q(nh(j,ind))-q(nh(i,ind))
c            enddo
c         enddo
c      enddo
c
c  calculate bond lengths
c
       rc1b=sqrt(tc1b(1)*tc1b(1)+tc1b(2)*tc1b(2)+tc1b(3)*tc1b(3))
       rc2b=sqrt(tc2b(1)*tc2b(1)+tc2b(2)*tc2b(2)+tc2b(3)*tc2b(3))
       rcc=sqrt(tcc(1)*tcc(1)+tcc(2)*tcc(2)+tcc(3)*tcc(3))
c CRD 2019
       rno=sqrt(tno(1)*tno(1)+tno(2)*tno(2)+tno(3)*tno(3))
c fin CRD 2019
       do i=1,3
         rc1h(i)=sqrt(tc1h(i,1)*tc1h(i,1)+tc1h(i,2)*tc1h(i,2)+
     *                tc1h(i,3)*tc1h(i,3))
c
         rc2h(i)=sqrt(tc2h(i,1)*tc2h(i,1)+tc2h(i,2)*tc2h(i,2)+
     *                tc2h(i,3)*tc2h(i,3))
c
       enddo
       do i=1,6
         rbh(i)=sqrt(tbh(i,1)*tbh(i,1)+tbh(i,2)*tbh(i,2)+
     *                tbh(i,3)*tbh(i,3))
       enddo
c calculamos las distancias H-H
c
c      do i=1,3
c         do j=4,6
c        rhh(i,j)=sqrt(thh(i,j,1)*thh(i,j,1)+thh(i,j,2)*thh(i,j,2)+
c    *                thh(i,j,3)*thh(i,j,3))
c         enddo
c      enddo
c
c crd 2013 modificaciones de las distancias
       argmax=19.d0
       P3=1.d0
       do i=1,3
          argp3=(w3*(rc1h(i)-w4))
          if (argp3.lt.argmax) then
            t3tmp=1-tanh(argp3)
          else
            t3tmp=0.d0
          endif
          P3=P3*t3tmp
       enddo
       r0c1h=P3*r0c1hr+(1-P3)*r0c1hp
c      r0c1h=r0c1hr
c
       P4=1.d0
       do i=1,3
          argp4=(w5*(rc2h(i)-w6))
          if (argp4.lt.argmax) then
            t4tmp=1-tanh(argp4)
          else
            t4tmp=0.d0
          endif
          P4=P4*t4tmp
       enddo
       r0c2h=P4*r0c2hr+(1-P4)*r0c2hp
c      r0c2h=r0c2hr
c      write(35,*) 'r0c1h', r0c1h, 'r0c2h',r0c2h
c
       return
       end
c
c******************************************************
c
c------------------------------------------------------
c CRD 2011 subrutina para calcular el potencial de torsion
       subroutine torsion_c2h6cn(vtor)
c------------------------------------------------------
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
       dimension t1(3), t2(3)
       double precision pi
       pi=4.0d0*atan(1.0d0)
c
c CRD 2011 incluimos un apartado para calcular los angulos dihedro
c para poder calcular despues vtorsion
c
       call dihedro_c2h6cn(3,1,2,6,dihed)
       thetahcch(1)=dihed*pi/180.0d0
       call dihedro_c2h6cn(3,1,2,7,dihed)
       thetahcch(2)=dihed*pi/180.0d0
       call dihedro_c2h6cn(3,1,2,8,dihed)
       thetahcch(3)=dihed*pi/180.0d0
       call dihedro_c2h6cn(4,1,2,6,dihed)
       thetahcch(4)=dihed*pi/180.0d0
       call dihedro_c2h6cn(4,1,2,7,dihed)
       thetahcch(5)=dihed*pi/180.0d0
       call dihedro_c2h6cn(4,1,2,8,dihed)
       thetahcch(6)=dihed*pi/180.0d0
       call dihedro_c2h6cn(5,1,2,6,dihed)
       thetahcch(7)=dihed*pi/180.0d0
       call dihedro_c2h6cn(5,1,2,7,dihed)
       thetahcch(8)=dihed*pi/180.0d0
       call dihedro_c2h6cn(5,1,2,8,dihed)
       thetahcch(9)=dihed*pi/180.0d0
c
c CRD 2011 calculamos los terminos switching t1 y t2
c
      do i=1,3
        t1(i)=0.5d0*(1-tanh(w1*(rc1h(i)-w2)))
        t2(i)=0.5d0*(1-tanh(w1*(rc2h(i)-w2)))
c       write(*,*) 't1',i,t1(i),'t2',i,t2(i)
c       t1(i)=1-tanh(w1*(rc1h(i)-w2))
c       t2(i)=1-tanh(w1*(rc2h(i)-w2))
      enddo
c
c CRD 2011 inicializamos vtor
c
       vtor=0.0d0
c
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(1)))*t1(1)*t2(1)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(2)))*t1(1)*t2(2)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(3)))*t1(1)*t2(3)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(4)))*t1(2)*t2(1)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(5)))*t1(2)*t2(2)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(6)))*t1(2)*t2(3)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(7)))*t1(3)*t2(1)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(8)))*t1(3)*t2(2)
       vtor=vtor+(V3/3.0d0)*(1+cos(3*thetahcch(9)))*t1(3)*t2(3)
       return
       end
c
c------------------------------------------------------
c      double precision function dihed(i,j,k,l)
       subroutine dihedro_c2h6cn(i,j,k,l, dihed)
c------------------------------------------------------
c
         implicit double precision (a-h,o-z)
         include 'c2h6cn.inc'
         integer i, j, k,l
         integer m
c        double precision x(300)
c        double precision dista
c        double precision d12,d13,d14,d23,d24,d34
c        double precision p,q
c        double precision p
         double precision pi
c
         double precision anorm_c2h6cn
         double precision dotprod_c2h6cn
         double precision rij(3),rjk(3),rkj(3),rkl(3)
         double precision zijjk(3),zijkj(3), zkjkl(3), zjkkl(3)
         double precision tmpsign(3), tmp, arg1,arg2
c
         pi=4.0d0*atan(1.0d0)
c
         do m = 1, 3
            rij(m)=q((i-1)*3+m)-q((j-1)*3+m)
            rjk(m)=q((j-1)*3+m)-q((k-1)*3+m)
            rkl(m)=q((k-1)*3+m)-q((l-1)*3+m)
         enddo
         call crosprod_c2h6cn(rjk,rkl,zjkkl)
         call crosprod_c2h6cn(rij,rjk,zijjk)
         arg2 = anorm_c2h6cn(rjk)*dotprod_c2h6cn(rij,zjkkl)
         arg1 = dotprod_c2h6cn(zijjk,zjkkl)
c
c
         if (abs(arg1).lt.1.d-6) then
             if (arg2.gt.0.d0) dihed = 90.d0
             if (arg2.lt.0.d0) dihed = -90.d0
             return
         endif
         if (arg1.gt.0.d0) then
                 dihed = atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
         if (arg2.ge.0.d0) then
                 dihed = 180.d0 + atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
         if (arg2.lt.0.d0) then
                 dihed = -180.d0 + atan(arg2/arg1)
     *           *180.d0/pi
                 return
         endif
c
         return
         end
c
c------------------------------------------------------
       subroutine crosprod_c2h6cn(x,y,z)
c------------------------------------------------------
c
         implicit none
         double precision x(3),y(3),z(3)
c
         z(1)=x(2)*y(3)-x(3)*y(2)
         z(2)=x(3)*y(1)-x(1)*y(3)
         z(3)=x(1)*y(2)-x(2)*y(1)
c
         return
         end
c
c------------------------------------------------------
         double precision function dotprod_c2h6cn(x,y)
c------------------------------------------------------
c
         implicit none
         integer k
         double precision x(3),y(3), tmp
c
         tmp = 0.d0
         do k = 1, 3
            tmp = tmp + x(k)*y(k)
         enddo
         dotprod_c2h6cn = tmp
c
         return
         end
c
c------------------------------------------------------
         double precision function anorm_c2h6cn(x)
c------------------------------------------------------
         implicit none
         integer k
         double precision x(3), tmp
c
         tmp = 0.d0
         do k = 1, 3
            tmp = tmp + x(k)*x(k)
         enddo
         anorm_c2h6cn = dsqrt(tmp)
c
         return
         end
c
c------------------------------------------------------
       subroutine refangles_c2h6cn
c------------------------------------------------------
c
c  subroutine calculates reference angles for the "in-plane" potential
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
c      dimension sumd2(4),sumd4(4),ddr(4,4)
       tau=acos(-1.0d0/3.0d0)
c
c CRD 2012 taucc para el angulo H-C-C
c
       taucc=acos(-0.3616d0)
c
       pi=4.0d0*atan(1.0d0)
       halfpi=0.5d0*pi
       twopi=2.0d0*pi
c jcc-2010
c valor para ajustar el angulo
c      anghch=0.649d0*pi
c      tausih = pi - asin ( sin(anghch/2.d0) / sin(pi/3.d0) )
c      tausih=acos(-0.4648d0)
c jcc-2010
c
c
c  set diagonal elements to zero
c
       do i=1,8
         theta0(i,i)=0.0d0
       enddo
c
c  calculate reference angles
c
c      theta0(1,2)=tau+(tau-tausih)*(sphi(1)*sphi(2)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
c      theta0(1,3)=tau+(tau-tausih)*(sphi(1)*sphi(3)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
c      theta0(1,4)=taucc+(tau-tausih)*(sphi(1)*sphi(4)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
c      theta0(2,3)=tau+(tau-tausih)*(sphi(2)*sphi(3)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
c      theta0(2,4)=taucc+(tau-tausih)*(sphi(2)*sphi(4)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
c      theta0(3,4)=taucc+(tau-tausih)*(sphi(3)*sphi(4)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c      theta0(5,6)=tau+(tau-tausih)*(sphi(5)*sphi(6)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(7)*stheta(8)-1.0d0)
c      theta0(5,7)=tau+(tau-tausih)*(sphi(5)*sphi(7)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(6)*stheta(8)-1.0d0)
c      theta0(5,8)=taucc+(tau-tausih)*(sphi(5)*sphi(8)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(6)*stheta(7)-1.0d0)
c      theta0(6,7)=tau+(tau-tausih)*(sphi(6)*sphi(7)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(8)-1.0d0)
c      theta0(6,8)=taucc+(tau-tausih)*(sphi(6)*sphi(8)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(7)-1.0d0)
c      theta0(7,8)=taucc+(tau-tausih)*(sphi(7)*sphi(8)-1.0d0)
c    *             +(tau-anghch/3.0d0)*(stheta(5)*stheta(6)-1.0d0)
c
c el atomo 4 es el c2
c
       theta0(1,2)=tau+(tau-halfpi)*(sphi(1)*sphi(2)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(3)*stheta(4)-1.0d0)
       theta0(1,3)=tau+(tau-halfpi)*(sphi(1)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(4)-1.0d0)
       theta0(1,4)=taucc+(tau-halfpi)*(sphi(1)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(2)*stheta(3)-1.0d0)
       theta0(2,3)=tau+(tau-halfpi)*(sphi(2)*sphi(3)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(4)-1.0d0)
       theta0(2,4)=taucc+(tau-halfpi)*(sphi(2)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(3)-1.0d0)
       theta0(3,4)=taucc+(tau-halfpi)*(sphi(3)*sphi(4)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(1)*stheta(2)-1.0d0)
c
c  el atomo 8 es el c1
c
       theta0(5,6)=tau+(tau-halfpi)*(sphi(5)*sphi(6)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(7)*stheta(8)-1.0d0)
       theta0(5,7)=tau+(tau-halfpi)*(sphi(5)*sphi(7)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(6)*stheta(8)-1.0d0)
       theta0(5,8)=taucc+(tau-halfpi)*(sphi(5)*sphi(8)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(6)*stheta(7)-1.0d0)
       theta0(6,7)=tau+(tau-halfpi)*(sphi(6)*sphi(7)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(5)*stheta(8)-1.0d0)
       theta0(6,8)=taucc+(tau-halfpi)*(sphi(6)*sphi(8)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(5)*stheta(7)-1.0d0)
       theta0(7,8)=taucc+(tau-halfpi)*(sphi(7)*sphi(8)-1.0d0)
     *             +(tau-twopi/3.0d0)*(stheta(5)*stheta(6)-1.0d0)
c
c  fill in the other half of the matrix
c
        do i=1,3
          do j=i+1,4
            theta0(j,i)=theta0(i,j)
          enddo
        enddo
        do i=5,7
          do j=i+1,8
            theta0(j,i)=theta0(i,j)
          enddo
        enddo
c      do i=1,3
c         do j=i+1,4
c          write(*,*) 'theta',i,j,'=', theta0(i,j)*180/pi
c         enddo
c      enddo
c      do i=5,7
c         do j=i+1,8
c          write(*,*) 'theta',i,j,'=', theta0(i,j)*180/pi
c         enddo
c      enddo
       return
       end
c
c******************************************************
c
c
c------------------------------------------------------
       subroutine stretch_c2h6cn(vstr)
c------------------------------------------------------
c
c  subroutine to calculate leps-type stretching potential and its
c  derivatives
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
       dimension vqbh(6),vjbh(6),vq(3),vj(3),
     *           achdc(3),achdh(4,3),achdx(3), factj(4)
     *           ,vqc1h(3), vqc2h(3), vjc1h(3), vjc2h(3)
c
c  calculate avergage bond length for the methane moiety
c
       rav1=(rc1h(1)+rc1h(2)+rc1h(3))/3.0d0
       rav2=(rc2h(1)+rc2h(2)+rc2h(3))/3.0d0
c
c  initialise:
c
       vstr=0.0d0
c
c  ach:
c
c  nb: in double precision tanh(19.0d0)=1.0d0 and we put the if statement
c  in to avoid overflow/underflow errors
c
       arga1=c1c1h*(rav1-r0c1h)
       arga2=c1c2h*(rav2-r0c2h)
       if(arga1.lt.19.0d0)then
         ach1=a1c1h+b1c1h*(tanh(arga1)+1.0d0)*0.5d0
       else
         ach1=a1c1h+b1c1h
       endif
       if(arga2.lt.19.0d0)then
         ach2=a1c2h+b1c2h*(tanh(arga2)+1.0d0)*0.5d0
       else
         ach2=a1c2h+b1c2h
       endif
c      write(6,*) r0c1h, rav1, r0c2h, rav2, ach1, ach2
c
c      argac=c1cc*(rcc-r0cc)
c      if(argac.lt.19.0d0)then
c        acc=a1cc+b1cc*(tanh(argac)+1.0d0)*0.5d0
c      else
c        acc=a1cc+b1cc
c      endif
c
c  calculate singlet: e1, triplet: e3 energies and vq and vj
c  terms for each bond
c
       e1=d1cb*(exp(-2.0d0*acb*(rc1b-r0cb))-2.0d0*exp(-acb*(rc1b-r0cb)))
       e3=d3cb*(exp(-2.0d0*acb*(rc1b-r0cb))+2.0d0*exp(-acb*(rc1b-r0cb)))
       vqc1b=(e1+e3)*0.5d0
       vjc1b=(e1-e3)*0.5d0
c
       e1=d1cb*(exp(-2.0d0*acb*(rc2b-r0cb))-2.0d0*exp(-acb*(rc2b-r0cb)))
       e3=d3cb*(exp(-2.0d0*acb*(rc2b-r0cb))+2.0d0*exp(-acb*(rc2b-r0cb)))
       vqc2b=(e1+e3)*0.5d0
       vjc2b=(e1-e3)*0.5d0
c      write(*,*) vqc1b, vqc2b, vjc1b, vjc2b
c
c      e1=d1cc*(exp(-2.0d0*acc*(rcc-r0cc))-2.0d0*exp(-acc*(rcc-r0cc)))
c      e3=d3cc*(exp(-2.0d0*acc*(rcc-r0cc))+2.0d0*exp(-acc*(rcc-r0cc)))
c      vqcc=(e1+e3)*0.5d0
c      vjcc=(e1-e3)*0.5d0
c
       do i=1,3
         e1=d1c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h))
     *              -2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
         e3=d3c1h*(exp(-2.0d0*ach1*(rc1h(i)-r0c1h))
     *              +2.0d0*exp(-ach1*(rc1h(i)-r0c1h)))
         vqc1h(i)=(e1+e3)*0.5d0
         vjc1h(i)=(e1-e3)*0.5d0
c
         e1=d1c2h*(exp(-2.0d0*ach2*(rc2h(i)-r0c2h))
     *              -2.0d0*exp(-ach2*(rc2h(i)-r0c2h)))
         e3=d3c2h*(exp(-2.0d0*ach2*(rc2h(i)-r0c2h))
     *              +2.0d0*exp(-ach2*(rc2h(i)-r0c2h)))
         vqc2h(i)=(e1+e3)*0.5d0
         vjc2h(i)=(e1-e3)*0.5d0
       enddo
c
       do i=1,6
         e1=d1hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              -2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         e3=d3hh*(exp(-2.0d0*ahh*(rbh(i)-r0hh))
     *              +2.0d0*exp(-ahh*(rbh(i)-r0hh)))
         vqbh(i)=(e1+e3)*0.5d0
         vjbh(i)=(e1-e3)*0.5d0
       enddo
c
c  calculate 3 body potential
c
       do i=1,3
         vq(i)=vqc1h(i)+vqc1b+vqbh(i)
         vj(i)=-sqrt(((vjc1h(i)-vjc1b)**2+(vjc1b-vjbh(i))**2
     *                 +(vjbh(i)-vjc1h(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
c
       do i=1,3
         vq(i)=vqc2h(i)+vqc2b+vqbh(i+3)
         vj(i)=-sqrt(((vjc2h(i)-vjc2b)**2+(vjc2b-vjbh(i+3))**2
     *                 +(vjbh(i+3)-vjc2h(i))**2)*0.5d0)
         vstr=vstr+vq(i)+vj(i)
       enddo
c
c CRD 2012 potencial morse para el enlace C-C
c
       dt=(rcc-r0cc)
       expterm=exp(-acc*dt)
       vcc=d1cc*(1.0d0-expterm)**2.0d0
       vstr=vstr + vcc
c
c JEG 2023
c
c C-N  new term (simple morse term)
c
       dt=(rno-1.172d0)
       expterm=exp(-0.80d0*dt)
       vno=80.0d0*(1.d0-expterm)**2.d0
c      dt=(rno-r0hh)
c      expterm=exp(-ahh*dt)
c      vno=d1hh*(1.d0-expterm)**2.d0
c
       vstr=vstr+vno
c fin CRD 2019
c
       return
       end
c
c******************************************************
c
c------------------------------------------------------
       subroutine opbend_c2h6cn(vop)
c------------------------------------------------------
c
c  subroutine calculates symmetrized vop potential and derivatives
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
       double precision norma
       dimension sumd2(8),sumd4(8)
       dimension in(3),a(3),b(3),axb(3),c(8,3),argd(8)
c
c
       vop=0.0d0
c
c  calculate force constants and their derivatives
c
       call opforce
c
c  calculate out-of-plane angle and derivatives
c
c CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C2
c
c      do i=1,4
       i=4
c
c      do i=1,3
         j=i+1
         if(j.gt.4)j=j-4
         k=j+1
         if(k.gt.4)k=k-4
         l=k+1
         if(l.gt.4)l=l-4
c
         call calcdelta1_c2h6cn(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
c      write(*,*) i,j,k,l,sumd2(i)*fdelta(i),sumd4(i)*hdelta(i),vop
c      enddo
c
c CRD 2012 cambiamos el bucle a 3 para no incluir que salga el C1
c
c      do i=5,8
       i=8
c
c      do i=4,6
         j=i+1
         if(j.gt.8)j=j-4
         k=j+1
         if(k.gt.8)k=k-4
         l=k+1
         if(l.gt.8)l=l-4
c
c        j=1
c        k=2
c        l=3
c
         call calcdelta2_c2h6cn(i,j,k,l,sum2,sum4)
         sumd2(i)=sum2
         sumd4(i)=sum4
         vop=vop+fdelta(i)*sumd2(i)+hdelta(i)*sumd4(i)
c      write(*,*) i,j,k,l,sumd2(i)*fdelta(i),sumd4(i)*hdelta(i),vop
c      enddo
c
       return
       end
c
c
c******************************************************
c
c
c------------------------------------------------------
       subroutine ipbend_c2h6cn(vip)
c------------------------------------------------------
c
c  subroutine calculates symmetrised in plane bend term
c  and its derivatives
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
       dimension costh(8,8),theta(8,8),dth(8,8)
c
c  initialise
c
       vip=0.0d0
c
c  calculate force constants: fk0(i,j), f1(i)
c  and derivatives wrt rch(k) and rbh(k): dfdc(i,j,k), dfdh(i,j,k)
c
       call ipforce_c2h6cn
c
c  calculate theta(i,j) and in plane bend potential
c
c
c bucle para el c1
c
       do i=1,3
         do j=i+1,4
           if(j.lt.4.0) then
           costh(i,j)=tc1h(i,1)*tc1h(j,1)+tc1h(i,2)*tc1h(j,2)
     *                       +tc1h(i,3)*tc1h(j,3)
           costh(i,j)=costh(i,j)/rc1h(i)/rc1h(j)
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
         else
           costh(i,j)=tc1h(i,1)*tcc(1)+tc1h(i,2)*tcc(2)
     *                       +tc1h(i,3)*tcc(3)
           costh(i,j)=costh(i,j)/rc1h(i)/rcc
           theta(i,j)=acos(costh(i,j))
           dth(i,j)=theta(i,j)-theta0(i,j)
           vip=vip+0.5d0*fk0(i,j)*f1(i)*f1(j)*dth(i,j)**2
         endif
       enddo
       enddo
c
c bucle para el c2 para los theta 5,6,7,8
c
       do i=1,3
         do j=i+1,4
           if(j.lt.4.0) then
           costh(i+4,j+4)=tc2h(i,1)*tc2h(j,1)+tc2h(i,2)*tc2h(j,2)
     *                       +tc2h(i,3)*tc2h(j,3)
           costh(i+4,j+4)=costh(i+4,j+4)/rc2h(i)/rc2h(j)
           theta(i+4,j+4)=acos(costh(i+4,j+4))
           dth(i+4,j+4)=theta(i+4,j+4)-theta0(i+4,j+4)
           vip=vip+0.5d0*fk0(i+4,j+4)*f1(i+4)*f1(j+4)*dth(i+4,j+4)**2
         else
c modificacion para corregir los angulos del c2
c          costh(i,j)=tc2h(i,1)*tcc(1)+tc2h(i,2)*tcc(2)
c    *                       +tc2h(i,3)*tcc(3)
           costh(i+4,j+4)=tc2h(i,1)*(-tcc(1))+tc2h(i,2)*(-tcc(2))
     *                       +tc2h(i,3)*(-tcc(3))
           costh(i+4,j+4)=costh(i+4,j+4)/rc2h(i)/rcc
           theta(i+4,j+4)=acos(costh(i+4,j+4))
           dth(i+4,j+4)=theta(i+4,j+4)-theta0(i+4,j+4)
           vip=vip+0.5d0*fk0(i+4,j+4)*f1(i+4)*f1(j+4)*dth(i+4,j+4)**2
        endif
       enddo
       enddo
c
c CRD 2019
c
c
c  Now, calculate h2o bending energies.
c
c  First, calculate the angles:
c
       do i = 1,6
         dot = 0.d0
         do j = 1,3
           dot = dot - tno(j)*tbh(i,j)
         enddo
         cosine = dot / (rno*rbh(i))
         cosine = min(1.d0, max(-1.d0, cosine))
         angh2o(i) = acos(cosine)
       enddo
c
c      write (*,*) "angulos:"
c      do i=1,4
c      write (*,*) i, angh2o(i)*360.d0/(2.d0*3.1415926d0)
c      enddo
c
c  Now, calculate each force constant fkh2o(i) as a function of the O-H(I)
c  distance, rbh(i)
c
       do i=1,6
          arga= alph2o* (rbh(i) - r0hh)
          if(arga.lt.19.0d0)then
            fkh2o(i) = fkh2oeq * (1 - tanh(arga))
          else
            fkh2o(i) = 0.0d0
          endif
          write (95+i,*) rbh(i),fkh2o(i),fkh2oeq*(1 - tanh(arga))
       enddo
c
c  Now calculate the contribution to the energy for each H(I)-O-H(O)
c  harmonic bending term
c
      do i=1,6
           dang = (angh2o(i) - angh2oeq)
           vip=vip+0.5d0* fkh2o(i) * dang * dang
      enddo
c
c fin CRD 2019
c
       return
       end
c
c*************************************************************************
c
c------------------------------------------------------
       subroutine calcdelta1_c2h6cn(i,j,k,l,sum2,sum4)
c------------------------------------------------------
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
       include 'c2h6cn.inc'
c
       dimension  delta(8),in(3),a(3),b(3),axb(3),c(8,3),argd(8),
     *            daxb(8,3,3),cdot(8,3,3),atemp2(3)
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
c CRD 2005 incluimos bucles if para asegurar que siempre
c que una de los valores j, k o l sea 4 las coordenadas
c corresponden al atomo nuevo.
c
c  vector a is rk-rj, vector b is rl-rj
c
       if(j.eq.4)then
        do ind=1,3
          a(ind)=q(nh(k,ind))-q(nc2(ind))
          b(ind)=q(nh(l,ind))-q(nc2(ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
       if(k.eq.4)then
        do ind=1,3
          a(ind)=q(nc2(ind))-q(nh(j,ind))
          b(ind)=q(nh(l,ind))-q(nh(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
       if(l.eq.4)then
        do ind=1,3
          a(ind)=q(nh(k,ind))-q(nh(j,ind))
          b(ind)=q(nc2(ind))-q(nh(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
c
c modificacion para que hiciese vxu
c      axb(1)=b(2)*a(3)-b(3)*a(2)
c      axb(2)=a(1)*b(3)-b(1)*a(3)
c      axb(3)=b(1)*a(2)-a(1)*b(2)
c
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
c CRD 2012 eliminamos porque seria la salida del C2
c
       if(k.lt.4.and.j.lt.4.and.l.lt.4)then
        do ind=1,3
          a(ind)=q(nh(k,ind))-q(nh(j,ind))
          b(ind)=q(nh(l,ind))-q(nh(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
       else
       endif
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
       if(in(ii).eq.4) then
         do ind=1,3
           c(in(ii),ind)=-tcc(ind)/rcc
         enddo
       else
         do ind=1,3
           c(in(ii),ind)=-tc1h(in(ii),ind)/rc1h(in(ii))
         enddo
       endif
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
       return
       end

c******************************************************
c
c------------------------------------------------------
       subroutine calcdelta2_c2h6cn(i,j,k,l,sum2,sum4)
c------------------------------------------------------
c
c  subroutine calculates out of plane angle delta, loops
c  through delta(i,j), delta(i,k), delta(i,l)
c
c   also calculates the derivatives wrt delta
c
       implicit double precision (a-h,o-z)
       double precision norma
       include 'c2h6cn.inc'
c
       dimension  delta(8),in(3),a(3),b(3),axb(3),c(8,3),argd(8),
     *            daxb(8,3,3),cdot(8,3,3),atemp2(3)
c
c  initialise
c
       sum2=0.0d0
       sum4=0.0d0
c
c  set j,k,l indices
c
c      in(1)=j
c      in(2)=k
c      in(3)=l
       in(1)=4
       in(2)=5
       in(3)=6
c
c CRD 2005 incluimos bucles if para asegurar que siempre
c que una de los valores j, k o l sea 4 las coordenadas
c corresponden al atomo nuevo.
c
c  vector a is rk-rj, vector b is rl-rj
c
c      if(j.eq.8)then
c       do ind=1,3
c         a(ind)=q(nh(k,ind))-q(nc1(ind))
c         b(ind)=q(nh(l,ind))-q(nc1(ind))
c       enddo
c      axb(1)=a(2)*b(3)-a(3)*b(2)
c      axb(2)=a(3)*b(1)-a(1)*b(3)
c      axb(3)=a(1)*b(2)-a(2)*b(1)
c      norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
c      norma=sqrt(norma)
c      else
c      endif
c
c      if(k.eq.8)then
c       do ind=1,3
c         a(ind)=q(nc1(ind))-q(nh(j,ind))
c         b(ind)=q(nh(l,ind))-q(nh(j,ind))
c       enddo
c      axb(1)=a(2)*b(3)-a(3)*b(2)
c      axb(2)=a(3)*b(1)-a(1)*b(3)
c      axb(3)=a(1)*b(2)-a(2)*b(1)
c      norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
c      norma=sqrt(norma)
c      else
c      endif
c
c      if(l.eq.8)then
c       do ind=1,3
c         a(ind)=q(nh(k,ind))-q(nh(j,ind))
c         b(ind)=q(nc1(ind))-q(nh(j,ind))
c       enddo
c      axb(1)=a(2)*b(3)-a(3)*b(2)
c      axb(2)=a(3)*b(1)-a(1)*b(3)
c      axb(3)=a(1)*b(2)-a(2)*b(1)
c      norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
c      norma=sqrt(norma)
c      else
c      endif
c
c CRD 2012 eliminamos porque seria la salida del C1
c
c      if(k.lt.8.and.j.lt.8.and.l.lt.8)then
        do ind=1,3
          a(ind)=q(nh(k,ind))-q(nh(j,ind))
          b(ind)=q(nh(l,ind))-q(nh(j,ind))
        enddo
       axb(1)=a(2)*b(3)-a(3)*b(2)
       axb(2)=a(3)*b(1)-a(1)*b(3)
       axb(3)=a(1)*b(2)-a(2)*b(1)
       norma=axb(1)*axb(1)+axb(2)*axb(2)+axb(3)*axb(3)
       norma=sqrt(norma)
c      else
c      endif
c
c  c is position vector of h(ii): calculate c(j),c(k),c(l)
c
       do ii=1,3
c      if(in(ii).eq.4) then
c        do ind=1,3
c          c(in(ii),ind)=-tcc(ind)/rcc
c        enddo
c      else
         do ind=1,3
           c(in(ii),ind)=-tc2h(in(ii),ind)/rc2h(in(ii))
         enddo
c      endif
       enddo
c
c  argd is the dot product axb dot c
c
       do ii=1,3
         argd(in(ii))=axb(1)*c(in(ii),1)+axb(2)*c(in(ii),2)
     *                                +axb(3)*c(in(ii),3)
         argd(in(ii))=argd(in(ii))/norma
         delta(in(ii))=acos(argd(in(ii)))-theta0(i,in(ii))
         sum2=sum2+delta(in(ii))**2
         sum4=sum4+delta(in(ii))**4
       enddo
c
       return
       end
c******************************************************
c
c------------------------------------------------------
       subroutine opforce_c2h6cn
c------------------------------------------------------
c
c  calculates the out-of-plane bending force constants
c  and their derivatives
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
       dimension switch(8)
c
c  calculate switching functions:
c
       switch(1)=(1.0d0-s3(1))*s3(2)*s3(3)*s3(4)
       switch(2)=(1.0d0-s3(2))*s3(3)*s3(4)*s3(1)
       switch(3)=(1.0d0-s3(3))*s3(4)*s3(1)*s3(2)
       switch(4)=(1.0d0-s3(4))*s3(1)*s3(2)*s3(3)
       switch(5)=(1.0d0-s3(5))*s3(6)*s3(7)*s3(8)
       switch(6)=(1.0d0-s3(6))*s3(7)*s3(8)*s3(5)
       switch(7)=(1.0d0-s3(7))*s3(8)*s3(5)*s3(6)
       switch(8)=(1.0d0-s3(8))*s3(5)*s3(6)*s3(7)
c      do i=1,8
c       write(6,*) 'switch',i,'=',switch(i)
c      enddo
c
c  calculate the force constants
c
       do i=1,8
         fdelta(i)=switch(i)*fch3
         hdelta(i)=switch(i)*hch3
       enddo
c
       return
       end
c
c******************************************************
c
c
c------------------------------------------------------
       subroutine ipforce_c2h6cn
c------------------------------------------------------
c
c  calculates the symmetrised in plane bend force constants and
c  all partial derivatives involving them
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
c  set force constant at asymptotes
c
       f0=fkinf+ak
       f2=fkinf
c
       fk0(1,2)=f0+f0*(s1(1)*s1(2)-1.0d0)+(f0-f2)*(s2(3)*s2(4)-1.0d0)
       fk0(1,3)=f0+f0*(s1(1)*s1(3)-1.0d0)+(f0-f2)*(s2(2)*s2(4)-1.0d0)
       fk0(1,4)=f0+f0*(s1(1)*s1(4)-1.0d0)+(f0-f2)*(s2(2)*s2(3)-1.0d0)
       fk0(2,3)=f0+f0*(s1(2)*s1(3)-1.0d0)+(f0-f2)*(s2(1)*s2(4)-1.0d0)
       fk0(2,4)=f0+f0*(s1(2)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(3)-1.0d0)
       fk0(3,4)=f0+f0*(s1(3)*s1(4)-1.0d0)+(f0-f2)*(s2(1)*s2(2)-1.0d0)
c
       fk0(5,6)=f0+f0*(s1(5)*s1(6)-1.0d0)+(f0-f2)*(s2(7)*s2(8)-1.0d0)
       fk0(5,7)=f0+f0*(s1(5)*s1(7)-1.0d0)+(f0-f2)*(s2(6)*s2(8)-1.0d0)
       fk0(5,8)=f0+f0*(s1(5)*s1(8)-1.0d0)+(f0-f2)*(s2(6)*s2(7)-1.0d0)
       fk0(6,7)=f0+f0*(s1(6)*s1(7)-1.0d0)+(f0-f2)*(s2(5)*s2(8)-1.0d0)
       fk0(6,8)=f0+f0*(s1(6)*s1(8)-1.0d0)+(f0-f2)*(s2(5)*s2(7)-1.0d0)
       fk0(7,8)=f0+f0*(s1(7)*s1(8)-1.0d0)+(f0-f2)*(s2(5)*s2(6)-1.0d0)
c
c f1 para los h del c1
c
       do i=1,3
c
c  calc derivatives of fk0 wrt each of the rch(i) bonds
c  calculate the terms f1(i)
c
         arga1=aa1*rbh(i)*rbh(i)
         arga2=aa4*(rbh(i)-r0hh)*(rbh(i)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rc1h(i)-r0c1h)**2)
       enddo
c
c f1 para el c2
c
c CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
c lo incluimos en el CONST
c
c        arga1c=aa1*rbc*rbc
c        arga2c=aa4*(rbc-r0cb)*(rbc-r0cb)
c        a1c=1.0d0-exp(-arga1c)
c        a2c=aa2+aa3*exp(-arga2c)
         f1(4)=a1c*exp(-a2c*(rcc-r0cc)**2)
         f1(8)=f1(4)
c
c f1 para los h del c2
c
       do i=5,7
c
c  calc derivatives of fk0 wrt each of the rch(i) bonds
c  calculate the terms f1(i)
c
         arga1=aa1*rbh(i-1)*rbh(i-1)
         arga2=aa4*(rbh(i-1)-r0hh)*(rbh(i-1)-r0hh)
         a1=1.0d0-exp(-arga1)
         a2=aa2+aa3*exp(-arga2)
         f1(i)=a1*exp(-a2*(rc2h(i-4)-r0c2h)**2)
       enddo
c
c f1 para el c2
c
c CRD 2012 eliminamos el calculo de a1c y a2c pq es constante
c lo incluimos en el CONST
c
c        arga1c=aa1*rcb*rcb
c        arga2c=aa4*(rcb-r0cb)*(rcb-r0cb)
c        a1c=1.0d0-exp(-arga1c)
c        a2c=aa2+aa3*exp(-arga2c)
c        f1(8)=a1c*exp(-a2c*(rcc-r0cc)**2)
c
       return
       end
c
c******************************************************
c
c------------------------------------------------------
       subroutine switchf_c2h6cn
c------------------------------------------------------
c
c  calculates switching functions: s3,sphi,stheta
c  and their derivatives ds3,dsphi,dstheta
c
       implicit double precision (a-h,o-z)
       include 'c2h6cn.inc'
c
c  nb remember that integration units are:
c  energy in   1.0d+05 j/mol
c  time in     1.0d-14 s
c
c original     a1s=1.5313681d-7
c      a1s=1.5313681d-8
c      b1s=-4.6696246d0
c original      a2s=1.0147402d-7
c      a2s=1.0147402d-8
c      b2s=-12.363798d0
c
c  use double precision criterion:
c
c  tanh(19.0d0)=1.0d0
c
       argmax=19.0d0
c
c  calculate s1 and ds1
c
       do i=1,3
         args1=a1s*(rc1h(i)-r0c1h)*(rc1h(i)-b1s)**8
         if(args1.lt.argmax)then
           s1(i)=1.0d0-tanh(args1)
         else
           s1(i)=0.0d0
         endif
       enddo
c
       do i=1,3
         args1=a1s*(rc2h(i)-r0c2h)*(rc2h(i)-b1s)**8
         if(args1.lt.argmax)then
           s1(i+4)=1.0d0-tanh(args1)
         else
           s1(i+4)=0.0d0
         endif
       enddo
c
c CRD 2012 todos los S de los carbonos son 1
c
c      s1(4)=1.0d0
c      s1(8)=1.0d0
         argsc1=a1scc*(rcc-r0cc)*(rcc-b1scc)**8
         if(argsc1.lt.argmax)then
           s1(4)=1.0d0-tanh(argsc1)
           s1(8)=s1(4)
         else
           s1(4)=0.0d0
           s1(8)=s1(4)
         endif
c
c  calculate s2 and ds2
c
       do i=1,3
         args2=a2s*(rc1h(i)-r0c1h)*(rc1h(i)-b2s)**6
         if(args2.lt.argmax)then
           s2(i)=1.0d0-tanh(args2)
         else
           s2(i)=0.0d0
         endif
       enddo
c
       do i=1,3
         args2=a2s*(rc2h(i)-r0c2h)*(rc2h(i)-b2s)**6
         if(args2.lt.argmax)then
           s2(i+4)=1.0d0-tanh(args2)
         else
           s2(i+4)=0.0d0
         endif
       enddo
c
c CRD 2005 s2(4)for the new atom
c
c CRD 2012 todos los S de los carbonos son 1
c
c      s2(4)=1.0d0
c      s2(8)=1.0d0
         argsc2=a2scc*(rcc-r0cc)*(rcc-b2scc)**6
         if(argsc2.lt.argmax)then
           s2(4)=1.0d0-tanh(argsc2)
           s2(8)=s2(4)
         else
           s2(4)=0.0d0
           s2(8)=s2(4)
         endif
c
c  calculate s3 and ds3
c
       do i=1,3
         args3=a3s*(rc1h(i)-r0c1h)*(rc1h(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i)=1.0d0-tanh(args3)
         else
           s3(i)=0.0d0
         endif
       enddo
c
       do i=1,3
         args3=a3s*(rc2h(i)-r0c2h)*(rc2h(i)-b3s)**2
         if (args3.lt.argmax)then
           s3(i+4)=1.0d0-tanh(args3)
         else
           s3(i+4)=0.0d0
         endif
       enddo
c
c CRD 2005 s3(4) for the new atom
c
c CRD 2012 todos los S de los carbonos son 1
c
c      s3(4)=1.0d0
c      s3(8)=1.0d0
         argsc3=a3scc*(rcc-r0cc)*(rcc-b3scc)**2
         if (argsc3.lt.argmax)then
           s3(4)=1.0d0-tanh(argsc3)
           s3(8)=s3(4)
         else
           s3(4)=0.0d0
           s3(8)=s3(4)
         endif
c CRD 2019 fijamos el S3
       s3(4)=1.0d0
       s3(8)=1.0d0
c
c  calculate sphi and dsphi
c
c  condition here is on the bondlength rch(i)
c  st argsphi is lt approx 19.0d0
c
       do i=1,3
         if(rc1h(i).lt.3.8d0)then
           argsphi=aphi*(rc1h(i)-r0c1h)*exp(bphi*(rc1h(i)-cphi)**3)
           sphi(i)=1.0d0-tanh(argsphi)
         else
           sphi(i)=0.0d0
         endif
       enddo
c
       do i=1,3
         if(rc2h(i).lt.3.8d0)then
           argsphi=aphi*(rc2h(i)-r0c2h)*exp(bphi*(rc2h(i)-cphi)**3)
           sphi(i+4)=1.0d0-tanh(argsphi)
         else
           sphi(i+4)=0.0d0
         endif
       enddo
c
c CRD 2012 todos los S de los carbonos son 1
c
c      sphi(4)=1.0d0
c      sphi(8)=1.0d0
         if(rcc.lt.3.8d0)then
           argsphic=aphi*(rcc-r0cc)*exp(bphi*(rcc-cphi)**3)
           sphi(4)=1.0d0-tanh(argsphic)
           sphi(8)=sphi(4)
         else
           sphi(4)=0.0d0
           sphi(8)=sphi(4)
         endif
c
c  calculate stheta and dstheta
c
       do i=1,3
         if(rc1h(i).lt.3.8d0)then
       argstheta=atheta*(rc1h(i)-r0c1h)*exp(btheta*(rc1h(i)-ctheta)**3)
           stheta(i)=1.0d0-tanh(argstheta)
         else
           stheta(i)=0.0d0
         endif
       enddo
c
       do i=1,3
         if(rc2h(i).lt.3.8d0)then
       argstheta=atheta*(rc2h(i)-r0c2h)*exp(btheta*(rc2h(i)-ctheta)**3)
           stheta(i+4)=1.0d0-tanh(argstheta)
         else
           stheta(i+4)=0.0d0
         endif
       enddo
c
c CRD 2012 todos los S de los carbonos son 1
c
c      stheta(4)=1.0d0
c      stheta(8)=1.0d0
         if(rcc.lt.3.8d0)then
           argsthetac=atheta*(rcc-r0cc)*exp(btheta*(rcc-ctheta)**3)
           stheta(4)=1.0d0-tanh(argsthetac)
           stheta(8)=stheta(4)
         else
           stheta(4)=0.0d0
           stheta(8)=stheta(4)
         endif
       return
       end
C

Cc---------------------------------------------------------------------
C      subroutine vander(evdw)
Cc---------------------------------------------------------------------
Cc
Cc  subroutine to calculate van der waals potential
Cc
C       implicit double precision (a-h,o-z)
C       include 'c2h6cn.inc'
Cc
Cc original de jose carlos
Cc
CCC       rav=(rch(1)+rch(2)+rch(3)+rch(4))/4.0d0
CCC       p1vdw = vdwc1*(vdwr/rcb)**6.d0
CCC       p2vdw = vdwc2*exp(-12.d0*(rcb/vdwr))
CCC       evdw1 = vdwe/0.03812d0/627.51 * (p1vdw+p2vdw)
CCC
CCC       devdw = vdwe/0.03812d0/627.51 *
CCC     *        ( -6*p1vdw/rcb - 12*p2vdw/vdwr)
CCC
CCCc
CCCc      tvdw = 1.d0+vdt1*tanh(vdt2*(rav-vdtr))
CCCc
CCC       arga=(vdt2*(rav-vdtr))
CCC       if(arga.lt.19.0d0)then
CCC         tvdw = 1.d0+vdt1*tanh(arga)
CCC         dtvdw = vdt1*vdt2*(1/cosh(arga))**2.d0
CCC       else
CCC         tvdw=1.d0+vdt1
CCC         dtvdw=0.0d0
CCC       endif
CCC       evdw = evdw1 * tvdw
Cc
Cc
Cc  CRD 2013
Cc
Cc  calculate avergage bond length for the ethane
Cc
C       rav1=(rc1h(1)+rc1h(2)+rc1h(3))/3.0d0
C       rav2=(rc2h(1)+rc2h(2)+rc2h(3))/3.0d0
Ccc
Ccc Calculate for C1
Ccc
Cc       p1vdw = vdwc1*(vdwr/rc1b)**6.d0
Cc       p2vdw = vdwc2*exp(-12.d0*(rc1b/vdwr))
Cc       evdw1 = vdwe/0.03812d0/627.51*(p1vdw+p2vdw)
Ccc
Cc       arga=(vdt2*(rav1-vdtr))
Ccc
Cc       if(arga.lt.19.0d0)then
Cc         tvdw = 1.d0+vdt1*tanh(arga)
Cc       else
Cc         tvdw=1.d0+vdt1
Cc       endif
Ccc
Cc       evdwc1 = evdw1 * tvdw
Ccc
Ccc Calculate for C2
Ccc
Cc       p1vdw = vdwc1*(vdwr/rc2b)**6.d0
Cc       p2vdw = vdwc2*exp(-12.d0*(rc2b/vdwr))
Cc       evdw1 = vdwe/0.03812d0/627.51*(p1vdw+p2vdw)
Ccc
Cc       arga=(vdt2*(rav2-vdtr))
Ccc
Cc       if(arga.lt.19.0d0)then
Cc         tvdw = 1.d0+vdt1*tanh(arga)
Cc       else
Cc         tvdw=1.d0+vdt1
Cc       endif
Ccc
Cc       evdwc2 = evdw1 * tvdw
Cc       evdw = evdwc1 + evdwc2
Cc       write(48,*) evdwc1,evdwc2,evdw*0.03812d0*627.52
Cc
Cc CRD 2014 VDW para todos los H
Cc
C      evdw=0.0
C      do i=1,3
C       p1vdw = vdwc1*(vdwr/rc1h(i))**6.d0
C       p2vdw = vdwc2*exp(-12.d0*(rc1h(i)/vdwr))
C       evdw1 = vdwe/0.03812d0/627.51*(p1vdw+p2vdw)
Cc
C       arga=(vdt2*(rav1-vdtr))
Cc
C       if(arga.lt.19.0d0)then
C         tvdw = 1.d0+vdt1*tanh(arga)
C       else
C         tvdw=1.d0+vdt1
C       endif
Cc
C       evdwci = evdw1 * tvdw
C       evdw=evdw+evdwci
C      enddo
C      do i=1,3
C       p1vdw = vdwc1*(vdwr/rc2h(i))**6.d0
C       p2vdw = vdwc2*exp(-12.d0*(rc2h(i)/vdwr))
C       evdw1 = vdwe/0.03812d0/627.51*(p1vdw+p2vdw)
Cc
C       arga=(vdt2*(rav2-vdtr))
Cc
C       if(arga.lt.19.0d0)then
C         tvdw = 1.d0+vdt1*tanh(arga)
C       else
C         tvdw=1.d0+vdt1
C       endif
Cc
C       evdwci = evdw1 * tvdw
C       evdw=evdw+evdwci
C      enddo
Cc      write(47,*) evdw, rav1, rav2
Cc
Cc  CRD 2014 eliminamos la contribucion vander
Cc
Cc     evdw=0.0
Cc
C      return
C      end
c
       SUBROUTINE SETUP_c2h6cn(N3TM)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      PARAMETER (N3TMMN = 30)
C
C  CHECK THE NUMBER OF CARTESIAN COORDINATES SET BY THE CALLING PROGRAM
C
      WRITE (6, 1300)
      IF (N3TM .LT. N3TMMN) THEN
          WRITE (6, 6000) N3TM, N3TMMN
          STOP 'SETUP 1'
      ENDIF
C
      RETURN
C
1300  FORMAT(/,2X,T5,'SETUP has been called for the C2H6H ',
     *               'potential energy surface')
6000  FORMAT(/,2X,T5,'Warning: N3TM is set equal to ','9',
     *                  ' but this potential routine',
     *          /,2X,T14,'requires N3TM be greater than or ',
     *                   'equal to ','9',/)
C
      END
