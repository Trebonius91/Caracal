       common /angles/ theta0(8,8),dtheta0(8,8,8)
       common /bonds/ r0ch,r0hh,r0cb
       common /bonds2/ rc1h(3),rc2h(3),rbh(6), rhh(3,6)
       common /bonds3/ r0cc,rcc, rc1b, rc2b
       common /bonds4/ r0c1hr,r0c1hp,w3,w4,r0c1h
       common /bonds5/ r0c2hr,r0c2hp,w5,w6,r0c2h
       common /coords1/ tcb(3),tc1h(3,3),tc2h(3,3),tbh(6,3)
       common /coords2/ tcc(3),thh(3,6,3)
       common /coords3/ tc1b(3),tc2b(3)
       common /delta1/ fdelta(8),hdelta(8)
       common /delta2/ dfdelta(8,8),dhdelta(8,8)
       common /force1/ fk0(8,8),f1(8),dfdc(8,8,8),dfdh(8,8,8)
       common /force2/ fkinf,ak,bk,aa1,aa2,aa3,aa4
       common /fsw1/ a1s,b1s,a2s,b2s,a3s,b3s
       common /fsw2/ a1scc,b1scc,a2scc,b2scc,a3scc,b3scc
       common /ip1/ s1(8),ds1(8),s2(8),ds2(8)
       common /ndx/ nnb,nnh(6),nc(3),nhb(3),nh(6,3)
       common /ndx2/ nnc1, nnc2,nc1(3), nc2(3)
       common /op1/ s3(8),ds3(8)
       common /params/ fch3,hch3
       common /qpdot/ q(30),pdot(30)
       common /sphi1/ aphi,bphi,cphi
       common /stheta1/ atheta,btheta,ctheta
       common /stret1/ d1ch,d3ch,a1cc,b1cc,c1cc
       common /stret2/ d1hh,d3hh,ahh
       common /stret3/ d1cb,d3cb,acb
       common /stret4/ d1cc,acc, a1c, a2c
       common /stret5/ a1c1h, b1c1h, c1c1h
       common /stret6/ a1c2h, b1c2h, c1c2h 
       common /stret7/ d1c1h,d3c1h,d1c2h,d3c2h
       common /switch1/ sphi(8),dsphi(8),stheta(8),dstheta(8)
       common /thetahcch/ thetahcch(9), w1,w2,V3
       common /vdwpot/ vdwe, vdwr, vdwc1, vdwc2, vdt1, vdt2, vdtr
       common /h2o/ rno,tno(3),nno,no(3)
       common /addh2o/angh2o(6),fkh2oeq,fkh2o(6),alph2o,angh2oeq

