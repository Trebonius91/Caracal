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
!     subroutine ff_bond: Determine all bond terms to be used 
!      in the new QMDFF
!
!     part of QMDFF
! 
subroutine ff_bond(n,at,nb,xyz,pair,q,wbo,cn,molist)
use qmdff
implicit none
integer::n,nb(20,n),at(n),pair(n*(n+1)/2),molist(n)
real(kind=8)::xyz(3,n),wbo(n,n),q(n),cn(n)

real(kind=8)::r,rmin,rco,be(6),dum,dum2,aa,bb,cc,aaa,b1,b2,thr
integer::nn(n),list(1000,n),jj,nlist(1000,n),i2,k
integer::i,ni,newi,ii,i1,ni1,iii,d,j,newatom,nj,m,pqn
integer::imem,jmem,nnn(n),l,tag,irow,jrow,lina,jjj
logical::da,dai,daj,l1,l2,l3,metal13

nn(1:n)=nb(20,1:n)

ade = -1 
!
!     EMPIRICAL
!     fitted to G2 data with freq at tpss-d3/def2-tzvp(-f) 
!     value for be(4) and larger 
!     adjusted to De of Br2/I2/At2 (45.9,35.9,17.7 from  J. Chem. Phys. 104 (22), 1996)
!     values for group 1 and 2 adjusted to li2,na2,be2,mg2
!     MAD(G2)=22 kcal!
!

be(1)=1.755
be(2)=2.465
be(3)=2.75 
be(4)=2.95
be(5)=3.15 
be(6)=3.80
bb=-0.164

ade(1 )=1.755
ade(5 )=2.287  
ade(6 )=2.463
ade(7 )=2.559 
ade(8 )=2.579 
ade(9 )=2.465  
ade(13:17)=ade(5:9)+0.221

ade(3 )=2.2
ade(11)=2.2 
ade(19)=2.2 
ade(37)=2.2 
ade(55)=2.2 
ade(87)=2.2 

ade(4 )=2.8
ade(12)=2.8
ade(20)=2.8
ade(38)=2.8
ade(56)=2.8
ade(88)=2.8



!     open(unit=66,file='~/.param')
!     read(66,*) dum,dum2,ade(1),ade(5:9),cc,bb
!     read(66,*) ade13(1),ade13(5:9),cc,bb
!     close(66)
!     ade(13:17)=ade(5:9)+cc   
!     if ade_12 and ade_13 are made different for H,B-F one gets an G2 MAD
!     of 18.5 kcal with these parameters:
!     2.52236995  0.54446115  1.88661996  2.26125328  2.41006715  2.50378046  
!     2.50123260  2.60063104  0.17181602 -0.20460162
!     1.49705258  2.28621936  2.74602255  2.94826838  3.06517993  2.48384362

do i=1,94
   if (ade(i).lt.0) ade(i)=be(pqn(i))
enddo

ade13 = ade

metal13=.true. 

pair=0
!
!     List with neighbors?
!
list=0
do i=1,n
   ni=nn(i)
   list(1:ni,i)=nb(1:ni,i)
end do

nlist=list
!
!     one bond, tag=1      
!
tag=1
call pairs(n,nn,list,pair,tag)
!
!     determine up to 4 bonds in between      
!
do d=1,3
!
!     loop over all atoms
!      
   do i=1,n
      ni=nn(i)
      newi=ni
!
!     all neighbors of i (for bond terms)
!       
      do ii=1,ni
         i1=list(ii,i)
         ni1=nb(20,i1)
!
!     all neighbors of neighbors of i (for 1-3 interaction terms)  
!     
         do iii=1,ni1
            newatom=nb(iii,i1)
            da=.false.
            do j=1,newi
               if(newatom.eq.list(j,i))da=.true.
            end do
            if (.not.da) then
               newi=newi+1
               nlist(newi,i)=newatom           
            end if
         end do
      end do
      nnn(i)=newi
   end do
   list=nlist
   nn  =nnn
!
!     one bond more
!
   tag=tag+1
   call pairs(n,nn,list,pair,tag)
end do
!
!     pair=5 tags five or more covalent bonds
!
do i=1,n*(n+1)/2
   if(pair(i).eq.0) pair(i)=5
enddo
!
!     1,3 interaction involving a metal are made NCI
!     tag: nci(3,*)=6 (RE+DISP scaled, no ES)
!     so: NO METAL 1-3 BONDS POSSIBLE!!! 
!
if (metal13) then
   thr=0.20
   do i=1,n
      if (metal(at(i)).eq.1.or.nb(20,i).gt.6) then
         do ii=1,nb(20,i)
            iii=nb(ii,i)
            if (iii.eq.i) cycle
            do jj=1,nb(20,i)-1
               jjj=nb(jj,i)
               if (jjj.eq.j) cycle
               if (jjj.eq.i) cycle
               if (jjj.eq.iii) cycle
               k=lina(iii,jjj)
               l1=pair(k).eq.2.or.pair(k).eq.3
               l2=wbo(iii,jjj).lt.thr
!
!      1,3 have to be in different fragments                  
!      these are defined when the metal is removed
!
               l3=molist(iii).ne.molist(jjj)
               if (l1.and.l2.and.l3) then
                  pair(k)=6
                  write(10, & 
                   & '("setting 1,3-metal interaction to screened NCI",2i4)') &
                   & iii,jjj
               end if
            end do
         end do
      end if
   end do
end if

k=0
m=0
l=0
if (details) then
   open(unit=38,file="qmdff_details.dat",status="replace")
end if
do i=1,n
   do j=1,i   
      k=k+1
      if(i.eq.j) cycle
!
! the direct stretch       
!        
      if (pair(k).eq.1) then
         r=sqrt((xyz(1,i)-xyz(1,j))**2 &
           &   +(xyz(2,i)-xyz(2,j))**2 &
           &   +(xyz(3,i)-xyz(3,j))**2)*0.52917726
         l=l+1
         bond(1,l)=i
         bond(2,l)=j
         b1 =ade(at(i))
         b2 =ade(at(j))
         aaa=b1*b2+bb*(en(at(i))-en(at(j)))**2
         if (details) then
            write(38,*) "bond(i-j):",i,j,"k_a(i):",b1,"k_a(j):",b2,"k_EN:",bb,"EN(i):",en(at(i)),"EN(j):",en(at(j))
         end if

         vbond(3,l)=min(aaa,20.0d0)            
      else if (pair(k).gt.2) then
!
!     NCI if # cov. bonds > 2 = 1,4 interactions      
!     
         m=m+1
         nci(1,m)=i
         nci(2,m)=j
         nci(3,m)=pair(k)
      end if
   end do
end do
if (details) close (38)

nnci=m
nbond=l

write(10,'('' List of all directly bonded atoms: '')')
write(10,*)
write(10,*)'-------------------------------------------------'
write(10,'(''  #1   #2    pair distance   Rcov      WBO'')')
write(10,*)'-------------------------------------------------'
do l=1,nbond
   i=bond(1,l)
   j=bond(2,l)
   r=sqrt((xyz(1,i)-xyz(1,j))**2 &
      &  +(xyz(2,i)-xyz(2,j))**2 &
      &  +(xyz(3,i)-xyz(3,j))**2)*0.52917726
   rco=rad(at(i))+rad(at(j))
!   this should fixed in the future...
   write(10,'(2I4,3x,F8.3,A,F7.2,A,F7.2)') bond(1:2,l),r,"       ",rco,"  ",wbo(i,j)
!   write(10,*) bond(1:2,l),r,rco,wbo(i,j)
   if (wbo(i,j).lt.0) call warn('WBO<0. check orca.out!')
end do
write(10,*)'-------------------------------------------------'
write(10,*)

rmin=10000.
imem=0
do k=1,nnci
   i1=nci(1,k)
   i2=nci(2,k)
   r=sqrt((xyz(1,i1)-xyz(1,i2))**2 &
     &   +(xyz(2,i1)-xyz(2,i2))**2 &
     &   +(xyz(3,i1)-xyz(3,i2))**2)
   if (r.lt.rmin) then
      rmin=r
      imem=i1
      jmem=i2
   end if
end do

if (imem.gt.0) then
     write(10,'("smallest distance considered as Noncovalent interaction (NCI)")')
     write(10,'(" for #1, #2, r(bohr), r(Ang) :",2i4,2F8.2)')  imem,jmem,rmin,rmin*0.529167
end if

return
end subroutine ff_bond

