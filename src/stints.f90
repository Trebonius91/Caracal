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
!    subroutine stints: calculate s/t integrals in cao(6d) basis
!    tm int version (what does that mean???)
!
!    part of QMDFF
!
subroutine stints(nat,nbf,xyz,s)  
use qmdff
implicit none 
integer::nat,nbf
real(kind=8)::xyz(3,nat)          
real(kind=4)::s(nbf,nbf)

integer::i,j,k,l,iprimcount,jprimcount
integer::npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn,n
integer::lll(20),iall(4,4)
integer::lin(84),min(84),nin(84)
integer::lmnexp(84),ib,ie
real(kind=8)::xyza(3),xyzb(3),rab,est,ss,sss,lmnfak(84),gama,arg
real(kind=8)::aa(10),bb(10),gm2,ttt(1),tt,intcut

intcut=25.0d0

lll=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)

lin=(/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0 &
  & ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3, &
  & 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/)
min=(/0,0,1,0,0, &
  & 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5, &
  & 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0, &
  & 0,2,4,4,2,1,3,1,3,2,1,4,1,2/)
nin=(/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0, &
  & 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1, &
  & 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3, &
  & 1,1,4,2/)

s=0

iall(1,1)=1
iall(1,2)=4
iall(2,1)=4
iall(1,3)=10
iall(3,1)=10
iall(2,2)=10
iall(2,3)=20
iall(3,2)=20
iall(2,4)=35
iall(4,2)=35
iall(3,3)=35
iall(3,4)=56
iall(4,3)=56
iall(4,4)=84

k=0
iprimcount=0
do i=1,nbf 
!
!     aufpunkt i (first coordinate)
!
   xyza(1:3)=xyz(1:3,aoat(i))
!     #prims
   npri=nprim(i)
   jprimcount=0
   li=lll(lao(i))
   do j=1,i
      lj=lll(lao(j))
      k=k+1
      nprj=nprim(j)
!     aufpunkt j (second coordinate)
      xyzb(1:3)=xyz(1:3,aoat(j))
      rab=(xyza(1)-xyzb(1))**2 &
        & +(xyza(2)-xyzb(2))**2 &
        & +(xyza(3)-xyzb(3))**2
!
!     precalc some overlap terms that depend only on lm
!
      do ll=1,iall(li,lj)
         call lmnpre(lin(ll),min(ll),nin(ll),lmnexp(ll),lmnfak(ll))
      enddo
      aa=0
      bb=0
      aa(lao(i))=1.0d0
      bb(lao(j))=1.0d0
!     prim loop
      ss=0.0d0
      do ii=1,npri
         iii=iprimcount+ii
         do jj=1,nprj
            jjj=jprimcount+jj
            gama=1.0d0/(alp(iii)+alp(jjj))
            gm2 =0.5d0*gama
            est=rab*alp(iii)*alp(jjj)*gama              
!    
!     cutoff criterion
!
            if (est.lt.intcut) then
               arg=(pi*gama)**1.50d0
               call pola(xyza,xyzb,alp(iii),alp(jjj), &
                  &        gama,gm2,lao(i),lao(j),iall(li,lj), &
                  &       aa,bb,lmnexp,lmnfak,est,arg,sss)
               ss=ss+sss*cont(iii)*cont(jjj)   
            end if
         end do
      end do
      s(i,j)=ss
      s(j,i)=ss
 42   jprimcount=jprimcount+nprj
   end do
   iprimcount=iprimcount+npri
end do
!
!    normalized?
!
do i=1,nbf 
!
!    no diagonal contribution to H0
!
   if(abs(1.d0-1.0d0/sqrt(s(i,i))).gt.1.d-6)then
      write(*,*) i,s(i,i)
      stop 'function not normalized inside stints'
   endif
enddo

end subroutine stints
