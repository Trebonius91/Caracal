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
!     subroutine eht: Extended Hückel calculation on torsional profile
!          for FF parameter optimization; similar to ehtfull
!
!     part of QMDFF
!
subroutine eht(n,at,xyz,q,nel,nbf,eel,gscal,list,pair)
use qmdff
implicit none
integer::n,at(n),nel,nbf,list(*),pair(*)
real(kind=8)::xyz(3,n),eel,wb(n,n),q(n),gscal
logical::dens

real(kind=8)::b,bb,temp,xsum,eh1,gab,rab,brow(6),hav,r
real(kind=4),allocatable::X(:,:)
real(kind=4),allocatable::P(:,:)
real(kind=4),allocatable::S(:,:)
real(kind=4),allocatable::H(:,:)
real(kind=4),allocatable::aux(:)
real(kind=4),allocatable::emo(:)
real(kind=4),allocatable::ep (:)
real(kind=4),allocatable::focc(:)
integer,allocatable::iwork(:)
integer,allocatable::iiwork(:)
integer::LWORK,LIWORK,INFO
integer::npr,ii,jj,i,j,kk,m,k,iat,jat,lina
integer::mmm(10),xi,xj,ihomo,irow,jrow,pqn
character(len=1)::flag
logical::centralbond

mmm=(/1,2,2,2,3,3,3,3,3,3/)
liwork = 3 + 5*nbf
lwork  = 1 + 6*nbf + 2*nbf**2
allocate(H(nbf,nbf),X(nbf,nbf),aux(lwork),P(nbf,nbf), &
  & iwork(liwork),emo(nbf),focc(nbf),ep(nbf),S(nbf,nbf))

call stints(n,nbf,xyz,X)
S = X
!
!     EMPIRICAL Tue Jun 10 21:39:28 CEST 2014
!
bb=1.75

do ii=1,nbf
   xi=aoat(ii)
   irow=pqn(at(xi))
   do jj=1,ii-1
     xj=aoat(jj)         
     jrow=pqn(at(xj))
     centralbond=xi+xj.eq.2               
     if (list(xi).ne.list(xj).and.(.not.centralbond)) then 
        X(jj,ii)=sscal(pair(lina(xi,xj)))*X(jj,ii)
     end if
     r=sqrt((xyz(1,xi)-xyz(1,xj))**2 &
       &  +(xyz(2,xi)-xyz(2,xj))**2 &
       &  +(xyz(3,xi)-xyz(3,xj))**2)
     hav=(hdiag(ii)+hdiag(jj))*0.5
!
!     EMPIRICAL Tue Jun 10 21:39:28 CEST 2014
!
     H(jj,ii)=X(jj,ii)*bb*hav*(r/(r-0.35))
     H(ii,jj)=0           
     X(ii,jj)=0           
   end do
   H(ii,ii)=hdiag(ii)
end do

if (abs(gscal).gt.1.d-6) then
   do i=1,nbf
      ii=aoat(i)
       do j=1,i  
         jj=aoat(j)
         eh1=0.0d0
         do kk=1,n
            eh1=eh1+q(kk)*(gab(n,xyz,at,kk,ii)+gab(n,xyz,at,kk,jj))
         end do
         H(j,i)=H(j,i)-X(j,i)*eh1*0.5d0*gscal
      end do
   end do
end if
!
!     diagonalize the hessian matrix
!
flag='N'
call SSYGVD(1,flag,'U',nbf,H,nbf,X,nbf,emo,aux, &
   &   LWORK,IWORK,LIWORK,INFO)

temp=1500.0d0

focc=0
do i=1,nel/2
   focc(i)=2.0d0
end do
if (2*(nel/2).ne.nel) then
   ihomo=nel/2+1
   ep=0
   ep(1:ihomo)=1.0
   call fermismear(.false.,nbf,ihomo,temp,emo,ep)
   focc=0
   if (ihomo.gt.1) then
      focc(1:ihomo-1)=1.0
      call fermismear(.false.,nbf,ihomo-1,temp,emo,focc)
   end if
   focc(1:nbf)=focc(1:nbf)+ep(1:nbf)
else
   ihomo=nel/2
   call fermismear(.false.,nbf,ihomo,temp,emo,focc)
   focc=focc*2.
end if

emo=emo/27.212
eel=0
do i=1,nbf 
   eel=eel+emo(i)*focc(i)
end do

deallocate(emo,ep,focc,H,X,P,S,aux,iwork)

return
end subroutine eht
