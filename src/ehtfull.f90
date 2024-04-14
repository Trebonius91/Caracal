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
!     subroutine ehtfull: If Wiberg-Mayer bond orders for QMDFF parameterization 
!        don´t exist, calculate them based on EHT!
!
!     part of QMDFF
!
subroutine ehtfull(n,at,xyz,q,nel,nbf,eel,wb,gscal)
use qmdff
implicit none
integer::n,at(n),nel,nbf
real(kind=8)::xyz(3,n),eel,wb(n,n),q(n),gscal
logical::dens

real(kind=8)::b,bb,temp,xsum,eh1,gab,rab,brow(7)
real(kind=4),allocatable::X(:,:)
real(kind=4),allocatable::P(:,:)
real(kind=4),allocatable::S(:,:)
real(kind=4),allocatable::H(:,:)
real(kind=4),allocatable::aux(:)
real(kind=4),allocatable ::emo(:)
real(kind=4),allocatable ::ep (:)
real(kind=4),allocatable ::focc(:)
integer,allocatable::iwork(:)
integer,allocatable::iiwork(:)
integer::LWORK,LIWORK,INFO
integer::npr,ii,jj,i,j,kk,m,k,iat,jat,lina
integer::mmm(10),xi,xj,ihomo,irow,jrow,pqn
character(len=1)::flag

mmm=(/1,2,2,2,3,3,3,3,3,3/)
liwork = 3 + 5*nbf
lwork  = 1 + 6*nbf + 2*nbf**2
allocate(H(nbf,nbf),X(nbf,nbf),aux(lwork),P(nbf,nbf), &
  & iwork(liwork),emo(nbf),focc(nbf),ep(nbf),S(nbf,nbf))
!
!     print out info message
!
write(10,*)'Doing full Extended Hückel calculation ...'
write(10,*)'Number of atoms    : ',n
write(10,*)'Number of electrons  : ',nel
write(10,*)'Number of basis functions  : ',nbf

!
!     calculate overlap integrals
!
call stints(n,nbf,xyz,X)
S = X

b=1.75
do irow=1,6
   brow(irow)=b
   b=b+0.20
end do

do ii=1,nbf
   xi=aoat(ii)
   irow=pqn(at(xi))
   do jj=1,ii-1
     xj=aoat(jj)         
     jrow=pqn(at(xj))
     bb=0.5*(brow(irow)+brow(jrow))
     H(jj,ii)=X(jj,ii)*bb*(hdiag(ii)+hdiag(jj))*0.5
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
!     diagonalize the Hamiltonian matrix
!     by calling the lapack routine 
!     (generalized symmetric definite eigenproblem)
!
flag='V'
call SSYGVD(1,flag,'U',nbf,H,nbf,X,nbf,emo,aux, &
   & LWORK,IWORK,LIWORK,INFO)
temp=1500.0d0
focc=0
do i=1,nel/2
   focc(i)=2.0d0
end do
if (2*(nel/2).ne.nel) then
   ihomo=nel/2+1
   ep=0
   ep(1:ihomo)=1.0
   call fermismear_eht(.false.,nbf,ihomo,temp,emo,ep)
   focc=0
   if (ihomo.gt.1) then
      focc(1:ihomo-1)=1.0
      call fermismear_eht(.false.,nbf,ihomo-1,temp,emo,focc)
   end if
   focc(1:nbf)=focc(1:nbf)+ep(1:nbf)
else
   ihomo=nel/2
   call fermismear_eht(.false.,nbf,ihomo,temp,emo,focc)
   focc=focc*2.
end if
emo=emo/27.212
eel=0
do i=1,nbf 
   eel=eel+emo(i)*focc(i)
end do

!     calc P and WBOs     
write(10,*) 
write(10,*)'Total electronic energy  : ',eel 
!
!     looking at the name it seems that it print out 
!     eigenvalues and/or eigenvectors
!
write(10,*)
write(10,*) "Print out of all EHT orbitals (eigenvalues and occupations)"
write(10,*)
call preig4(10,focc,emo,nbf)
do m=1,nbf  
   do i=1,nbf
      X(i,m)=H(i,m)
      H(i,m)=H(i,m)*focc(m)
   end do
end do
!
!     call LAPACK routine for operations with three matrices..
!
call sgemm('N','T',nbf,nbf,nbf,1.0,X, &
   &  nbf,H,nbf,0.0,P,nbf)
call sgemm('N','N',nbf,nbf,nbf,1.0,P, &
   &  nbf,S,nbf,0.0,H,nbf)
wb=0
do i=1,n 
   do j=1,i-1
   xsum=0.0
   rab=(xyz(1,i)-xyz(1,j))**2 &
    & +(xyz(2,i)-xyz(2,j))**2 &
    & +(xyz(3,i)-xyz(3,j))**2
   if(rab.lt.100.0)then
      do k=fila(1,i),fila(2,i)
         do m=fila(1,j),fila(2,j)
            xsum=xsum+H(k,m)*H(m,k)
         end do
      end do
   end if
!
!     this scales WBO to about 1,2,3 for ethane,ethene,ethyne            
!     the model WBO are NOT used, this is just for convenience
!     and in strange cases when the DFT (orca.out) WBO do not
!     exist
!
      wb(i,j)=xsum*1.0
      wb(j,i)=xsum*1.0
   end do
end do

deallocate(emo,ep,focc,H,X,P,S,aux,iwork)
return
end subroutine ehtfull
