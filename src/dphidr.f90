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
!     subroutine dphidr: calculate derivative of torsional angle!
!
!     part of QMDFF
!
subroutine dphidr(natoms,xyz,i,j,k,l,phi,dphidri, &
   &  dphidrj,dphidrk,dphidrl)
use qmdff
use pbc_mod
implicit none
external::vecnorm_loc
integer::ic,i,j,k,l,natoms

real(kind=8)::sinphi,cosphi,onenner,thab,thbc
real(kind=8)::ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3)
real(kind=8)::raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3)
real(kind=8)::rapb(3),rbpc(3),na(3),nb(3),nan,nbn
real(kind=8)::dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3)
real(kind=8)::xyz(3,natoms),phi,vecnorm_loc,nenner,eps
parameter (eps=1.d-14)

cosphi=cos(phi)
sinphi=sin(phi)
!
!     at first calculate bondlengths
!
do ic=1,3
   ra(ic)=xyz(ic,j)-xyz(ic,i)
   rb(ic)=xyz(ic,k)-xyz(ic,j)
   rc(ic)=xyz(ic,l)-xyz(ic,k)
end do
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(ra)
   call box_image(rb)
   call box_image(rc)
end if

do ic=1,3
   rapb(ic)=ra(ic)+rb(ic)
   rbpc(ic)=rb(ic)+rc(ic)
end do


call crossprod(ra,rb,na)
call crossprod(rb,rc,nb)
nan=vecnorm_loc(na,3,0)
nbn=vecnorm_loc(nb,3,0)

nenner=nan*nbn*sinphi

if (abs(nenner).lt.eps) then
   dphidri=0
   dphidrj=0
   dphidrk=0
   dphidrl=0
   return
else
   onenner=1.d0/nenner
endif

call crossprod(na,rb,rab)
call crossprod(nb,ra,rba)
call crossprod(na,rc,rac)
call crossprod(nb,rb,rbb)
call crossprod(nb,rc,rbc)
call crossprod(na,ra,raa)

call crossprod(rapb,na,rapba)
call crossprod(rapb,nb,rapbb)
call crossprod(rbpc,na,rbpca)
call crossprod(rbpc,nb,rbpcb)
!
!     derivatives of the angle to all four possible distances:
!     dphidri
!
do ic=1,3
   dphidri(ic)=onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))
!
!     dphidrj
!
   dphidrj(ic)=onenner*(cosphi*(nbn/nan*rapba(ic) &
          &   +nan/nbn*rbc(ic))-(rac(ic)+rapbb(ic)))
!
!     dphidrk
!
  dphidrk(ic)=onenner*(cosphi*(nbn/nan*raa(ic) &
          &   +nan/nbn*rbpcb(ic))-(rba(ic)+rbpca(ic)))
!
!     dphidrl
!
   dphidrl(ic)=onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
end do

return
end subroutine dphidr
