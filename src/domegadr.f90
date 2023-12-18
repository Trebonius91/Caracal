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
!     calculate derivative of inversion angle
!     similar as dphidr..
!
!     part of QMDFF
! 
subroutine domegadr(natoms,xyz,i,j,k,l,omega, &
    &   domegadri,domegadrj,domegadrk,domegadrl)
use qmdff
use pbc_mod
implicit none
external::vecnorm_loc
integer::ic,i,j,k,l,natoms
real(kind=8)::omega,sinomega
real(kind=8)::xyz(3,natoms),onenner,vecnorm_loc,rnn,rvn
real(kind=8)::rn(3),rv(3),rd(3),re(3),rdme(3),rve(3)
real(kind=8)::rne(3),rdv(3),rdn(3)
real(kind=8)::rvdme(3),rndme(3),nenner
real(kind=8)::domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3),eps
parameter(eps=1.d-14)

sinomega=sin(omega)

do ic=1,3
   rv(ic)=xyz(ic,l)-xyz(ic,i)
   rd(ic)=xyz(ic,k)-xyz(ic,j)
   re(ic)=xyz(ic,i)-xyz(ic,j)
end do

!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(rv)
   call box_image(rd)
   call box_image(re)
end if

do ic=1,3
   rdme(ic)=rd(ic)-re(ic)
end do

call crossprod(re,rd,rn)
rvn=vecnorm_loc(rv,3,0)
rnn=vecnorm_loc(rn,3,0)

call crossprod(rv,re,rve)
call crossprod(rn,re,rne)
call crossprod(rd,rv,rdv)
call crossprod(rd,rn,rdn)
call crossprod(rv,rdme,rvdme)
call crossprod(rn,rdme,rndme)

nenner=rnn*rvn*cos(omega)
if (abs(nenner).gt.eps) then
   onenner=1.d0/nenner
   do ic=1,3
!
!     domega/dri
!
      domegadri(ic)=onenner*(rdv(ic)-rn(ic)- &
            &  sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))
!
!     domega/drj
!
      domegadrj(ic)=onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))
!
!     domega/drk
!
      domegadrk(ic)=onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))
!
!     domega/drl
!
      domegadrl(ic)=onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
   end do
else
   do ic=1,3
      domegadri(ic)=0.d0
      domegadrj(ic)=0.d0
      domegadrk(ic)=0.d0
      domegadrl(ic)=0.d0
   end do
end if

return
end subroutine domegadr
