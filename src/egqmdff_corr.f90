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
!     subroutine egqmdff_corr: for an 2-dimensional evb-qmdff-forcefield, 
!     the energies and gradients are calculated for the 2 QMDFFs
!     HERE: for the RP-EVB coupling term, the QMDFF energies and gradients 
!     are corrected to imply a bias energy/force at the TS region to correct 
!     possible errors offside the path
!
!     part of EVB
!
subroutine egqmdff_corr (xyz2,s_act,z_act,deltq,internal,e_qmdff1,e_qmdff2,g1,g2)
use evb_mod
use general
implicit none
integer::i,j,n2
real(kind=8)::s_act,z_act  ! RP-EVB path variables
real(kind=8)::xyz2(3,natoms),e_evb,g_evb(3,natoms)
real(kind=8)::g1(3,natoms),g2(3,natoms)
real(kind=8)::deltq(nat6)   ! shift vector to the IRC
real(kind=8)::internal(nat6)  ! actual coordinates in internal coordinates
real(kind=8)::grad_add(3,natoms)  ! shift vector in cartesian coordinates
real(kind=8)::e1_shifted,e2_shifted,e1,gnorm_two
real(kind=8)::ediff,offdiag,off4,root2,deldiscr,delsqrt
real(kind=8)::e,gnorm,e_qmdff1,e_qmdff2,e2
real(kind=8)::corr_stren  ! force constant for QMDFF correction
n=natoms
!
!
!     First QMDFF      
!
call ff_eg(n,at,xyz2,e1,g1)
call ff_nonb(n,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e1,g1)
call ff_hb(n,at,xyz2,e1,g1)
!
!     Second QMDFF
!
call ff_eg_two(n,at,xyz2,e2,g2)
call ff_nonb_two(n,at,xyz2,q_two,r0ab,zab,r094_mod,sr42, &
  &             c6xy_two,e2,g2)
call ff_hb_two(n,at,xyz2,e2,g2)

!
!     Correct the energies and gradients
!
!     First, calculate the actual value of the correction force constant
!write(*,*) "s,z",s_act,z_act
!if ((s_act .lt. s_ts-deltp) .or. (s_act .gt. s_ts+deltp)) then
!   corr_stren=0.d0
!else 
!   corr_stren=corr_max-((s_act-s_ts)*sqrt(corr_max)/deltp)*((s_act-s_ts)*sqrt(corr_max)/deltp)
!end if
   if ((s_act .lt. s_ts-deltp) .or. (s_act .gt. s_ts+deltp)) then
      corr_stren=0.d0
   else
      corr_stren=corr_max*(1.d0-(1.d0/deltp)*abs(s_act-s_ts))
   end if
!corr_stren=0.d0
!
!    then apply the correction to the energies of both QMDFFs
!
!write(*,*) corr_stren*z_act*z_act
e1=e1+corr_stren*z_act*z_act
e2=e2+corr_stren*z_act*z_act
!
!    also apply the correction to the gradients of both QMDFFs
!
!    first, convert the deltq vector to cartesian coordinates
!
call int2grad(xyz2,internal,deltq,grad_add)
g1=g1+2.d0*corr_stren*grad_add!/bohr/bohr
g2=g2+2.d0*corr_stren*grad_add!/bohr/bohr
!
!     Shift the energies
!
e1_shifted=e1+E_zero1  !E(qmdff1)
e2_shifted=e2+E_zero2  !E(qmdff2)
e_qmdff1=e1_shifted
e_qmdff2=e2_shifted


!write(*,*) s_act,corr_stren
return
end subroutine egqmdff_corr

