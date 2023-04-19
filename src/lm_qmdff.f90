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
!     subroutine lm_qmdff: calculate Levenberg-Marquardt fitness
!     for a QMDFF parameter optimization calculation
!
!     Part of EVB
!

subroutine lm_qmdff(m_ind,n_ind,x_var,fvec,iflag)
use evb_mod
use lm_module
use qmdff
implicit none
integer::m_ind,n_ind
real(kind=8)::fvec(m_ind),fun
real(kind=8)::dg_func,E1,E2,V12,deltaE,root
integer::i,iflag,j
real(kind=8)::x_var(n_ind),xdat(m_ind),ydat(m_ind)
real(kind=8)::e_one 
real(kind=8)::xyz_i(3,n_one)
real(kind=8)::e_vec(m_ind)
!
!     the alpha values are now stored into the x(:) array 
!     n=number of alphas, m=number of path points
!
mn_par(1)=x_var(1)
mn_par(2:n_ind)=x_var(2:n_ind)
mn_par(n_ind+1:2*n_ind)=x_var(2:n_ind)
do i=1,m_ind
   e_one=0.d0
   xyz_i(1:3,1:n_one)=xyz_ref_lm(1:3,1:n_one,i)
   call ff_eg(n_one,at,xyz_i,e_one,g_one)
   call ff_nonb(n_one,at,xyz_ref_lm(:,:,i),q,r0ab,zab,r094_mod,sr42,c6xy,e_one,g_one)
   call ff_hb(n_one,at,xyz_ref_lm(:,:,i),e_one,g_one)
   e_vec(i)=e_one
!   fvec(i)=(e_one-energies_ref_lm(i))*(e_one-energies_ref_lm(i))
end do
!
!    Shift the QMDFF energies such that the minimum values coincide
!
e_vec=e_vec-(minval(e_vec)-minval(energies_ref_lm))
do i=1,m_ind
   fvec(i)=(e_vec(i)-energies_ref_lm(i))*(e_vec(i)-energies_ref_lm(i))
end do
!
!    Add parabolic corrections if parameters seem to diverge to huge/small values
!
!do i=1,n_ind
!   if (x_var(i) .lt. 1.d0-2*lower_bond) then
!      f_vec(1)=f_vec(1)+(x_var(i)-(1.d0-2*lower_bond))**2*0.1d0
!   else if (x_var(i) .gt. 1.d0+2*lower_bond) then
!      f_vec(1)=f_vec(1)+(x_var(i)-(1.d0+2*lower_bond))**2*0.1d0
!   end if
!end do

return
end subroutine lm_qmdff


