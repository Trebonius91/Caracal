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
!     subroutine lm_function: calculate Levenberg-Marquardt fitness
!     for a DG-EVB calculation
!

subroutine lm_function(m_ind,n_ind,x_var,fvec,iflag)
use evb_mod
use lm_module
implicit none
integer::m_ind,n_ind
integer::slope_sign,sign_old ! for wiggle scoring
real(kind=8)::fvec(m_ind),fun
real(kind=8)::dg_old  ! for wiggle scoring
real(kind=8)::dg_func,E1,E2,V12,deltaE,root
integer::i,iflag
real(kind=8)::x_var(n_ind),xdat(m_ind),ydat(m_ind)
!
!     the alpha values are now stored into the x(:) array 
!     n=number of alphas, m=number of path points
!
call build_dmat(dg_evb_mode_lm,mat_size_lm,x_var)
call mat_diag(mat_size_lm)
extrema=0  ! number of EVB energy maxima/minima along the path
do i=1,m_ind
   V12=0
   E1=ff_e1_lm(i)
   E2=ff_e2_lm(i)
   deltaE=abs(E1-E2)
   call sum_v12(dg_evb_mode_lm,mat_size_lm,x_var,path_int_lm(:,i),V12)
   root=(0.5d0*(E1-E2))*(0.5d0*(E1-E2))+V12
   if (root .le. 0) then
      dg_func=0.5*(E1+E2)
   else
      dg_func=0.5*(E1+E2)-sqrt(root)
   end if
!
!     To avoid many artifical wiggles, score them and apply penalty
!
!   if (i .gt. 1) then
!      if (dg_func-dg_old .lt. 0d0) then
!     positive slope
!         slope_sign=+1 
!      else 
!     negative slope
!         slope_sign=-1
!      end if 
!
!     increment the number of extrema if the slope has changed 
!
!      if (slope_sign .ne. sign_old) then
!         extrema=extrema+1
!      end if
!   end if
!   sign_old=slope_sign
!   dg_old=dg_func
   fvec(i)=(dg_func-energies_ref_lm(i))*(dg_func-energies_ref_lm(i))
!
!     if more than three extrema (maxima and TS) are present, apply 
!     the penalty; the higher the more extrema are present!
!
!   if (extrema .gt. 3) then
!       write(*,*) extrema
!      f_vec=f_vec+0.001d0*(extrema-3)
!   end if
end do

return
end subroutine lm_function


