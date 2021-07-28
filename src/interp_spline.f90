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
!     The subroutine interp_spline needs a certain value of the parameter 
!     s to obtain the internal coordinates and the coupling strength
!     for that respective point on the 1D reaction path.
!
!     part of EVB
!
subroutine interp_spline(dim,s_act,path_s,klo,khi)
use evb_mod
implicit none
integer::intdim,dim   ! number of internal coordinates and total dimension
integer::n  ! number of points om the path
integer::i,j,k  ! loop indices
integer::j1,klo,khi  ! for interpolation routine
real(kind=8)::h,a,b    ! for interpolation routine
real(kind=8)::path_s(dim)   ! the resulting array with the internal coordinates 
                                 ! and coupling strength for the s_act value
real(kind=8)::s_act  ! the actual s value that shall be interpolated


!
!     Determine, in which interval between two reference points the 
!     current s-value is located. Perform the local cubic spline from 
!     this two points
!
n=rp_irc_points
do k=1,n-1
   if (s_act .ge. rp_spl_s(k) .and. s_act .le. rp_spl_s(k+1)) then
      j1=k
      exit
   end if
end do
klo=j1
khi=j1+1
h=rp_spl_s(khi)-rp_spl_s(klo)
!
!     The x values need to have distinct values!
!
if (h .eq. 0) then
   write(*,*)  "Error in interp_spline.f90: two points are equal!"
   write(*,*)  "Control your inpur for the IRC!"
   call fatal
end if
a=(rp_spl_s(khi)-s_act)/h
b=(s_act-rp_spl_s(klo))/h
!
!     Now determine the coordinate values for all of the dimensions!
!
do j=1,dim
   path_s(j)=a*rp_spl_ref(j,klo)+b*rp_spl_ref(j,khi)+((a*a*a-a)*& 
           & rp_spl_d2(klo,j)+(b*b*b-b)*rp_spl_d2(khi,j))*(h*h)/6.d0
end do
return
end subroutine interp_spline
