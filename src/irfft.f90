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
!     subroutine irfft: compute the inverse real fast fourier transform 
!      for a given array of data
!      Rewritten 09.10.2023 (updated fo Fortran93/C interface)
!    
!     part of EVB
!
subroutine irfft(x,N)
use fftw_mod
implicit none
integer, intent(in) :: N
double precision, intent(inout) :: x(N)
complex(C_DOUBLE_COMPLEX),dimension(N)::ain,aout

ain=x(1:N)

!
!     If this is the first execution, generate the FFT plan (optimized code for 
!     local machine)
!     If the plan already was generated for a different number of beads, destroy
!     it in advance
!
if (N .ne. Np) then
    if (Np .ne. 0) call fftw_destroy_plan(plan)
    plan=fftw_plan_dft_1d(N,ain,aout,FFTW_FORWARD,FFTW_ESTIMATE)
    factor = dsqrt(1.d0/N)
    Np = N
end if

call fftw_execute_dft(plan,ain,aout)

x = factor * real(aout(1:N))

return
end subroutine irfft

