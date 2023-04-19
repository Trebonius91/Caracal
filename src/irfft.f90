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
!    
!     part of EVB
!
subroutine irfft(x,N)
implicit none

integer, intent(in) :: N
double precision, intent(inout) :: x(N)

integer, parameter :: Nmax = 1024
integer :: Np
double precision :: copy(Nmax), factor
integer*8 :: plan

data Np /0/
save copy, factor, plan, Np

if (N .ne. Np) then
!
!     The input array is a different length than the last array, so we
!     must generate a new FFTW plan for the transform
!     First delete the previous plan
!
    if (Np .ne. 0) call dfftw_destroy_plan(plan)
    call dfftw_plan_r2r_1d(plan,N,copy,copy,1,64)
    factor = dsqrt(1.d0/N)
    Np = N
end if

copy(1:N) = x
call dfftw_execute(plan)
x = factor * copy(1:N)


return
end subroutine irfft
