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
!     "dftmod" computes the modulus of the discrete Fourier transform
!     of "bsarray" and stores it in "bsmod"
!
!
subroutine dftmod (bsmod,bsarray,nfft,order)
implicit none
integer::i,j,k
integer::nfft,jcut
integer::order,order2
real(kind=8)::eps,zeta
real(kind=8)::arg,factor
real(kind=8)::sum1,sum2
real(kind=8)::bsmod(*)
real(kind=8)::bsarray(*)
real(kind=8)::pi
pi=3.1415926535897932384626433832795029d0
!
!     get the modulus of the discrete Fourier transform
!
factor = 2.0d0 * pi / dble(nfft)
do i = 1, nfft
   sum1 = 0.0d0
   sum2 = 0.0d0
   do j = 1, nfft
      arg = factor * dble((i-1)*(j-1))
      sum1 = sum1 + bsarray(j)*cos(arg)
      sum2 = sum2 + bsarray(j)*sin(arg)
   end do
   bsmod(i) = sum1**2 + sum2**2
end do
!
!     fix for exponential Euler spline interpolation failure
!
eps = 1.0d-7
if (bsmod(1) .lt. eps)  bsmod(1) = 0.5d0 * bsmod(2)
do i = 2, nfft-1
   if (bsmod(i) .lt. eps)  &
        &   bsmod(i) = 0.5d0 * (bsmod(i-1)+bsmod(i+1))
end do
if (bsmod(nfft) .lt. eps)  bsmod(nfft) = 0.5d0 * bsmod(nfft-1)
!
!     compute and apply the optimal zeta coefficient
!
jcut = 50
order2 = 2 * order
do i = 1, nfft
   k = i - 1
   if (i .gt. nfft/2)  k = k - nfft
   if (k .eq. 0) then
      zeta = 1.0d0
   else
      sum1 = 1.0d0
      sum2 = 1.0d0
      factor = pi * dble(k) / dble(nfft)
      do j = 1, jcut
         arg = factor / (factor+pi*dble(j))
         sum1 = sum1 + arg**order
         sum2 = sum2 + arg**order2
      end do
      do j = 1, jcut
         arg = factor / (factor-pi*dble(j))
         sum1 = sum1 + arg**order
         sum2 = sum2 + arg**order2
      end do
      zeta = sum2 / sum1
   end if
   bsmod(i) = bsmod(i) * zeta**2
end do

return
end subroutine dftmod


