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
!     subroutine maxwell: randomly selected speed from 
!              3D-Maxwell-Boltzmann distribution
!
!       literature reference:
!     P. W. Atkins, "Physical Chemistry, 4th Edition", W. H. Freeman,
!     New York, 1990; see section 24.2 for general discussion
!
!     part of EVB
!
function maxwell(mass1,temper)
use general
implicit none
real(kind=8)::maxwell
real(kind=8)::mass1,temper
real(kind=8)::rho,beta
real(kind=8)::random,erfinv
real(kind=8)::xspeed,yspeed
real(kind=8)::zspeed
external::erfinv
!
!
!     set normalization factor for cumulative velocity distribution
!
beta = sqrt(mass1 / (2.0d0*boltzmann*temper))
!
!     pick a randomly distributed velocity along each of three axes
!
rho = random ()
xspeed = erfinv(rho) / beta
rho = random ()
yspeed = erfinv(rho) / beta
rho = random ()
zspeed = erfinv(rho) / beta
!
!     set the final value of the particle speed in 3-dimensions
!
maxwell = sqrt(xspeed**2 + yspeed**2 + zspeed**2)

return
end function maxwell
