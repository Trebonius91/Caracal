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
!     subroutine temper: applies a velocity correction at the half time step
!     as needed for the Nose-Hoover extended system thermostat
!
!     literature references:
!
!     D. Frenkel and B. Smit, "Understanding Molecular Simulation,
!     2nd Edition", Academic Press, San Diego, CA, 2002; see Appendix
!     E.2 for implementation details
!
!     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
!     "Explicit Reversible Integrators for Extended Systems Dynamics",
!     Molecular Physics, 87, 1117-1157 (1996)
!
!     part of EVB
!
subroutine temper (dt)
use general
implicit none
integer::i,j,nc,ns
real(kind=8)::dt,dtc,dts
real(kind=8)::dt2,dt4,dt8
real(kind=8)::eksum,ekt
real(kind=8)::scale,expterm
real(kind=8)::w(3)
real(kind=8)::ekin(3,3)
!
!     make half-step velocity correction for Nose-Hoover system
!

return
end subroutine temper
