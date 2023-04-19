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
!     subroutine egrad_mueller: Calculating energy and gradient for the famous
!       Mueller-Brown analytical 2D surface. This can be used for test cases.
!
!     part of EVB
!
subroutine egrad_mueller(coords,Epot,grad,info)
implicit none
!     The input coordinates
real(kind=8), intent(in) :: coords(3)
!     The energy
real(kind=8), intent(out) :: Epot
!     The gradient
real(kind=8), intent(out) :: grad(3)
!     Info parameter (not used)
integer, intent(out) :: info
!     Single coordinate components
real(kind=8)::x,y
!     Mueller-Brown parameters
real(kind=8)::A1_par(4)
real(kind=8)::a2_par(4)
real(kind=8)::x0_par(4)
real(kind=8)::y0_par(4)
real(kind=8)::b_par(4)
real(kind=8)::c_par(4)
!     Auxiliary array: sum components
real(kind=8)::e_parts(4)
!     Loop indices
integer::i
!
!      Fill the parameter vectors
!
A1_par=[-200.d0,-100.d0,-170.d0,15.d0]
A2_par=[-1.d0,-1.d0,-6.5d0,0.7d0]
x0_par=[1.d0,0.d0,-0.5d0,-1.d0]
y0_par=[0.d0,0.5d0,1.5d0,1.d0]
b_par=[0.d0,0.d0,11.d0,0.6d0]
c_par=[-10.d0,-10.d0,-6.5d0,0.7d0]
!
!      determine current coordinates
!
x=coords(1)
y=coords(2)
!
!      calculate the energy for the actual point
!
Epot=0.d0
do i=1,4
   e_parts(i)=A1_par(i)*exp(a2_par(i)*(x-x0_par(i))*(x-x0_par(i))+&
        &  b_par(i)*(x-x0_par(i))*(y-y0_par(i))+c_par(i)*&
        &  (y-y0_par(i))*(y-y0_par(i)))
end do
Epot=sum(e_parts)
!
!       calculate the gradient for the actual point
!
grad=0.d0
do i=1,4
   grad(1)=grad(1)+A1_par(i)*(2.d0*a2_par(i)*(x-x0_par(i))+b_par(i)*&
        & (y-y0_par(i)))*e_parts(i)
   grad(2)=grad(2)+A1_par(i)*(2.d0*c_par(i)*(y-y0_par(i))+b_par(i)*&
        & (x-x0_par(i)))*e_parts(i)
end do


return
end subroutine egrad_mueller
