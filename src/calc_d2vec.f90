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
!     subroutine calc_d2vec: Calculate the second derivative tensor of 
!     a normalized three dimensional difference bond vector as it is 
!     needed for calculation of Wilson matrix derivatives of internal
!     coordinates 
!     Input: the normalized vector, the empty matrix and the vector length
!
!     part of EVB
!

subroutine calc_d2vec(d2_mat,v,len)
implicit none 
!     The normalized vector to be derivated 
real(kind=8), intent(in) :: v(3)
!     The length of the given vector 
real(kind=8), intent(in) :: len
!     The calculated second derivative matrix 
real(kind=8), intent(out) :: d2_mat(3,3,3)

!
!     Calculate entries of second derivative matrix 
!     --> the matrix is diagonal!
!
!     dv/dx^2
d2_mat(1,1,1)=-3d0*v(1)*(v(2)**2+v(3)**2)
d2_mat(1,1,2)=v(2)*(2d0*v(1)**2-v(3)**2-v(2)**2)
d2_mat(1,1,3)=v(3)*(2d0*v(1)**2-v(2)**2-v(3)**2)
!     dv/dxdy
d2_mat(1,2,1)=-v(2)*(-2d0*v(1)**2+v(2)**2+v(3)**2)
d2_mat(1,2,2)=-v(1)*(v(1)**2-2d0*v(2)**2+v(3)**2)
d2_mat(1,2,3)=3d0*v(1)*v(2)*v(3)
!     dv/dxdz
d2_mat(1,3,1)=-v(3)*(-2d0*v(1)**2+v(2)**2+v(3)**2)
d2_mat(1,3,2)=3d0*v(1)*v(2)*v(3)
d2_mat(1,3,3)=-v(1)*(v(1)**2+v(2)**2-2d0*v(3)**2)
!     dv/dydx
d2_mat(2,1,:)=d2_mat(1,2,:)
!     dv/dydy
d2_mat(2,2,1)=v(1)*(-v(1)**2+2d0*v(2)**2-v(3)**2)
d2_mat(2,2,2)=-3d0*v(2)*(v(1)**2+v(3)**2)
d2_mat(2,2,3)=v(3)*(-v(1)**2+2d0*v(2)**2-v(3)**2)
!     dv/dydz
d2_mat(2,3,1)=3d0*v(1)*v(2)*v(3)
d2_mat(2,3,2)=-v(3)*(v(1)**2-2d0*v(2)**2+v(3)**2)
d2_mat(2,3,3)=-v(2)*(v(1)**2+v(2)**2-2d0*v(3)**2)
!     dv/dzdx
d2_mat(3,1,:)=d2_mat(1,3,:)
!     dv/dzdy
d2_mat(3,2,:)=d2_mat(2,3,:)
!     dv/dzdz
d2_mat(3,3,1)=v(1)*(-v(1)**2-v(2)**2+2d0*v(3)**2)
d2_mat(3,3,2)=v(2)*(-v(1)**2-v(2)**2+2d0*v(3)**2)
d2_mat(3,3,3)=-3d0*v(3)*(v(1)**2+v(2)**2)
!
!     apply length prefactor 
!
d2_mat=d2_mat/(len**5)

return
end subroutine calc_d2vec
