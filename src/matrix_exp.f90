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
!     Subroutine matrix_exp: Return the expoential of a (square) matrix 
!       using the scale and square algorithm.
!    Parameters:
!   M - The square matrix to exponentiate
!   n - The size of the matrix
!   j - The number of terms in the exponential to evaluate
!   k - The exponent of the power of two used to scale the matrix
!    Returns:
!   EM - The expoential of the given matrix
!    --> taken from RPMDrate
!
!    part of EVB
!
subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real*8, intent(in)   :: M(n,n)
    real*8, intent(out)   :: EM(n,n)

    real *8 :: tc(j+1), SM(n,n)
    integer p, i
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo

    !scale
    SM=M*(1./2.**k)
    EM=0.
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo

    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo

    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo

end subroutine matrix_exp

