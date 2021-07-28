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
!     subroutine cholesky: Compute the brute-force stabilized 
!        Cholesky decomposition of a square matrix.
!     Parameters:
!   SST - The matrix to determine the Cholesky decomposition of
!   n - The size of the matrix
!     Returns:
!   S - The computed Cholesky decomposition
!     ---> taken from RPMDrate 
!
!     part of EVB
!
subroutine cholesky(SST, S, n)

    integer, intent(in)  :: n
    real*8, intent(in)   :: SST(n,n)
    real*8, intent(out)   :: S(n,n)
    real*8 :: D(n), L(n,n)

    integer i,j,k

    D = 0.d0
    L = 0.d0
    do i = 1,n
        L(i,i)=1.
        D(i)=SST(i,i)
        do j=1,i-1
            L(i,j)=SST(i,j);
            do k=1,j-1
                L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
            end do
            if (D(j).ne.0.) L(i,j)=L(i,j)/D(j)
        end do
        do k=1,i-1
            D(i)=D(i)-L(i,k)*L(i,k)*D(k)
        end do
    end do
    S=0.
    do i=1,n
        do j=1,i
            if (D(j)>0.) S(i,j)=S(i,j)+L(i,j)*sqrt(D(j))
        end do
    end do

end subroutine cholesky

