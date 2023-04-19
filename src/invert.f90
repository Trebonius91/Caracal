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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine invert: inverts a matrix using the Gauss-Jordan method
!     needed for dynamics
!     variables and parameters:
!     n     dimension of the matrix to be inverted
!     a     matrix to invert; contains inverse on exit
!
!     part of EVB
!
subroutine invert (n,a)
implicit none
integer::i,j,k,n
integer::icol,irow
integer,allocatable::ipivot(:)
integer,allocatable::indxc(:)
integer,allocatable::indxr(:)
real(kind=8)::big,temp
real(kind=8)::pivot
real(kind=8)::a(n,*)
!
!     perform dynamic allocation of some local arrays
!
allocate (ipivot(n))
allocate (indxc(n))
allocate (indxr(n))
!
!     perform matrix inversion via the Gauss-Jordan algorithm
!
do i = 1, n
   ipivot(i) = 0
end do
do i = 1, n
   big = 0.0d0
   do j = 1, n
      if (ipivot(j) .ne. 1) then
         do k = 1, n
            if (ipivot(k) .eq. 0) then
               if (abs(a(j,k)) .ge. big) then
                  big = abs(a(j,k))
                  irow = j
                  icol = k
               end if
            else if (ipivot(k) .gt. 1) then
               write (*,'(/," INVERT  --  Cannot Invert a Singular Matrix")')
               call fatal
            end if
         end do
      end if
   end do
   ipivot(icol) = ipivot(icol) + 1
   if (irow .ne. icol) then
      do j = 1, n
         temp = a(irow,j)
         a(irow,j) = a(icol,j)
         a(icol,j) = temp
      end do
   end if
   indxr(i) = irow
   indxc(i) = icol
   if (a(icol,icol) .eq. 0.0d0) then
      write (*,'(/," INVERT  --  Cannot Invert a Singular Matrix")')
      call fatal
   end if
   pivot = a(icol,icol)
   a(icol,icol) = 1.0d0
   do j = 1, n
      a(icol,j) = a(icol,j) / pivot
   end do
   do j = 1, n
      if (j .ne. icol) then
         temp = a(j,icol)
         a(j,icol) = 0.0d0
         do k = 1, n
            a(j,k) = a(j,k) - a(icol,k)*temp
         end do
      end if
   end do
end do
do i = n, 1, -1
   if (indxr(i) .ne. indxc(i)) then
      do k = 1, n
         temp = a(k,indxr(i))
         a(k,indxr(i)) = a(k,indxc(i))
         a(k,indxc(i)) = temp
      end do
   end if
end do
!
!     perform deallocation of some local arrays
!
deallocate (ipivot)
deallocate (indxc)
deallocate (indxr)

return
end subroutine invert
