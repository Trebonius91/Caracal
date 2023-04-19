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
!     subroutine dsyprj: projection of a real symmetric matrix asym 
!        packed in an upper triangular form:
!        asym = (1 - bmat*bmat')*asym*(1 - bmat*bmat')
!        where (1 - bmat*bmat') is a projector constructed from bmat
!
!     part of QMDFF
!
subroutine dsyprj(nbdim,m,bmat,n,asym)
implicit none

integer,intent(in)::nbdim,m,n
real(kind=8),dimension(nbdim,m)::bmat
real(kind=8),dimension(n*(n+1)/2)::asym
integer::i,j,ij
real(kind=8),dimension(n,m)::scrb
real(kind=8),dimension(n,n)::scra
!
!     Expand trigonal matrix asym to full matrix on scra
!
call htosq(n,scra,asym)
!
!     Calculate scrb = asym*bmat
!
call dsymm('l','u',n,m,1.0d0,scra,n,bmat,nbdim,0.0d0,scrb,n)
!
!     Calculate scra = scrb*bmat'
!
call dgemm('n','t',n,n,m,1.0d0,scrb,n,bmat,nbdim,0.0d0,scra,n)
!
!     Calculate asym = asym - scra - scra'
!
do i=1,n
   do j=1,i
      ij = i*(i-1)/2 + j
      asym(ij) = asym(ij) - scra(i,j) - scra(j,i)
   end do
end do
!
!     Calculate scrb = scra'*bmat
!
call dgemm('t','n',n,m,n,1.0d0,scra,n,bmat,nbdim,0.0d0,scrb,n)
!
!     Calculate scra = bmat*scrb'
!
call dgemm('n','t',n,n,m,1.0d0,bmat,nbdim,scrb,n,0.0d0,scra,n)
!
!     Calculate asym = asym + scra
!
do i=1,n
   do j=1,i
      ij = i*(i-1)/2 + j
      asym(ij) = asym(ij) + scra(i,j)
   end do
end do

return
end subroutine dsyprj
