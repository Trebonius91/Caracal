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
!    subroutine getc6: Interpolate C6 parameters according to standard D3
!
!    part of QMDFF
!

subroutine getc6(maxc,max_elem,c6ab,mxc,iat,jat,nci,ncj,c6)
implicit none
integer::maxc,max_elem
integer::iat,jat,i,j,mxc(max_elem)
real(kind=8)::nci,ncj,c6,c6mem
real(kind=8)::c6ab(max_elem,max_elem,maxc,maxc,3)
!
! the exponential is sensitive to numerics
! when nci or ncj is much larger than cn1/cn2
!
real(kind=8)::cn1,cn2,r,rsum,csum,tmp,tmp1
real(kind=8)::r_save
real(kind=8)::k1,k3
parameter (k1     =16)
parameter (k3     =-4)

c6mem=-1.d+99
rsum=0.0
csum=0.0
c6  =0.0
r_save=10000.0   ! changed from 1000.0 for test reasons!!! (12.05.2017)
do i=1,mxc(iat)
   do j=1,mxc(jat)
      c6=c6ab(iat,jat,i,j,1)
      if (c6.gt.0)then
         cn1=c6ab(iat,jat,i,j,2)
         cn2=c6ab(iat,jat,i,j,3)
!
!     distance measure
!
         r=(cn1-nci)**2+(cn2-ncj)**2
         if (r.lt.r_save) then
            r_save=r
            c6mem=c6
         end if
         tmp1=exp(k3*r)
         rsum=rsum+tmp1
         csum=csum+tmp1*c6
      end if
   end do
end do

if (rsum.gt.1.0d-99) then
   c6=csum/rsum
else
   c6=c6mem
end if

end subroutine getc6
