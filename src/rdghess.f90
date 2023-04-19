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
!     subroutine rdghess: Read the gaussian reference hessian, 
!       analogous to gaurd section for QMDFF
!
!     part of QMDFF
!

subroutine rdghess(n3,h,fname,prefix1,finished)
implicit none
integer::n3
integer::finished ! if formchk should be executed
real(kind=8)::h(n3,n3)
real(kind=8)::hess1d(((n3)*(n3+1))/2)
real(kind=8)::hess1d2(n3*n3)
character(len=*)::fname
character(len=*)::prefix1
integer::k1,k2,j1,ir,n,nkpb,ibl,j2,kk,kd,k1s,j,i,k,iuout,jdum,nn
character(len=80)::adum,a80
real(kind=8)::xx(10)
real(kind=8)::r(n3*(n3+1)/2)
logical::hread
integer::fchkstat ! for gaussian output: status of fchk output conversion!
integer::nat3,maxcol,mincol,natsum,sum1,sum2

n=n3
nat3=n3
nkpb=5
iuout=1
r=0

open(unit=42,file=prefix1 // ".fchk")
12  read(42,'(a)',end=32)a80
!
!   Read in the hessian (upper triangular matrix)
!
   if(index(a80,'Cartesian Force Constants').ne.0)then
      do i=1,int(((nat3*(nat3+1))/10)+1)
         mincol = maxcol + 1
         maxcol = min(maxcol+5,(nat3*(nat3+1))/2)
         read(42,*) (hess1d(j),j=mincol,maxcol)

      enddo
      do i=1,nat3
         do j=1,nat3
            if (j .le. i) then
               natsum=((i-1)*i)/2
               sum1=(i-1)*nat3+j
               sum2=(j-1)*nat3+i
               hess1d2(sum1)=hess1d(natsum+j)
               hess1d2(sum2)=hess1d2(sum1)
            end if
         end do
      end do
   end if

   goto 12
   stop
32    close(11)

!   write(*,*) grad1d

   do i=1,nat3
      do j=1,nat3
         h(i,j)=hess1d2((i-1)*nat3+j)
      end do
   end do
close(42)

65 continue
open(iuout,file=fname)
do
   read(iuout,'(a)',end=100) adum
   if(index(adum,'Force constants in Cartesian').ne.0)then
      hread=.true.
      ibl=n/nkpb
      ir=n-ibl*nkpb
      j1=1
      k1s=1
      kd=0
      if(ibl.eq.0) go to 50
      j2=nkpb
      do i=1,ibl
         read(iuout,'(a)') adum
         k1=k1s
         k2=k1
         kk=0
         do j=j1,j2
            read(iuout,1003)jdum,(r(k),k=k1,k2)
            kk=kk+1
            k1=k1+kd+kk
            k2=k1+kk
         end do
         j1=j1+nkpb
         if (j1.gt.n) goto 999
         j2=j2+nkpb
         k2=k1-1
         k1=k2+1
         k2=k1+(nkpb-1)
         k1s=k2+1
         kk=kd+nkpb
         do j=j1,n
            read(iuout,1003) jdum,(r(k),k=k1,k2)
            kk=kk+1
		    k1=k1+kk
            k2=k2+kk
         end do
          kd=kd+nkpb
      end do
   50 if(ir.eq.0) go to 70
      k1=k1s
      j2=j1+ir-1
      kk=0
      k2=k1
      read(iuout,'(a)') adum
      do  j=j1,j2
         read(iuout,1003) jdum,(r(k),k=k1,k2)
         kk=kk+1
         k1=k1+kd+kk
         k2=k1+kk
      end do
   70 goto 999
  999 continue
      k=0
      do i=1,n
         do j=1,i
            k=k+1
            h(j,i)=r(k)
            h(i,j)=r(k)
         end do
      end do
   end if
!
!     end of loop over file lines
!
end do
100  close(iuout)
1003 format(i7,5d14.6)
return
end subroutine rdghess
