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
!     subroutine gaurd: Read in coordinates and atom types from gaussian 
!      output for QMDFF reference 
!
!     part of QMDFF
!
subroutine gaurd(fname,nat,at,h,xyz,q,wbo,check_coord)
implicit none
integer::nat
character(len=*) fname
real(kind=8)::h(3*nat,3*nat),q(nat),xyz(3,nat),wbo(nat,nat)
integer::at(nat)

integer::k1,k2,j1,ir,n,nkpb,ibl,j2,kk,kd,k1s,j,i,k,iuout,jdum,nn
character(len=80)::adum
real(kind=8)::xx(10)
real(kind=8)::r(3*nat*(3*nat+1)/2)
logical::hread,check_coord

hread=.false.
n=3*nat
nkpb=5
iuout=1
r=0
xyz=0
q=0
at=0
wbo=0

open(iuout,file=fname)
      
do
   read(iuout,'(a)',end=100) adum
!
!     read coordinates
!           
   if (index(adum,'Standard orientation:').ne.0) then
      read(iuout,'(a)',end=100) adum
      read(iuout,'(a)',end=100) adum
      read(iuout,'(a)',end=100) adum
      read(iuout,'(a)',end=100) adum
      do i=1,nat
         read(iuout,'(a)',end=100) adum
         call readl(adum,xx,nn)
         xyz(1:3,i)=xx(4:6)/0.52917726
         at(i)=idint(xx(2))
      end do
   end if
!
!     read Hirshfeld charges
! 
   if (index(adum,'Hirshfeld charges, spin densities, dipoles').ne.0) then
      read(1,'(a)',end=100) adum
      do i=1,nat
         read(iuout,'(a)',end=100) adum
         call readl(adum,xx,nn)
  !       write(*,*) "XX",xx
         q(i)=xx(2)
      end do
   end if
!
!     Read Wiberg-Mayer bond orders
!      
   if (index(adum,'Mayer Atomic Bond Orders:').ne.0) then
      ibl=nat/6
      if (ibl*6.ne.nat) ibl=ibl+1
      do j=1,ibl
         read(iuout,'(a)',end=100)adum
         call readl(adum,xx,nn)
         do i=1,nat
            read(iuout,'(11x,6F11.6)') (wbo(idint(xx(kk)),i),kk=1,nn)
         end do
      end do
   end if
!
!     Read the hessian matrix     
! 
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
!
!     try to read the specific fchk file   
!   
if (.not.hread .and. .not. check_coord) then
   do i=1,len(fname)
      if(fname(i:i).ne.' ')j=i
   end do
   fname(j-2:j-2)='f'
   fname(j-1:j-1)='c'
   fname(j  :j  )='h'
   fname(j+1:j+1)='k'
   write(*,*) 'trying to read ',trim(fname)
   open(33,file=fname)
   do
      read(33,'(a)',end=200)adum
      if (index(adum,'Cartesian Force Constants').ne.0) then
         read(33,'(5E16.8)')r(1:n*(n+1)/2)
         exit
      end if
   end do
   200 close(33)
   k=0
   do i=1,n
      do j=1,i
         k=k+1
         h(j,i)=r(k)
         h(i,j)=r(k)
      end do
   end do
end if
1003 format(i7,5d14.6)
end subroutine gaurd
