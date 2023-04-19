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
!     subroutine hfit: fit QMDFF hessian to reference hessian and determine 
!            force constants
!
!     part of QMDFF
!
subroutine hfit(n,nb,at,xyz,q,r00ab,zab,r0094,sr42,c6xy, &
  &     h,torsfit,hneglect)
use qmdff
implicit none
integer::elem 
integer::n,at(n),nb(20,n)
real(kind=8)::xyz(3,n)
logical::torsfit
real(kind=8)::r00ab(94,94),zab(94,94),r0094(94,94)
real(kind=8)::sr42(94,94),c6xy(n,n),q(n)
real(kind=8)::h(3*n,3*n)
real(kind=8)::hneglect   
real(kind=8)::step,det,sumx,inkre,s0,s1,thr,dsq,de1,de2,dum
real(kind=8)::we,erep,mad,escal,mem
integer::n3,i,j,k,ic,jc,ia,ja,ii,jj,iiter,ns,js,is,kk
integer::L,np,errorflag,iii,jjj,nproc,ise(2,99)
integer::TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
character(len=20)::atmp
logical::ex
integer,allocatable::hessatoms(:,:)
integer,allocatable::iadr(:,:)
real(kind=8),allocatable::fx(:),y(:),pp(:),dp(:)
real(kind=4),allocatable::Nmat(:,:),tt(:,:),NM(:,:),aa(:,:),aaa(:)

write(10,*)
write(10,*)'======================='
write(10,*)'Levenberg-Marquardt fit'
write(10,*)'======================='
write(10,*)
call getpar0(np)
write(10,'(I5,A)') np,' parameters to determine'

CALL OMP_SET_DYNAMIC(.false.)
!$OMP PARALLEL PRIVATE(TID)
TID = OMP_GET_THREAD_NUM()
IF (TID .EQ. 0) THEN
   nproc = OMP_GET_NUM_THREADS()
   write(10,*) ' # OMP threads =', nproc
END IF
!$OMP END PARALLEL 
if(nproc.gt.np) stop ' # parameters < # procs!'


!
!     determine number of hessian elements that remain
!
kk=0
do i=1,n*3
   do j=1,i 
      dum=(h(i,j)+h(j,i))*0.5
      if (abs(dum).gt.hneglect) then
      kk=kk+1
      end if
   end do
end do
n3=kk              

j=3*n
i=j*(j+1)/2
write(10,'(I6,A)') n3,' data points'
write(10,'('' This is'',f6.2,'' % of the full Hessian'')')100.*dble(n3)/dble(i)
allocate(y(n3),fx(n3),iadr(2,n3))

kk=0
do i=1,n*3
   do j=1,i 
      dum=(h(i,j)+h(j,i))*0.5
      if (abs(dum).gt.hneglect) then
      kk=kk+1
      y(kk)=dum
      iadr(1,kk)=i
      iadr(2,kk)=j
      end if
   end do
end do
!
!     damp 
!     
inkre=0.01
!
!     exit
!              
thr=1.d-5
!
!     calculate memory for hessian
!
mem=2*float(n3)*float(np)+float(np)*float(np)
write(10,'(A,f11.6)')' --> memory requirements (Gb)',4*mem/1024**3
write(10,*)
allocate(hessatoms(5,np),pp(np),dp(np))
call procload(nproc,np,0,0.0d0,ise)
call getparini(pp,.false.)
!
!     parameter loop for num gradient
!     because the Hessian is linear in the FC,
!     the matrix aa is constant in the LM iterations!
!     --> one must calculate derivatives only at the beginning!
!
write(10,*)'computing dH/dFC ...'
!
!     calculate derivatives of hessian after force constants
!     for Levenberg-Marquardt step
!
!!!!$OMP PARALLEL private (i) copyin ( /ffdata/ ) 
!$OMP PARALLEL private (i) 
!$OMP DO
do i=1,nproc
   call pderiv(ise(1,i),ise(2,i),n,at,xyz,q,r00ab,zab,r0094,sr42,c6xy, &
      &  n3,np,iadr)
end do
!$OMP END DO
!$OMP END PARALLEL


allocate(aa(n3,np),aaa(n3))
allocate(tt(np,n3))
!
!     open file with unformatted derivative matrix
!
do tid=0,nproc-1
   if (tid.lt.100) write(atmp,'(''qmdfftmp.'',i2)')tid   
   if (tid.lt.10 ) write(atmp,'(''qmdfftmp.'',i1)')tid   
   open(unit=44,file=atmp)!,form='unformatted')
   do i=ise(1,tid+1),ise(2,tid+1)
      read(44,*)aaa
!      write(*,*) "aa", aaa
      aa(1:n3,i)=aaa(1:n3)
   end do
   close(44,status='delete')
end do


do i=1,n3
   do j=1,np
      tt(j,i)=aa(i,j)
   end do
end do
!
!     (*Matrizenmultiplikation *)
!     (* T^*A^ = Nmat          *)
!
allocate(NM(np,np))
call sgemm('T','T',np,np,n3,1.0e0,aa,n3,tt,np,0.0e0,NM,np)
deallocate(aa,aaa)

allocate(Nmat(np,np))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     iterations for Levenberg Maquardt!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
s1=1.d+42
write(10,*)
write(10,*) "Enter Levenberg Marquardt iterations:"
write(10,*)
write(10,*) "-------------------------------------------------"
do iiter=1,50
!
!     calculate QMDFF hessian with actual parameters
!
   call ffhess(n,at,xyz,q,r00ab,zab,r0094,sr42,c6xy,n3,iadr, &
       &   hessatoms,fx,0)

   s0=s1
   s1=0 
!
!     calculate change of error between reference and FF
!         
   do i=1,n3
      s1=s1+(y(i)-fx(i))**2
!     write(*,*) y(i)-fx(i),y(i),fx(i)
   end do
   dsq=s1-s0

   if (dsq.lt.0) then
      inkre=inkre/1.5
      inkre=min(inkre,1.0d0)
   else
      inkre=0.2          
   end if
   write(10,'('' iter '',i4,'' RMSD ='',F12.6,'' change '',F12.4)') &
     &  iiter,s1,dsq/s1
   if (abs(dsq/s1).lt.thr.or.s1.lt.thr) exit
   Nmat = NM
!
!     damping of diagonal
!
   do i=1,np
      Nmat(i,i)=Nmat(i,i)+inkre*Nmat(i,i)
   end do
!
!     (*  T^ * residuen = dp*)
!
   do i = 1, np 
      sumx = 0
      do L = 1, n3 
         sumx = sumx+tt(i, L)*(y(L)-fx(L))
      end do
      dp(i) = sumx
   end do
!
!     diagonalization of matrix..
!
   call spotrf('U',np, Nmat, np, errorflag)
   if (errorflag.ne.0) then
      write(10,*) 'matrix singular!!!'
      write(10,*) 'EXIT from hfit'
      goto 999
   end if
   call spotri('U',np, Nmat, np, errorflag)
   if (errorflag.ne.0) then
      write(10,*) 'matrix singular!!!'
      write(10,*) 'EXIT from hfit'
      goto 999
   end if

   do i = 1, np 
      do L = 1, np 
         Nmat(L, i)=Nmat(i,L)
      end do
   end do
!
!      (* N^-1 * (transp.*residuen) *)
!
   do i = 1, np 
      sumx = 0
      do L = 1, np 
         sumx = sumx+Nmat(i, L)*dp(L)
      end do
!
!      (* new Parameter *)
!
      pp(i)=pp(i)+sumx
   end do

   call putpar(pp,torsfit)
end do

99  continue

call putpar(pp,torsfit)

999 continue 
deallocate(tt,Nmat,NM,dp,pp)
write(10,*) "-------------------------------------------------"
write(10,*)
write(10,*) 'Fit of force constants to hessian matrix completed!'
if (iiter.gt.15) &
    & call warn('large number of iterations required. consider -hthr')

return
end subroutine hfit
