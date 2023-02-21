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
!     Subroutine geoopt: optimize a local minimum of a structure using 
!     the BFGS algorithm
!
!     part of EVB
!

subroutine geoopt(coord)
use evb_mod
use general
                          
implicit none
!integer at(natoms)
real(kind=8)::h(3*natoms,3*natoms)
real(kind=8)::epot
real(kind=8)::coord(3,natoms)
real(kind=8)::grd(3,natoms)
real(kind=8)::g_dummy(3,natoms)
logical::restart
integer np,errorflag,i,j,k,icycle,nmax,lwork,liwork,nvar,nvar1
integer ifound,info,lina,ii,jj,ia,ja,ic,jc,nline,m
real(kind=8) ,allocatable :: hesinv(:)        
real(kind=8) ,allocatable :: d4(:)        
real(kind=8) ,allocatable :: g4(:)        
real(kind=8) ,allocatable :: gold(:)        
real(kind=8) ,allocatable :: xvar(:),gvar(:),xparam(:),xlast(:),gg(:)
real(kind=8) alpha,eold,f,emin,alp,el,xyz2(3,n),yhy,sy,ang,pnorm
real(kind=8) ggi,xvari,ddot,gnorm,alp0
character(len=2) asym
logical::converged
!
!     max. # of steps
!
nmax=maxiter
!
!     number of line searches
!
nline=10
!
!     start value for alpha
!
alp0=0.1


restart=.false.
alpha=0
ang  =0
np=3*natoms
nvar=np

allocate(d4(np),g4(np),gold(np),xparam(np),hesinv(np*(np+1)/2), &
   & xvar(np),gvar(np),xlast(np),gg(np))

!
!     initialization: calculate hessian and invert it for starting    
!

call hessevb(coord,h)
ii=0
hesinv=0
do i=1,np
   ii=ii+i
   hesinv(ii)=1./h(i,i)
enddo

epot=0
g4  =0
xparam=0
xlast =0
yhy   =1
open(unit=47,file="geoopt.xyz",status="unknown")
!
!     Enter the optimization loop
!
converged=.false.
write(15,*) "Do the geometry optimization with the BFGS algorithm!"
write(15,*) "Enter loop till convergence.."
do icycle=1,nmax
   eold=epot
!
!     calculate energy and gradient      
!
   
   call gradient(coord,epot,grd,1)
!
!     Write current structure to trajectory file
!
   write(47,*) nats
   write(47,*)
   do i=1,nats
      write(47,*) name(indi(i)),coord(:,(indi(i)))*bohr
   end do

  ! if(mod(icycle,10).eq.0)then
  !    write(85,*)n
  !    write(85,*) epot
  !    do i=1,n
  !       write(85,'(1x,a2,3F12.4)')asym(at(i)),coord(1:3,i)*0.52917726
  !    end do
  ! end if
   gnorm=sum(abs(grd))
!
!     normal exit if convergence was reached
!
   if (abs(eold-epot).lt.ethr .and. gnorm.lt.gthr &
         &  .and. sqrt(yhy).lt.dthr) then
      converged=.true.
      exit
   end if
!
!     absolutely no progress (optimization stagnates)
!
   if(sqrt(yhy).lt.1.d-6) then
       write(15,*) "The geometry optimization seems to be stagnating!"
       write(15,*) "Since no further progress can be made, the optimization"
       write(15,*) "will be ended here..."
       converged=.true.
       exit
   endif

   k=0
   do i=1,n
      do j=1,3
         k=k+1
         gold(k)=g4(k)
         g4(k)=grd(j,i)
         xvar(k)=xparam(k)-xlast(k)
         gvar(k)=g4(k)-gold(k)
      end do
   end do
!
!     Print info for this cycle
!
   if (mod(icycle,10).eq.0.or.icycle.eq.1) then
      write(15,'('' iteration : '',i4,''  E='',F14.8, &
          &      ''  Gnorm='',F8.5,  &
          &      ''  displacement='',f8.5, &
          &      ''  alpha='',F6.2, &
          &      ''  angle='',F6.2)')  &
          &      icycle,epot,gnorm,sqrt(yhy),alpha,ang
   end if
!
!     update the BFGS every icycle steps!
!
   if (icycle.gt.1) then
      call supdot(gg,hesinv,gvar,nvar,1)
      yhy=ddot(nvar,gg,1,gvar,1)
      sy =ddot(nvar,xvar,1,gvar,1)+1.d-12
      yhy=1.0d0 + yhy/sy
      k=0
      do i=1,nvar
         xvari=xvar(i)/sy
         ggi=gg(i)/sy
         do j=1,i
            k=k+1
            hesinv(k)=hesinv(k)-gg(j)*xvari-xvar(j)*ggi + & 
               &     yhy*xvar(j)*xvari
         end do
      end do
      call supdot(d4,hesinv,g4,nvar,1)
   else
      d4 = g4
   end if
!
!     line search
!
   alp=alp0
   alpha=alp
   emin=1.d+42
   do m=1,nline
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            xyz2(j,i)=coord(j,i)-d4(k)*alp
         end do
      end do
!
!     energy only     
! 
      call gradient(xyz2,el,g_dummy,1)   
      if (el .lt. emin) then
         alpha=alp
         emin=el   
      end if
      alp=alp+1./dble(nline)
   end do


   pnorm=sqrt(ddot(nvar,alp*d4,1,alp*d4,1))
   gnorm=sqrt(ddot(nvar,g4,1,g4,1))
   ang=ddot(nvar,alp*d4,1,g4,1)/(pnorm*gnorm)

   if(ang.lt.-0.2)then
      ii=0
      hesinv=0
      do i=1,np
         ii=ii+i
         hesinv(ii)=1.
      end do
      call warn('resetting inverse Hessian')
      d4 = g4
      alpha=0.02
   end if

   k=0
   do i=1,n
      do j=1,3
         k=k+1
         xlast(k)=coord(j,i)
         coord(j,i)=coord(j,i)-d4(k)*alpha
         xparam(k)=coord(j,i)
      end do
   end do

   yhy=ddot(nvar,d4,1,d4,1)

end do
if (.not. converged) then
   write(*,*) "WARNING: The geometry optimization did not reached convergence" 
   write(*,*) " in the given steps! Increase MAXITER!"
   write(15,*) "WARNING: The geometry optimization did not reached convergence"
   write(15,*) " in the given steps! Increase MAXITER!"
   
else 
   write(15,*) "The geometry optimization has converged!"
end if
!
!     last printout after the optimization is done 
!
write(15,'('' iteration : '',i4,''  E='',F14.8, &
      &    ''  Gnorm='',F8.5, &
      &    ''  displacement='',f8.5, &
      &    ''  alpha='',F6.2)') &
      &    icycle,epot,gnorm,sqrt(yhy),alpha

close(85)

deallocate(d4,g4,gold,xparam,hesinv,xvar,gvar,xlast,gg)

open(unit=49,file="opt_final.xyz",status="replace")
write(49,*) nats
write(49,*)
do i=1,nats
   write(49,*) name(indi(i)),coord(:,(indi(i)))*bohr
end do
close(49)


end subroutine geoopt


