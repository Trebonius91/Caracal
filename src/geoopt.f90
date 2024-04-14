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
!     Subroutine geoopt: optimize a local minimum of a structure using 
!     the BFGS algorithm or the CG algorithm!
!
!     part of EVB
!

subroutine geoopt(coord)
use evb_mod
use general
use pbc_mod
                          
implicit none
!integer at(natoms)
real(kind=8)::h(3*natoms,3*natoms)
real(kind=8)::epot
real(kind=8)::coord(3,natoms),coord_old(3,natoms)
real(kind=8)::grd(3,natoms)
real(kind=8)::g_dummy(3,natoms)
logical::restart
logical::conv_ethr,conv_gthr,conv_dthr,conv_gmaxthr,conv_dmaxthr
integer np,errorflag,i,j,k,l,icycle,nmax,lwork,liwork,nvar,nvar1
integer ifound,info,lina,ii,jj,ia,ja,ic,jc,nline,m
real(kind=8) ,allocatable :: hesinv(:)        
real(kind=8) ,allocatable :: d4(:)        
real(kind=8) ,allocatable :: g4(:)        
real(kind=8) ,allocatable :: g_old(:)        
real(kind=8) ,allocatable :: search_dir(:)
real(kind=8) ,allocatable :: xvar(:),gvar(:),xparam(:),xlast(:),gg(:)
real(kind=8) alpha,e_old,f,emin,alp,el,xyz2(3,natoms),yhy,sy,ang,pnorm
real(kind=8) ggi,xvari,ddot,gnorm,alp0,step_norm,gmax_act,dmax_act,de_act
real(kind=8),allocatable::grd_1d(:),grd_old_1d(:)
real(kind=8) :: cell_mat(3,3),cell_mat_inv(3,3)  ! for VASP files
real(kind=8) :: x_frac,y_frac,z_frac  ! direct coordinates (CONTCAR)
character(len=2) asym
logical::converged
logical::act_fix
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

if (geoopt_algo .eq. "bfgs") then

   restart=.false.
   alpha=0
   ang  =0
   np=3*natoms
   nvar=np

   allocate(d4(np),g4(np),g_old(np),xparam(np),hesinv(np*(np+1)/2), &
      & xvar(np),gvar(np),xlast(np),gg(np))

!
!     initialization: calculate hessian and invert it for starting    
!
   write(*,*) "Approximate the Hessian for the initial structure ..."
 !  call hessevb(coord,h)
   h=0.5d0
   ii=0
   hesinv=0.d0
   do i=1,np
      ii=ii+i
      hesinv(ii)=1./h(i,i)
   end do

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
   write(*,*) "------------------------------------------------------&
                &-------------------------------------------"
   write(*,*) "  Step        Etot         DeltaE        grad. norm  &
               &   step norm      grad. max      step max"
   write(*,*) "------------------------------------------------------&
                &-------------------------------------------"

   do icycle=1,maxiter
      e_old=epot
!
!     calculate energy and gradient      
!
      call gradient(coord,epot,grd,1,1)
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
!      gnorm=sum(abs(grd))/natoms
!
!     normal exit if convergence was reached
!



!      if (abs(e_old-epot).lt.ethr .and. gnorm.lt.gthr &
!            &  .and. sqrt(yhy).lt.dthr) then
!         write(*,*) "converged"
!         converged=.true.
!         exit
!      end if

!
!     Check for convergence, all five criteria
!
!     A: energy change
!
      conv_ethr = .false.
      de_act=abs(epot-e_old)
      if (de_act .lt. ethr) conv_ethr = .true.
!
!     B: gradient norm (per atom)
!
      gnorm=0.d0
      do i=1,n
         do j=1,3
            k=k+1
            gnorm=gnorm+grd(j,i)*grd(j,i)
         end do
      end do
      gnorm=sqrt(gnorm)/natoms
      conv_gthr = .false.
      if (gnorm .lt. gthr) conv_gthr = .true.
!
!     C: geometry step norm (per atom)
!      
      step_norm=0.d0
      do i=1,n
         do j=1,3
            k=k+1
            step_norm=step_norm+(coord(j,i)-coord_old(j,i))**2
         end do
      end do
      conv_dthr = .false.
      if (step_norm .lt. dthr) conv_dthr = .true.
!
!     D: gradient maximum component
!
      conv_gmaxthr = .false.
      gmax_act=maxval(abs(grd))
      if (gmax_act .lt. gmaxthr) conv_gmaxthr = .true.
!
!     E: geometry step largest component
!
      conv_dmaxthr = .false.
      dmax_act=maxval(abs(coord-coord_old))
      if (dmax_act .lt. dmaxthr) conv_dmaxthr = .true.

      write(*,'(i6,f17.9,a,es10.3,a,L1,a,es10.3,a,L1,a,es10.3,a,L1,a, &
             & es10.3,a,L1,a,es10.3,a,L1,a)') icycle,epot," ",de_act," (",conv_ethr, &
             &   ") ",gnorm," (",conv_gthr,") ",step_norm," (",conv_dthr,") ", &
             & gmax_act," (",conv_gmaxthr,") ",dmax_act," (",conv_dmaxthr,")"
!
!     If all criteria are fulfilled, declare the optimization as successful
!      
      if (conv_ethr .and. conv_gthr .and. conv_dthr .and. conv_gmaxthr &
           & .and. conv_dmaxthr) then
         converged = .true.
         write(*,*) "------------------------------------------------------&
                &-------------------------------------------"

         write(*,*) "Hurray! The geometry optimization has converged!"
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
            g_old(k)=g4(k)
            g4(k)=grd(j,i)
            xvar(k)=xparam(k)-xlast(k)
            gvar(k)=g4(k)-g_old(k)
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
         call gradient(xyz2,el,g_dummy,1,1)   
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
      coord_old=coord
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

   deallocate(d4,g4,g_old,xparam,hesinv,xvar,gvar,xlast,gg)
!
!     The conjugate gradient algorithm
!
else if (geoopt_algo .eq. "cg") then

   write(*,*) "------------------------------------------------------&
                &-------------------------------------------"
   write(*,*) "  Step        Etot         DeltaE        grad. norm  &
               &   step norm      grad. max      step max"
   write(*,*) "------------------------------------------------------&
                &-------------------------------------------"

   converged=.false.
   np=3*natoms
   allocate(g_old(np),search_dir(np))
   allocate(grd_1d(np))
!
!     Initialization: calculate initial gradient and set search direction to -g
!
   call gradient(coord,epot,grd,1,1)
    
   k=0
   do i=1,n
      do j=1,3
         k=k+1
         search_dir(k)=-grd(j,i)
      end do
   end do

!
!     Inter the main optimization loop
!
   do icycle=1,maxiter
!
!     Store old energy and gradient and structure (last step)
!
      e_old=epot
      k=0
      do i=1,n 
         do j=1,3 
            k=k+1 
            g_old(k)=grd(j,i)
         end do
      end do
      coord_old=coord

!
!     Perform the line search along the current direction
! 
      alpha=0.3d0
      do l=1,nline
!
!     Update the coordinates
!
         k=0
         do i=1,n
            do j=1,3
               k=k+1
               xyz2(j,i)=coord(j,i)+search_dir(k)*alpha
            end do
         end do
!
!     Calculate new energy 
!          
         call gradient(xyz2,el,g_dummy,1,1)
!
!     Check if energy decreased sufficiently
!
         if (el < e_old) then
            exit
         else 
            alpha=alpha/2.d0
         end if
      end do
!
!     Update the coordinates
!
      coord=xyz2
!
!     Calculate new gradient
!
      call gradient(coord,epot,grd,1,1)

!
!     Calculate beta factor with Polack-Ribiere formula
!
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            grd_1d(k)=grd(j,i)
         end do
      end do
      beta=dot_product(grd_1d,grd_1d-g_old) / dot_product(g_old,g_old)
!
!     Update the srach direction
!
      k=0
      do i=1,n
         do j=1,3
            k=k+1
            search_dir(k)=-grd(j,i)+search_dir(k)*beta
         end do
      end do

!
!     Check for convergence, all five criteria
!
!     A: energy change
!
      conv_ethr = .false.
      de_act=abs(epot-e_old)
      if (de_act .lt. ethr) conv_ethr = .true.  
!
!     B: gradient norm (per atom)
!
      gnorm=0.d0
      do i=1,n
         do j=1,3
            k=k+1
            gnorm=gnorm+grd(j,i)*grd(j,i)
         end do
      end do
      gnorm=sqrt(gnorm)/natoms
      conv_gthr = .false.
      if (gnorm .lt. gthr) conv_gthr = .true.
!
!     C: geometry step norm (per atom)
!      
      step_norm=0.d0
      do i=1,n
         do j=1,3
            k=k+1
            step_norm=step_norm+(coord(j,i)-coord_old(j,i))**2
         end do
      end do
      conv_dthr = .false.
      if (step_norm .lt. dthr) conv_dthr = .true. 
!
!     D: gradient maximum component
!
      conv_gmaxthr = .false.
      gmax_act=maxval(abs(grd))
      if (gmax_act .lt. gmaxthr) conv_gmaxthr = .true.
!
!     E: geometry step largest component
!
      conv_dmaxthr = .false.
      dmax_act=maxval(abs(coord-coord_old))
      if (dmax_act .lt. dmaxthr) conv_dmaxthr = .true.

      write(*,'(i6,f17.9,a,es10.3,a,L1,a,es10.3,a,L1,a,es10.3,a,L1,a, &
             & es10.3,a,L1,a,es10.3,a,L1,a)') icycle,epot," ",de_act," (",conv_ethr, &
             &   ") ",gnorm," (",conv_gthr,") ",step_norm," (",conv_dthr,") ", &
             & gmax_act," (",conv_gmaxthr,") ",dmax_act," (",conv_dmaxthr,")" 
!
!     If all criteria are fulfilled, declare the optimization as successful
!      
      if (conv_ethr .and. conv_gthr .and. conv_dthr .and. conv_gmaxthr &
           & .and. conv_dmaxthr) then
         converged = .true.
         write(*,*) "------------------------------------------------------&
                &-------------------------------------------"

         write(*,*) "Hurray! The geometry optimization has converged!"
         exit
      end if
   end do

!
!     Give final information about calculation
!
   if (.not. converged) then
      write(*,*) "WARNING: The geometry optimization did not reached convergence"
      write(*,*) " in the given steps! Increase MAXITER!"
      write(15,*) "WARNING: The geometry optimization did not reached convergence"
      write(15,*) " in the given steps! Increase MAXITER!"

   else
      write(15,*) "The geometry optimization has converged!"
   end if

end if
!
!     Print the final optimized geometry
!
open(unit=49,file="opt_final.xyz",status="replace")
write(49,*) nats
write(49,*)
do i=1,nats
   write(49,*) name(indi(i)),coord(:,(indi(i)))*bohr
end do
close(49)

!
!    In the case of VASP type input coordinates (POSCAR), 
!      also print out a CONTCAR file!
!

if (coord_vasp) then
   open(unit=50,file="CONTCAR",status="replace")
   write(50,*) "CONTCAR file written by Caracal (explore.x)"
   write(50,*) vasp_scale
   write(50,*) vasp_a_vec
   write(50,*) vasp_b_vec
   write(50,*) vasp_c_vec
   do i=1,nelems_vasp
      write(50,'(a,a)',advance="no") " ",trim(vasp_names(i))
   end do
   write(50,*)
   write(50,*) vasp_numbers(1:nelems_vasp)
   write(50,*) "Direct"
   if (vasp_selective) then
      write(50,*) "Selective dynamics"
   end if
!
!     As in usual CONTCAR files, give the positions in direct coordinates!
!     convert them back from cartesians, by using the inverse matrix
!
   cell_mat(1,:)=vasp_a_vec
   cell_mat(2,:)=vasp_b_vec
   cell_mat(3,:)=vasp_c_vec
   call matinv3(cell_mat,cell_mat_inv)
   do i=1,nats
      x_frac=dot_product(coord(:,indi(i))*bohr,cell_mat_inv(1,:))
      y_frac=dot_product(coord(:,indi(i))*bohr,cell_mat_inv(2,:))
      z_frac=dot_product(coord(:,indi(i))*bohr,cell_mat_inv(3,:))
      if (vasp_selective) then
         act_fix=.false.
         do j=1,fix_num
            if (fix_list(j) .eq. indi(i)) then
               act_fix=.true.
            end if
         end do
         if (act_fix) then
            write(50,*) x_frac,y_frac,z_frac,"   F   F   F "
         else 
            write(50,*) x_frac,y_frac,z_frac,"   T   T   T "
         end if
      else
         write(50,*) x_frac,y_frac,z_frac
      end if  
   end do
   close(50)
end if

end subroutine geoopt

