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
!     Subroutine opt_irc: optimize the minimum energy reaction for a 
!         chemical reaction starting from its transition state calculated
!         before in opt_ts.f90. The simple Euler method is used, more
!         sophisticated methods might be calculated in later versions...
!
!     part of EVB
!
subroutine opt_irc
use general
use evb_mod
implicit none
!     loop indices
integer :: i,j,k,l
!     local structure array
real(kind=8) :: xyz2(3,natoms)
!     the actual structure
real(kind=8) :: act_struc(3,natoms)
!     the actual gradient of structure 
real(kind=8) :: grad(3,natoms)
!     potential energy of structure (and of TS)
real(kind=8) :: epot,e_ts
!     potential energy of last structure 
real(kind=8) :: e_last
!     distances to TS of actual and next structure 
real(kind=8) :: dist_act,dist_next
!     the actual step (small)
real(kind=8) :: act_step(3,natoms)
!     the actual value of the path coordinate
real(kind=8) :: s_act
!     the absolute value of the actual gradient
real(kind=8) :: grad_abs
!     the cartesian hessian matrix
real(kind=8) :: hess(3*natoms,3*natoms)
!     the cartesian mass weighted eigenvalues
real(kind=8) :: eigval_cart(3*natoms)
!     the cartesian mass weighted eigenvectors
real(kind=8) :: eigvecs(3*natoms,3*natoms)
!     the intrinsic small euler step to guarantee a save ongoing 
!     of the algorithm for the Euler method
real(kind=8) :: euler_step 
!     how many Euler steps shall be done for one IRC step?
integer :: num_euler
!     total maximum number of gradient steps per side
integer :: steps_total
!     if the TS was too bad: steps needed for first correction
integer :: steps_corr
!     Arrays with resulting structures for forward/backward IRC
real(kind=8) :: irc_forward(irc_maxstep,3,natoms)
real(kind=8) :: irc_backward(irc_maxstep,3,natoms)
!     Arrays with resulting energies for forward/backward IRC
real(kind=8) :: ens_forward(irc_maxstep),ens_backward(irc_maxstep)
!     Arrays with important parameter values for forward/backward IRC
real(kind=8) :: s_forward(irc_maxstep),s_backward(irc_maxstep)
!     Number of entries for printout of arrays
integer :: steps_forward,steps_backward
!     Read in the start structure from global array
!     and convert it into bohr
!
do i=1,natoms
   xyz2(1,i)=x(i)/bohr
   xyz2(2,i)=y(i)/bohr
   xyz2(3,i)=z(i)/bohr
end do
write(*,*) "STEP2: Starting the IRC optimization..."
write(*,*)
write(15,*) " --- PART (B): IRC-OPTIMIZATION: --- "
write(15,*)

!
!     calculate actual cartesian hessian as starting point 
!
call hessevb(xyz2,hess)
!
!     Calculate mass-weighted eigenvalues and eigenvectors for 
!     first step of IRC!
!
call calc_freq(hess,eigval_cart,eigvecs,.false.)
!
!     FIRST POSSIBILITY: use the simple Euler method!
!
if (irc_method .eq. "euler") then
   write(15,*) "The simple Euler method will be used for numerical integration."
!
!     set hard coded intrinsic stepsize: the user defined step will 
!     be a multiplicity of it!
!   
   euler_step=irc_eulerlen
!
!     calculate total number of steps and how many Euler steps are
!     needed per full IRC step
!
   steps_total=int(real(irc_maxstep)/euler_step*irc_steplen)
   num_euler=int(irc_steplen/euler_step)
!
!     Calculate potential energy of starting point
!
   call gradient(xyz2,e_ts,grad,1)

!
!     Do big loop with two steps: follow both directions of the IRC
!     k=1: along the imaginary eigenmode
!     k=2: along the reversed imaginary eigenmode
!
   do k=1,2
      act_struc=xyz2
      e_last=0.d0
!
!     First step: Go along (or reverse) the imaginary eigenmode
!   
      if (k .eq. 1) then
         write(15,*) 
         write(15,*) "(B-1) Follow the forward direction of the reaction path!"
         write(15,*)
         do i=1,natoms
            do j=1,3
               act_step(j,i)=eigvecs((i-1)*3+j,1)
               s_act=euler_step
            end do
         end do
      else 
         write(15,*)
         write(15,*) "(B-2) Follow the backward direction of the reaction path!"
         write(15,*)
         do i=1,natoms
            do j=1,3
               act_step(j,i)=-eigvecs((i-1)*3+j,1)
               s_act=-euler_step
            end do
         end do            
      end if
!
!     calculate abolute value of the gradient and update structure
!
      grad_abs=sqrt(dot_product(act_step(1,:),act_step(1,:))+ &
                  & dot_product(act_step(2,:),act_step(2,:))+ &
                  & dot_product(act_step(3,:),act_step(3,:)))

!
!     write information about starting point
!
      write(15,*) "POINT          S(Q)-VALUE            ENERGY              GRADIENT-NORM"
      write(15,'(i6,a,f18.10,a,e18.10,a,e18.10)') 1,"    ",0.d0,"    ",e_ts,"    ",grad_abs


!
!     If the TS wasn't located exactly (numerical errors), the first step
!     could be too small in order to get in the right direction of the IRC; 
!     then, the energy would rise after the first step
!     Enlarge the size of the first step, until the energy decays!
!
      steps_corr=0
      do 
         act_struc=act_struc-euler_step*act_step/grad_abs
         call gradient(act_struc,epot,grad,1)
!
!     if goal was archived, leave the loop
! 
         if (e_ts-epot .gt. 0.d0) exit
!
!     if too many steps were needed, terminate program
!
         if (k .eq. 1) then
            s_act=s_act+euler_step
         else 
            s_act=s_act-euler_step
         end if
!
!     if number of correction steps is larger than one big Euler step, store 
!     needed data
!
         if (mod(steps_corr+2,num_euler) .eq. 0) then
            write(15,'(i6,a,f18.10,a,e18.10,a,e18.10)') (steps_corr+2)/num_euler+1,"(C) ",s_act,& 
                   & "    ",epot,"    ",grad_abs
            if (k .eq. 1) then
               irc_forward(int((steps_corr+2)/num_euler),:,:)=act_struc
               steps_forward=(steps_corr+2)/num_euler
               ens_forward(int((steps_corr+2)/num_euler))=epot
               s_forward(int((steps_corr+2)/num_euler))=s_act
            else
               irc_backward(int((steps_corr+2)/num_euler),:,:)=act_struc
               steps_backward=(steps_corr+2)/num_euler
               ens_backward(int((steps_corr+2)/num_euler))=epot
               s_backward(int((steps_corr+2)/num_euler))=s_act
            end if
         end if

         steps_corr=steps_corr+1
         if (steps_corr .gt. 100) then
            write(15,*) "The first step of the IRC couldn't be done for this system!"
            write(15,*) "It seems that the TS optimization archieved a bad outcome!"
            write(15,*) "Restart the calculation with new settings!"
            write(*,*) "The first step of the IRC couldn't be done for this system!"
            write(*,*) "It seems that the TS optimization archieved a bad outcome!"
            write(*,*) "Restart the calculation with new settings!"
            call fatal
         end if
      end do
      if (steps_corr .gt. 0) then
         write(15,'(i3,a)') steps_corr,"  correction steps were needed to push to the right &
                          & side of the TS!"
      end if
!
!     perform the loop over up to irc_maxstep steps in order to determine 
!     the remaining steps of the IRC (first step was eigenvector)
!
      do i=2+steps_corr,steps_total
!
!     calculate actual distance to TS
!
         dist_act=0.d0
         do j=1,natoms
            do l=1,3 
               dist_act=dist_act+(xyz2(l,j)-act_struc(l,j))*(xyz2(l,j)-act_struc(l,j))
            end do
         end do
         dist_act=sqrt(dist_act)
!
!     calculate gradient of actual structure
!         
         call gradient(act_struc,epot,grad,1)
!
!     mass-weight the gradient!--> actual step
!
         do j=1,natoms
            act_step(:,j)=grad(:,j)*sqrt(mass(indi(j)))
         end do
!
!     calculate norm of mass-weighted gradient vector
!
         grad_abs=sqrt(dot_product(act_step(1,:),act_step(1,:))+ &
                     & dot_product(act_step(2,:),act_step(2,:))+ &
                     & dot_product(act_step(3,:),act_step(3,:)))
!
!     update the actual structure 
!         
         act_struc=act_struc-euler_step*act_step/grad_abs
!
!     update the actual value of the reaction coordinate s
!
         if (k .eq. 1) then
            s_act=s_act+euler_step
         else 
            s_act=s_act-euler_step
         end if
!
!     
!
         if (mod(i,num_euler) .eq. 0) then
            write(15,'(i6,a,f18.10,a,e18.10,a,e18.10)') i/num_euler+1,"    ",s_act,"    ",&
                   & epot,"    ",grad_abs
            if (k .eq. 1) then
               irc_forward(int(i/num_euler),:,:)=act_struc 
               steps_forward=i/num_euler
               ens_forward(int(i/num_euler))=epot
               s_forward(int(i/num_euler))=s_act
            else 
               irc_backward(int(i/num_euler),:,:)=act_struc
               steps_backward=i/num_euler
               ens_backward(int(i/num_euler))=epot
               s_backward(int(i/num_euler))=s_act
            end if
         end if
!
!     Check if any of the convergene criteria is fulfilled!
! 
        ! write(15,*) epot,e_last,abs(epot-e_last)
         if (abs(epot-e_last) .lt. irc_ethr) then
            write(15,*) "-------------------------------------------------------"
            write(15,*) "Energy change fulfills convergence criterion!"
            write(15,'(a,e15.8,a,e15.8)') " DeltaE = ",abs(epot-e_last)," Criterion = ",irc_ethr
            write(15,*) "Convergence will be therefore signaled now for this branch!"
            exit
         end if       

!
!     Check if the minimum was reached and the optimization oscillates in it:
!     terminate optimization if distance to TS has grown with last step
!     Only, if the first 4 steps are finished
! 
         dist_next=0.d0
         do j=1,natoms
            do l=1,3
               dist_next=dist_next+(xyz2(l,j)-act_struc(l,j))*(xyz2(l,j)-act_struc(l,j))
            end do
         end do
         dist_next=sqrt(dist_next)
         if ((dist_next .le. dist_act) .and. (i .gt. 4)) then
            write(15,*) "-------------------------------------------------------"
            write(15,*) "It seems that the optimization has reached a minimum and"
            write(15,*) "is now oscillating in it! Therefore this branch will be"
            write(15,*) "will be stopped here and convergence is signaled!"
            exit
         end if

 
         e_last=epot

!
!
!     Print info message if all cycles are finished, but no convergence was reached
!
         if (i.eq.steps_total) then
            write(*,*) abs(epot-e_last) ,irc_ethr
            write(15,*) "-------------------------------------------------------"
            write(15,*) "WARNING: IRC-optimization hasn't converged for this branch!"
            write(15,*) "It will be written out, but maybe more than",steps_total," points"
            write(15,*) "along each branch will be needed!"
         end if

      end do

   end do
end if
write(15,*)
write(15,*) "IRC calculation FINISHED!"
write(*,*) "IRC calculation FINISHED!"

!
!     finally, write all IRC structures with TS in one coordinate file!
!
open(unit=23,file="irc.xyz",status="unknown")
!
!     first, all negative s values
!
do i=steps_backward,1,-1
   write(23,*) natoms; write(23,*)
   do j=1,natoms
      write(23,*) name(j),irc_backward(i,:,j)*bohr
   end do
end do
!
!     then the TS itself
!
write(23,*) natoms; write(23,*)
do j=1,natoms
   write(23,*) name(j),xyz2(:,j)*bohr
end do
!
!     last, all positive s values
!
do i=1,steps_forward
   write(23,*) natoms; write(23,*)
   do j=1,natoms
      write(23,*) name(j),irc_forward(i,:,j)*bohr
   end do
end do
close(23)
!
!     Info message in logfile
!
write(15,*)
write(15,*) "IRC structures written to irc.xyz"

!
!     Write also all IRC energies to file(s)!
!     Two files: irc_ens.dat (only energies) and irc_parens.dat
!       (energies and values of parameter s) 
!

open(unit=24,file="irc_ens.dat",status="unknown")
open(unit=25,file="irc_parens.dat",status="unknown")
!
!     first, all negative s values
!
do i=steps_backward,1,-1
   write(24,*) ens_backward(i)
   write(25,*) s_backward(i),ens_backward(i)
end do
!
!     then the TS itself
!
write(24,*) e_ts
write(25,*) 0.d0,e_ts
!
!     last, all positive s values
!
do i=1,steps_forward
   write(24,*) ens_forward(i)
   write(25,*) s_forward(i),ens_forward(i)
end do
close(24)
close(25)
!
!     Info message in logfile
!
write(15,*)
write(15,*) "IRC energies written to irc_ens.dat"

write(15,*)
write(15,*) "IRC energies with s-values written to irc_parens.dat"



return
end subroutine opt_irc
