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
!     subroutine gradient: calculates the potential EVB-QMDFF-energy
!     and first derivatives with respect to Cartesian coordinates
!     NOTE (21.12.2016): for the evb_qmdff program, gradients in hartree/bohr
!     are needed for numerical hessians! Therefore no convertion to kcal/ang
!     is done here.
!     NOTE (07.06.2018): Now the structure is directly given as dummy argument
!
!     part of EVB
!
subroutine gradient (xyz2,e_evb,g_evb,act_bead)
use general
use evb_mod
implicit none
integer i,j,t,n2 
integer::act_bead  ! the actual bead number/index to calculate
real(kind=8)::xyz2(3,natoms),e_evb,g_evb(3,natoms),g_test(3,natoms) 
real(kind=8)::e1_shifted,e2_shifted,e3_shifted,e_two,gnorm_two
real(kind=8)::e_three,gnorm_three
real(kind=8)::ediff,offdiag,off4,root2,deldiscr,delsqrt
real(kind=8)::e,gnorm,e_qmdff1,e_qmdff2,e_qmdff3,delta_Q,delta_Q2
real(kind=8)::dE12,dE13,dE23,gnorm_evb,dVdQ,dq_sqrt
real(kind=8)::e_cov_local
!     For dQ
real(kind=8),dimension(:),allocatable::int_ts2,int_coord
!     For Sonennberg-Schlegel:
real(kind=8),dimension(3*natoms,3*natoms)::pC
real(kind=8),dimension(3*natoms,3*natoms)::d1d2,d2d1
real(kind=8)::b_dq,dq_c_dq,dq_dq,alph,coupl,pA
real(kind=8),dimension(3*natoms)::dq,dq2,sosc_dq
real(kind=8),dimension(3*natoms)::pB,pC_dq
real(kind=8)::e_lower,e_upper,step,str_e1,str_e2
integer::nat,nat3,k,l,m,o
real(kind=8),dimension(:,:),allocatable::U,Ut,g_mat,UtA
!     For dg_evb/DG-EVB
real(kind=8),dimension(nat6)::q_qts
real(kind=8),dimension(nat6)::int_path
real(kind=8)::V12,root,d_p,expo,V12_0
integer::dg_evb_mode,mat_size 
!     For rp_evb/RP-EVB
real(kind=8)::denom,diff,eterm,num,s_act,z_act
real(kind=8)::zeta   ! the damping function
real(kind=8),dimension(nat6+2)::path_s   ! interpolated strcutures/energies along the path
real(kind=8),dimension(rp_spl2_dim/2)::path_s2 ! interpolated gradient/hessians
real(kind=8),dimension(nat6)::int_min  ! coordinates of the nearest IRC point
integer::igh_left,igh_right  ! indices of nearest gradient/frequency ref. point
real(kind=8)::wgh_left,wgh_right  ! weights of gradients/frequencies 
real(kind=8)::dv12_act(nat6),d2v12_act(nat6,nat6)  ! the actual taylor coefficients
real(kind=8)::grad_act(nat6),hess_act(nat6,nat6)  ! the actual taylor coefficients
real(kind=8)::deltq(nat6),deltq2(nat6)  ! difference vector between IRC and actual point
real(kind=8)::dq_dv12_act(nat6),dq_hess_act(nat6)  ! intermediary result...
real(kind=8)::dvds  ! for ts region and energy interpolation
real(kind=8)::en_evb,en_ref   ! energies for EVB and reference expansion in mixing region
real(kind=8)::grad_int(nat6),grad_taylor(nat6)  ! for ts region and energy interpolation
real(kind=8)::grad_v_xyz(3,natoms)  ! cartesian gradient of energy interpolation
real(kind=8)::Vref   ! interpolated reference energy for TS region
real(kind=8)::s_upper,s_lower,s_slope ! TEST 06.09.2018
real(kind=8)::q_square   ! vector norm of the difference vector to IRC
real(kind=8)::V12_const  ! constant taylor expansion term for analytical derivative
real(kind=8)::grad_norm,hess_norm   !TEST
real(kind=8)::xyz_backup(3,natoms)   ! store the actual structure 
real(kind=8)::int_backup(nat6)   ! store the actual internal path projection
real(kind=8)::s_backup   ! store the actual s (progress) value
real(kind=8)::g_int(nat6)  ! TEST for internal gradient
real(kind=8)::x_val,smoothstep   ! for smoothstep function
real(kind=8)::g_V12_int(nat6)  ! internal gradient of the coupling term
real(kind=8)::g_V12(3,natoms)   ! cartesian gradient of the coupling term
real(kind=8)::apar,b,dV12ds,h   ! for interpolation gradient
real(kind=8)::dqds(nat6),dsdq(nat6)  ! for interpolation gradient (path direction)
real(kind=8)::abs_dqds  ! absolute value of the interpolation gradient (length)
real(kind=8)::dq_taylor(nat6)  ! for interpolation gradient (taylor derivative)
real(kind=8)::taylor_proj(nat6),numgrad_int(nat6)   ! TEST 17.09.2018
logical::unset_V12  ! if the coupling is negative, avoid NaN values of roots
logical::calc_evb,calc_ens  ! if EVB-QMDFF amd/or interpolation of energies shall be calculated
logical::rp_transition  ! if both switches from above are turned on
logical::trans_left,trans_right  ! if we are in the right or left part of the transtion
real(kind=8)::deltxyz(3*natoms)  ! the cartesian analogue to deltq
integer::khi,klo    ! for interpolation gradient
real(kind=8)::pi
character(len=100)::adum  ! for bug prevention
!    complex numbers for analytical DG-EVB gradient (V12*dV12)
complex(kind=16)::rootV,V_lower,V_upper,V12_comp,dV12_comp
real(kind=8),allocatable::alph1(:)
real(kind=8),allocatable::VdV_i(:),VdV_k(:) ! internal and cartesian coordinate vectors
!
!     variables for lapack-dsyev-Subroutine
!
character(len=1)::JOBZ,UPLO
integer::INFO,LDA,LWORK,Nd
!real(kind=8),dimension(:,:),allocatable::A,A1
real(kind=8),dimension(:),allocatable::W,WORK
real(kind=8),dimension(3,3)::evb_mat
!
!    Variables for analytical potential energy surface submissions
!
real(kind=8)::pot_geo(3,natoms,1)
real(kind=8)::pot_grad(3,natoms,1)
n=natoms
nat3=3*natoms
nat=natoms

pi=3.141592653589793238462643d0

!    TEST
!error_vec_ens=0.d0
!error_vec_v12=0.d0
! 
!     convert input-coordinates from angström to bohr 
!     for RPMD program do not convert! (all is in bohr)
!
if (use_rpmd) then
!
!     for ring polymers: read in the structure of the current bead
!
else 

end if
!stop "hggf"

!write(*,*) "max",maxval(xyz2(:,:))
!
!     If for test reasons an analytical potential shall be calculated, 
!     skip everything else 
!
if (pot_ana) then
   
!
!     The input array has one additional dimension for RPMD beads!
!     for usual egrad.x calculations, convert structures to bohr
!
!   if (use_rpmd) then
!      pot_geo(:,:,1)=xyz2
!   else 
!      pot_geo(:,:,1)=xyz2/bohr
!   end if
!

pot_geo(:,:,1)=xyz2

!
!     Call one of the different potential energy surfaces
!
   if (pot_type .eq. "mueller_brown") then
      call egrad_mueller(pot_geo,e_evb,pot_grad,info)
   else if (pot_type .eq. "ch4oh") then
      call egrad_ch4oh(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "ch4h") then
      call egrad_ch4h(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "nh3oh") then
      call egrad_nh3oh(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "ch4cn") then
      call egrad_ch4cn(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "h3") then
      call egrad_h3(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "brh2") then
      call egrad_brh2(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "o3") then
      call egrad_o3(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "clnh3") then
      call egrad_clnh3(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "oh3") then
      call egrad_oh3(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "geh4oh") then
      call egrad_geh4oh(pot_geo,natoms,1,e_evb,pot_grad,info)
   else if (pot_type .eq. "c2h7") then
      call egrad_c2h7(pot_geo,natoms,1,e_evb,pot_grad,info)
   end if
!
!     Convert back to smaller gradient array
!
      g_evb=pot_grad(:,:,1)

   return
end if

!
!     Invoke orca if ab-initio MD is to be used
!
if (orca) then
   call orca_grad(xyz2,pot_grad,e_evb)
   g_evb=pot_grad(:,:,1)
   return
end if


!
!     convert cartesian coordinates to bohr, if no analytical PES is used
!
!xyz2=xyz2/bohr

delta_Q=0 !if no dQ couplingterm is used
!
!     Calculate the energies and gradients of the first QMDFF
!
n2=natoms
if (.not. rp_evb) then
   call ff_eg(n_one,at,xyz2,e,g_one)
   if (energysplit) then
      e_cov_local=e
      e_cov_split=e_cov_split+e
   end if   
 !  write(*,*) "e1",e
   call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e,g_one)
!   write(*,*) "e2",e
   call ff_hb(n_one,at,xyz2,e,g_one)
   if (energysplit) then
      e_noncov_split=e_noncov_split+e-e_cov_local
   end if
end if
!
!     In case you have only one QMDFF, these are already the results
!
if (nqmdff.eq.1) then

   e_evb=e+E_zero1      
   g_evb=g_one

   return
!
!     In case of evb-qmdff, energy and gradients of the second QMDFF are determined
!
else if (nqmdff.eq.2) then
   if (.not. rp_evb) then
      call ff_eg_two(n_one,at,xyz2,e_two,g_two)
      call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,r094_mod,sr42,&
                 &c6xy_two,e_two,g_two)
      call ff_hb_two(n_one,at,xyz2,e_two,g_two)

      gnorm=sqrt(sum(g_one**2))
      gnorm_two=sqrt(sum(g_two**2))
   end if
!
!   Calculate shifted energies of the QMDFFs (corrected with ab-initio-energies)
!
   e1_shifted=e+E_zero1  
   e2_shifted=e_two+E_zero2  
   e_qmdff1=e1_shifted
   e_qmdff2=e2_shifted
   
!  Maybe this is needed for EVB optimizsation: abs(...)
   ediff=e1_shifted-e2_shifted
!
!   The Distributed Gaussian EVB coupling term
!
   if (dg_evb) then
!
!   calculate the gradient (analytical and mumerical):
!
      dg_evb_mode=dg_mode
      if (dg_evb_mode .eq. 1) then
         mat_size=dg_evb_points
      else if (dg_evb_mode .eq. 2) then
         mat_size = dg_evb_points*(1+nat6)
      else if (dg_evb_mode .eq. 3) then
         mat_size=dg_evb_points*(1+nat6+(nat6)*(nat6+1)/2)
      end if
   !   num_grad=.false.
      if (num_grad .eqv. .true.) then
!
!   The numerical DG-EVB gradient
!
         step=num_grad_step
         do i=1,natoms
            do j=1,3
               do k=1,2
                  if (k.eq.1) then
                     xyz2(j,i)=xyz2(j,i)+step
                  else if (k.eq.2) then
                     xyz2(j,i)=xyz2(j,i)-2d0*step
                  end if
                  call ff_eg(n_one,at,xyz2,e,g_one)
                  call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,&
                       &sr42,c6xy,e,g_one)
                  call ff_hb(n_one,at,xyz2,e,g_one)
                  call ff_eg_two(n_one,at,xyz2,e_two,g_two)
                  call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,&
                       &r094_mod,sr42,c6xy_two,e_two,g_two)
                  call ff_hb_two(n_one,at,xyz2,e_two,g_two)

                  call xyz_2int(xyz2,int_path,nat)

                  e1_shifted=e+E_zero1
                  e2_shifted=e_two+E_zero2
                  e_qmdff1=e1_shifted
                  e_qmdff2=e2_shifted
!
                  V12=0
 
!
!   call subroutine for calculating the DG-EVB-resonance integral 
!
                  call sum_v12(dg_evb_mode,mat_size,alph_opt,int_path,V12)
!
!   prohibit negative coupling-element values
!
                  root=(0.5d0*(e_qmdff1-e_qmdff2))*(0.5d0*(e_qmdff1-e_qmdff2))+V12
                  if (root .le. 0) then
                     e=0.5d0*(e_qmdff1+e_qmdff2)
                  else
                     e=0.5d0*(e_qmdff1+e_qmdff2)-sqrt(root)
                  end if
                  if (k.eq.1) then
                     e_upper=e
                  else if (k.eq.2) then
                     e_lower=e
                  end if
               end do
               xyz2(j,i)=xyz2(j,i)+step
               g_evb(j,i)=(e_upper-e_lower)/(2d0*step)
               g_evb(j,i)=g_evb(j,i)
            end do
         end do
      else 
!
!     The analytical DG-EVB gradient
!     Newly implemented: 22.07.2018!
!     Now the calculation is done complete analytically, like in 
!     the case of RP-EVB

         call xyz_2int(xyz2,int_path,nat)      
!
!     calculate the coupling strengh
!  
         call sum_v12(dg_evb_mode,mat_size,alph_opt,int_path,V12)
!
!     calculate the coupling gradient in internal coordinates!
!
         call sum_dv12(dg_evb_mode,mat_size,alph_opt,int_path,g_V12_int)

!
!     Convert the squared coupling gradient to cartesian coordinates
!
         call int2grad(xyz2,int_path,g_V12_int,g_V12)
!
!     Now, calculate the EVB-QMDFF gradient using the well known formula
!     Using some mathematical tricks, we can avoid the root of the coupling term!        
! 
         offdiag=V12
         off4=4.d0*offdiag
!
!     If the root argument becomes negative, set the coupling to zero..
!
         unset_V12=.false.
         if (ediff*ediff+off4 .lt. 0.d0) unset_V12=.true.
         if (unset_V12 .eqv..false.) then
            root2=sqrt(ediff*ediff+off4)
         end if
         do i=1,n_one
            do j=1,3
               deldiscr=ediff*(g_one(j,i)-g_two(j,i))+2.d0*g_V12(j,i)
               if (unset_V12 .eqv. .false.) then
                  delsqrt=deldiscr/root2
               else
                  delsqrt=0.d0
               end if
               g_evb(j,i)=0.5d0*(g_one(j,i)+g_two(j,i)-delsqrt)
            end do
         end do
      end if

!
!     Optional test: write components of internal gradients to file
!
      if (int_grad_plot) then
         g_int=0.d0
         call grad2int(xyz2,int_path,g_int,g_evb)
         write(192,*) sum(g_int**2)
      end if

!
!      calculate energy of the undistorted structure
!

      call xyz_2int(xyz2,int_path,nat)
      V12=0
!
!    Write internal coordinates to file if option is activated
!    Print a third column for splot an convert them into Angstrom
!
      if (int_coord_plot) then       
         write(191,*) int_path*bohr,0.0d0
      end if

      call sum_v12(dg_evb_mode,mat_size,alph_opt,int_path,V12)
      root=(0.5d0*(e_qmdff1-e_qmdff2))*(0.5d0*(e_qmdff1-e_qmdff2))+V12
      if (root .le. 0) then
         e_evb=0.5d0*(e_qmdff1+e_qmdff2)
      else
         e_evb=0.5d0*(e_qmdff1+e_qmdff2)-sqrt(root)
      end if 
!
!   jump directly to back to the main program
!
      return
   end if
!
!   The Reaction-Path(RP)-EVB coupling term
!   corrected with pure energy interpolation in the TS region
!
   if (rp_evb) then
   s_act=0.d0
z_act=0.d0
!
!   Convert structure to internal coordinates 
!
call xyz_2int(xyz2,int_path,natoms)
!   
!   Calculate the s-value along the IRC, which has the lowest distance 
!   to the coordinates of the given point 
!   Then, calculate the values for the internal coordinates and the 
!   coupling strength on the s-parameterized path for this s value
!
!   If the structure is outside the coupling region, set V12 and zeta to zero
!   to turn off the coupling
!

call spline_dist(nat6,s_act,int_path,int_min,act_bead)

!
!   If we are near the TS, calculate also/only the directly interpolated 
!   reference energy: determine here if it shall be done
!   Both variants can/should be active together for some regions, because 
!   a smooth transition is going to happen there
!
calc_ens=.false.
calc_evb=.false.
rp_transition=.false.
!
!     If the NO_EVB option is activated, activate always calc_ens, else,
!     check, which s-values are there
!
!     If the total interval of direct energy interpolation is too small, deactivate 
!     it completely in order to avoid bad behavior
!
if (.not. no_evb) then
   if (rp_mid_tot .lt. 0.02d0) then
      calc_evb=.true.
   else  
      if ((s_act .ge. trans_r_lo) .or. (s_act .le. trans_l_hi)) then
         calc_evb=.true.
      end if

      if ((s_act .ge. trans_l_lo) .and. (s_act .le. trans_r_hi)) then
         calc_ens=.true.
      end if
!
!     if we are in the transition region turn on local information boolean
!
      if (calc_ens .and. calc_evb) then
         rp_transition=.true.
      end if
   end if     
else 
!
!     for the NO_EVB option, activate the direct interpolation, unless we are in the 
!     pure QMDFF regions
!
   if ((s_act .ge. iq_l_lo) .and. (s_act .le. iq_r_hi)) then
      calc_ens=.true.
   end if
end if
!
!     calculate single QMDFF energies if EVB is activated
!     If no RPMD calculation shall be done, write out single QMDFF energies
!
      ! write(*,*) "s_act,z_act",s_act,z_act,calc_evb
if (calc_evb) then
   call ff_eg(n_one,at,xyz2,e,g_one)
   call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,&
        &sr42,c6xy,e,g_one)
   call ff_hb(n_one,at,xyz2,e,g_one)
   call ff_eg_two(n_one,at,xyz2,e_two,g_two)
   call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,&
        &r094_mod,sr42,c6xy_two,e_two,g_two)
   call ff_hb_two(n_one,at,xyz2,e_two,g_two)

   e_qmdff1=e+E_zero1
   e_qmdff2=e_two+E_zero2
   if (.not. use_rpmd) write(99,*) e_qmdff1,e_qmdff2
else 
   if (.not. use_rpmd) write(99,*) E_zero1,E_zero2
end if
call xyz_2int(xyz2,int_path,natoms)

!
!     If we are outside the region where the reference path is defined, set
!     the coupling to zero! In practice, all far distant points will have 
!     s values of 0.9999... oder 0.0000.., therefore, we set these also to zero
!     The IQ_R_HI and IQ_L_LO values are only different from 1 and 0 if the 
!     NO_EVB option is activated
!
if ((s_act .gt. iq_r_hi) .or. (s_act .lt. iq_l_lo)) then
   V12=0.d0
   zeta=0.d0
else 
   V12=0.d0 
   z_act=0.d0
   do i=1,nat6
      z_act=z_act+(int_path(i)-int_min(i))**2
   end do
   z_act=sqrt(z_act)
!
!    Calculate the projection vector of the actual structure with respect 
!    to the nearest structure on the IRC
!
   call int_diff(int_path,int_min,deltq)
!
!    Calculate the difference of both QMDFF energies for the analyctial gradient!
!
   if (calc_evb) ediff=e_qmdff1-e_qmdff2
        
!
!    Call the spline routine to estimate the coupling strength on the IRC
!
   call interp_spline(nat6+2,s_act,path_s,klo,khi)
!
!    Constant term of the upfollowing taylor expansion
!
   if (calc_evb) V12=V12+path_s(nat6+1)
!
!    Test: energy of reference method, interpolated
!
      
   if (calc_ens) Vref=path_s(nat6+2)
!
!    Store this value of analytical gradient
!
   if (calc_evb) V12_const=V12
!
!    If the structure is not located on the actual minimum path; calculate 
!    the gradient and hessian values for the taylor expansion
!    first, calculate in which spline region for gradients/hessians we are
!     
   do i=1,rp_evb_points
      if (s_act .ge. (rp_point_s(i)) .and. (s_act .le. rp_point_s(i+1))) then
         igh_left=i
         igh_right=i+1
      end if
   end do
!
!    Test: linear interpolation
!
   wgh_left=(rp_point_s(igh_right)-s_act)/(rp_point_s(igh_right)-&
         & rp_point_s(igh_left))
   wgh_right=1.d0-(rp_point_s(igh_right)-s_act)/(rp_point_s(igh_right)-&
         & rp_point_s(igh_left))
!     
!    For the usage of EVB-QMDFF coupling gradients and hessians
!
   if (calc_evb) then
       
      dv12_act=dv12(igh_left,:)*wgh_left+dv12(igh_right,:)*wgh_right
      d2v12_act=d2v12(igh_left,:,:)*wgh_left+d2v12(igh_right,:,:)*wgh_right
!
!    Now calculate the taylor expansion to the actual structure originating 
!    at the found nearest structure on the IRC
!
!
!    the linear term
!
      V12=V12+dot_product(dv12_act,deltq)
!
!    the quadratic term
!
      dq_dv12_act=matmul(d2v12_act,deltq)
      V12=V12+0.5d0*dot_product(deltq,dq_dv12_act) 
!
!    The exponential damping term
!
!    calculate the vector norm of the difference vector only once 
!    (also needed for analytical gradient)
!
      q_square=dot_product(deltq,deltq)
      V12=V12*exp(-q_square*pre_exp)
!
!    Finally, correct the coupling strengh at the areas near s=0 and s=1, to
!    avoid jumps there
!          
      if (s_act .le. pi/pareta) then
         zeta=0.5d0*(sin(s_act*pareta-pi/2.d0)+1.d0)
      else if (s_act .ge. 1-pi/pareta) then
         zeta=0.5d0*(sin(s_act*pareta-pareta-pi/2.d0)+1.d0)
      else
         zeta=1.d0
      end if 
      V12=V12*zeta  
   end if
!
!    Taylor expansion in the near-TS region for the reference energies 
!       
  if (calc_ens) then
     
     grad_act=all_grad(:,igh_left)*wgh_left+all_grad(:,igh_right)*wgh_right
     hess_act=all_hess(:,:,igh_left)*wgh_left+all_hess(:,:,igh_right)*wgh_right
!
!     the linear term
!         
     Vref=Vref+dot_product(grad_act,deltq)
!
!     the quadratic term
!            
      dq_hess_act=matmul(hess_act,deltq)
      Vref=Vref+0.5d0*dot_product(deltq,dq_hess_act)    
   end if
end if  
!
!   Calculate the RP-EVB-QMDFF energy
! 
!
!    Write internal coordinates to file if option is activated
!    Print a third column for splot an convert them into Angstrom
!
if (int_coord_plot) then
   write(191,*) int_path*bohr,0.0d0 
end if

if (calc_evb) then
   root=(0.5d0*(e_qmdff1-e_qmdff2))*(0.5d0*(e_qmdff1-e_qmdff2))+V12
   if (root .le. 0) then
      e_evb=0.5d0*(e_qmdff1+e_qmdff2)
   else
      e_evb=0.5d0*(e_qmdff1+e_qmdff2)-sqrt(root)
   end if
end if
!
!     For the mixing region, where reference energy and EVB-QMDFF are calculated,
!     apply the smoothstep function!
!
if (rp_transition) then
   if (s_act .lt. s_ts) then
      x_val=(s_act-trans_l_lo)/(trans_l_hi-trans_l_lo)

   else 
      x_val=1.d0-(s_act-trans_r_lo)/(trans_r_hi-trans_r_lo)
   end if
!
!     Calculate value of simple smoothstep function (3x^2-2x^3)
!
!         smoothstep=3.d0*x_val*x_val-2*x_val*x_val*x_val
!         smoothstep=10d0*x_val**3-15d0*x_val**4+6d0*x_val**5
!          smoothstep=-252d0*x_val**11+1386d0*x_val**10-3080d0*x_val**9+3465d0*x_val**8-1980d0*x_val**7+462d0*x_val**6
   smoothstep=x_val
!
!     Calculate final energy for this structure
!
!
   en_evb=e_evb
   en_ref=Vref
   e_evb=smoothstep*en_ref+(1-smoothstep)*en_evb
else if (.not. calc_evb .and. calc_ens) then
    e_evb=Vref
end if
 
!
!     If the NO_EVB option was activated, check if we are in the transition QMDFF/direct interpolation,
!     else, only use direct interpolation
!     For deactivated NO_EVB option, only set the direct interpolation
!
if (no_evb) then
   if ((s_act .le. iq_l_hi) .or. (s_act .ge. iq_r_lo)) then
!
!     calculate the QMDFF energies, no coupling
!
   call ff_eg(n_one,at,xyz2,e,g_one)
   call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,&
        &sr42,c6xy,e,g_one)
   call ff_hb(n_one,at,xyz2,e,g_one)
   call ff_eg_two(n_one,at,xyz2,e_two,g_two)
   call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,&
        &r094_mod,sr42,c6xy_two,e_two,g_two)
   call ff_hb_two(n_one,at,xyz2,e_two,g_two)

   e_qmdff1=e+E_zero1
   e_qmdff2=e_two+E_zero2

   root=(0.5d0*(e_qmdff1-e_qmdff2))*(0.5d0*(e_qmdff1-e_qmdff2))
   en_evb=0.5d0*(e_qmdff1+e_qmdff2)-sqrt(root)
!
!     calculate the interpolation weights between QMDFF and direct interpolation (simple 
!     linear interpolation will be used in the first try)
! 
!     If we are in the outside pure QMDFF region, set the weight for direct interpolation to zero
!
   smoothstep=0.d0
   x_val=0.d0
!
!     This code line is crucial!! Without this command, NaN errors occur! (nobody understands this)
!
            write(adum,*) 66.0!s_ts
   if (s_act .lt. s_ts) then
      if (s_act .ge. iq_l_lo) then
         x_val=(s_act-iq_l_lo)/(iq_l_hi-iq_l_lo)
      else 
         x_val=0.d0
      end if
   else
      if (s_act .le. iq_r_hi) then
         x_val=1.d0-(s_act-iq_r_lo)/(iq_r_hi-iq_r_lo)
      else 
         x_val=0.d0
      end if
   end if 
!   TEST switch off  
   if (xi_test .le. no_evb_xi) then
   smoothstep=0d0
   else 
   smoothstep=x_val
   end if
  !          smoothstep=3.d0*x_val*x_val-2*x_val*x_val*x_val
      e_evb=smoothstep*Vref+(1-smoothstep)*en_evb
   else 
      e_evb=Vref
   end if
!
!     apply linear interpolation 
!
end if
!
!     if the gen_test option was activated: increment the total energy 
!
if (gen_test) then
   gen_energies=gen_energies+e_evb 
end if
!
!   Now, calculate the numerical gradient!
!
if (num_grad .eqv. .true.) then
!
!   The numerical RP-EVB gradient
!
   write(*,*) "For RP-EVB, NUM_GRAD is currently not availiable!"
   call fatal
!
!    Do the analytical RP-EVB gradient
!    We need to calculate the gradient of the coupling dV12/dq!
!    This gradient will be done in three parts, regarding to the three 
!    coefficients of the taylor expansion of the coupling term
!
else 
!
!    The coupling gradient (usual used, may be replaced in TS region (see below)
!
   if (calc_evb) then
!
!    Calculate the difference of both QMDFF energies for the analyctial gradient!
!
      ediff=e_qmdff1-e_qmdff2
!            write(*,*) "ediff",ediff
!
!    set coupling to zero if we are (almost) outside the defined area
!
      if ((s_act .gt. iq_r_hi) .or. (s_act .lt. iq_l_lo)) then
         V12=0.d0
         g_V12=0.d0
      else 

!
!    First, the derivative with respect to the IRC (constant term)
!    Calculate the gradient vector of the interpolated curve, regarding to
!    the spline functions! Multiply derivative of the coupling and of the 
!    coordinate!         
!
!    Analytical derivatives of spline curve: tangential vector of coordinates, relative
!    to total gradient in energy
!

         h=rp_spl_s(khi)-rp_spl_s(klo)
         apar=(rp_spl_s(khi)-s_act)/h
         b=(s_act-rp_spl_s(klo))/h
         j=nat6+1

         dvds=-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+1.d0/6.d0*h*(1.d0-3.d0*apar*apar)*&
              & rp_spl_d2(klo,j)+1.d0/6.d0*h*(-1.d0+3.d0*b*b)*rp_spl_d2(khi,j)
!
!     Derivatives of all coordinates on the curve
!
         do j=1,nat6
            dqds(j)=-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+(1.d0/6.d0)*(h*h*(rp_spl_d2(khi,j)*&
                & (-3.d0*apar*apar/h+1.d0/h)+rp_spl_d2(klo,j)*(3.d0*b*b/h-1.d0/h)))
            dqds(j)=(-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+1.d0/6.d0*h*(1.d0-3.d0*apar*apar)*&
                & rp_spl_d2(klo,j)+1.d0/6.d0*h*(-1.d0+3.d0*b*b)*rp_spl_d2(khi,j))
         end do
         abs_dqds=dot_product(dqds,dqds)
         g_V12_int=dvds*dqds/(sqrt(abs_dqds))/s_tot
!
!     Now, calculate the derivative of the orthognal part, which is a taylor 
!     expansion with exponential damping around it
!     Sum up all terms!
!
!     First
!
         dq_taylor=pre_exp*q_square*matmul(d2v12_act,deltq)
! 
!     Second
!
         dq_taylor=dq_taylor+2.d0*pre_exp*q_square*dv12_act
!
!     Third
!
         dq_taylor=dq_taylor+2.d0*V12_const*pre_exp*deltq
!
!     Fourth
!
         dq_taylor=dq_taylor-matmul(deltq,d2v12_act)
!
!     Fifth
!
         dq_taylor=dq_taylor-dv12_act
!
!     Add minus sign and exponential damping
!
         dq_taylor=-dq_taylor*exp(-pre_exp*q_square)
!
!     Subtract the projected part of the taylor gradient on the path from the full
!     gradient to obtain the orthogonal part
!
         taylor_proj=(dq_taylor-(dot_product(dq_taylor,dqds)/abs_dqds)*dqds)
         
         g_V12_int=g_V12_int+taylor_proj   
!
!     Convert the squared coupling gradient to cartesian coordinates
!
         call int2grad(xyz2,int_path,g_V12_int,g_V12)

!
!     Scale the coupling gradient with the zeta parameter to archive a smooth transition
!     to the uncoupled regions and to avoid artificial asymptotic behavior of the 
!     spline functions
!
         g_V12=zeta*g_V12
      end if
!
!     Now, calculate the EVB-QMDFF gradient using the well known formula
!     Using some mathematical tricks, we can avoid the root of the coupling term!        
! 
      offdiag=V12
      off4=4.d0*offdiag
  
!
!     If the root argument becomes negative, set the coupling to zero..
! 
           
      unset_V12=.false.
      if (ediff*ediff+off4 .lt. 0.d0) unset_V12=.true.
      if (unset_V12 .eqv..false.) then
         root2=sqrt(ediff*ediff+off4)
      end if
      do i=1,n_one
         do j=1,3
            deldiscr=ediff*(g_one(j,i)-g_two(j,i))+2.d0*g_V12(j,i)
            if (unset_V12 .eqv. .false.) then
               delsqrt=deldiscr/root2 
            else 
               delsqrt=0.d0
            end if
            g_evb(j,i)=0.5d0*(g_one(j,i)+g_two(j,i)-delsqrt)
         end do
      end do   
   end if
!
!     For the direct interpolation of reference energies/gradients etc near the TS
!
   if (calc_ens) then
!
!     Forst, derivatives with respect to the IRC (constant term). Calculate the 
!     potential energy gradient vector of the interpolated curve, regarding to the 
!     spline functions!
!
!
!    Analytical derivatives of spline curve: tangential vector of coordinates, relative
!    to total gradient in energy
!
      h=rp_spl_s(khi)-rp_spl_s(klo)
      apar=(rp_spl_s(khi)-s_act)/h
      b=(s_act-rp_spl_s(klo))/h
      j=nat6+2
 
      dvds=-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+1.d0/6.d0*h*(1.d0-3.d0*apar*apar)*&
           & rp_spl_d2(klo,j)+1.d0/6.d0*h*(-1.d0+3.d0*b*b)*rp_spl_d2(khi,j)
!
!     Derivatives of all coordinates on the curve
!
      do j=1,nat6
          dqds(j)=-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+(1.d0/6.d0)*(h*h*(rp_spl_d2(khi,j)*&
              & (-3.d0*apar*apar/h+1.d0/h)+rp_spl_d2(klo,j)*(3.d0*b*b/h-1.d0/h)))
         dqds(j)=(-rp_spl_ref(j,klo)/h+rp_spl_ref(j,khi)/h+1.d0/6.d0*h*(1.d0-3.d0*apar*apar)*&
              & rp_spl_d2(klo,j)+1.d0/6.d0*h*(-1.d0+3.d0*b*b)*rp_spl_d2(khi,j))
      end do
      abs_dqds=dot_product(dqds,dqds)
      grad_int=dvds*dqds/(sqrt(abs_dqds))/s_tot!/2
!
!     Now, calculate the derivative of the orthogonal part, which is a taylor 
!     expansion of the energy
!
!     1. derivative of linear term           
!
      grad_taylor=grad_act
!
!     2. derivative of quadratic term
!
      grad_taylor=grad_taylor+matmul(hess_act,deltq)
      
      taylor_proj=(grad_taylor-(dot_product(grad_taylor,dqds)/abs_dqds)*dqds)
      grad_int=grad_int+taylor_proj    
!
!     convert the gradient to cartesian cordinates!
!

      call int2grad(xyz2,int_path,grad_int,grad_v_xyz)
   end if     
!
!     Calculate final gradient for this structure: same procedure as above for energies..
!
   if (rp_transition) then
      g_evb=smoothstep*grad_v_xyz+(1-smoothstep)*g_evb
   else if (.not. calc_evb .and. calc_ens) then
      g_evb=grad_v_xyz
   end if
!
!     If the NO_EVB option was activated, apply the same scheme as for energies 
!
!
   if (no_evb) then
      if ((s_act .le. iq_l_hi) .or. (s_act .ge. iq_r_lo)) then
!
!     calculate the QMDFF gradient, no coupling
!
         root2=sqrt(ediff*ediff)
         do i=1,n_one
            do j=1,3
               deldiscr=ediff*(g_one(j,i)-g_two(j,i))
                  delsqrt=deldiscr/root2
               g_evb(j,i)=0.5d0*(g_one(j,i)+g_two(j,i)-delsqrt)
            end do
         end do
!
!     calculate the interpolation weights between QMDFF and direct interpolation (simple 
!     linear interpolation will be used in the first try)
! 
!     If we are in the outside pure QMDFF region, set the weight for direct interpolation to zero
!
         g_evb=smoothstep*grad_v_xyz+(1-smoothstep)*g_evb
      end if
   end if

!
!     calculate gradient norm
!
   gnorm=sqrt(sum(g_evb**2)) 
!
!     write out energy, s-value, z-value and structure for testing during generation
!
   if (gen_test) then
      if (s_act .gt. 0.00001d0 .and. s_act .lt. 0.9999d0) then
      gen_pr_act=gen_pr_act+1
      if (mod(gen_pr_act,gen_pr_frac) .eq. 0) then 
!
!     write energy and s/z values
!
      write(296,*) e_evb,s_act,z_act
!
!     write structure of first bead
! 
         write(295,*) natoms
         write(295,*) 
         do i=1,natoms
            write(295,*) name(i),q_i(:,i,1)*bohr
         end do
      end if
      end if
   end if 
! ######################################################################
!      DEBUG SECTION: COMMENT IN TO APPLY IT
!      ATTENTION: LARGE OUTPUT WILL BE WRITTEN!
!
!           if ((s_act .le. iq_l_hi) .or. (s_act .ge. iq_r_lo)) then
!             write(269,*) s_act,z_act,e_evb*2600,smoothstep!+1230000
!           else 
!             write(269,*) s_act,z_act,e_evb*2600
!           end if
!          if (e_evb*2600 .gt. -840000) then
!             write(189,*) "ERROR-REPORT:"
!             do j=1,nat6
!                write(189,*) coord_def(j,:),deltq(j)
!             end do
!             write(190,*) deltq(:)
!          end if
!
!      END DEBUG SECTION
! ######################################################################
end if
!
!     Optional test: write components of internal gradients to file
!
if (int_grad_plot) then
!   g_int=0.d0
!   call grad2int(xyz2,int_path,g_int,g_evb)
!   write(192,*) -g_int!1000 !-g_V12_int
end if   
   return
   end if
!
!   The dQ-couplingterms (2x2-EVB)
!
   if (use_dq .eqv. .true.) then 
      delta_Q=0
      allocate(int_ts2(3*natoms-6))
      allocate(int_coord(3*natoms-6))
      xyz2=xyz2*bohr
      call xyz_2int(xyz2,int_coord,nat)
      int_ts2=ts_coordinates
      do i=1,3*nat-6
         dq(i)=int_coord(i)-int_ts2(i)
         if (i .gt. 2*nat-3) then
            dq(i)=sqrt((int_coord(i)-int_ts2(i))*(int_coord(i)-&
                 &int_ts2(i))+(int_coord(i+nat-3)-int_ts2(i+nat-3))*&
                 &(int_coord(i+nat-3)-int_ts2(i+nat-3)))
         end if
      end do
      xyz2=xyz2/bohr

!   modded, if there are new problems, look here
      do i=1,3*natoms-6
         delta_Q=delta_Q+(dq(i)*dq(i))
      end do
      delta_Q=dsqrt(delta_Q)/10
      if (num_grad .eqv. .true.) then
         call gradnum(xyz2,g_evb,delta_Q,delta_Q2)
      else
         if (off_basis=="const") then
            offdiag=offa      !a
         else if (off_basis=="1g") then
            offdiag=offa*exp(-offb*delta_Q*delta_Q)  !a*exp(-b*deltaE^2)
         else if (off_basis=="3g") then
            offdiag=offa*exp(-offb*delta_Q*delta_Q)+offc*exp(-offd*(delta_Q+offm)&
                 &*(delta_Q+offm))+offe*exp(-offf*(delta_Q+offn)*(delta_Q+offn))
         else if (off_basis=="sd2") then
            offdiag=offa*exp(-offb*delta_Q*delta_Q)+offc*delta_Q*delta_Q*&
                 &exp(-offd*delta_Q*delta_Q)+offe*delta_Q*delta_Q*&
                 &exp(-offf*delta_Q*delta_Q)
         end if
!
!     Now calculate the analytical evb-qmdff-gradient!
!     in the dQ-case, the derivatives of the off-diagonal elements
!     (dV12/dQ) vanishes, because they have no dependance on energy!
!

   
         off4=4.d0*offdiag*offdiag
         root2=sqrt(ediff*ediff+off4)
         e_evb=0.5d0*(e1_shifted+e2_shifted-root2)
         do i=1,n_one
            do j=1,3
               if (off_basis=="const") then
                  deldiscr=ediff*(g_one(j,i)-g_two(j,i))
               else if (off_basis=="1g") then
                  deldiscr=ediff*(g_one(j,i)-g_two(j,i))
               else if (off_basis=="3g") then
                  deldiscr=ediff*(g_one(j,i)-g_two(j,i))
               else if (off_basis=="sd2") then
                  deldiscr=ediff*(g_one(j,i)-g_two(j,i))
               end if
               delsqrt=deldiscr/root2
               g_evb(j,i)=0.5d0*(g_one(j,i)+g_two(j,i)-delsqrt)
!
!    Convert from hartree/bohr to kcal/mol*angstöm; not for EVB-QMDFF!
!
               g_evb(j,i)=0!g_evb(j,i)!*hartree*bohr
            end do
         end do
      end if
!  
!     diagonalize 2x2 matrix and take the lower of the two eigenvalues;
!     we do this by simply evaluating the ubiquitous analytic expression 
!     for the lower eigenvalue of a symmetric 2x2 matrix ("Zweizustandsproblem"):
!     All the different implemented coupling-terms are recognized
!

!
!     for the dE-couplingterms
!
   else
      if (off_basis=="const") then
         offdiag=offa      !a
      else if (off_basis=="1g") then
         offdiag=offa*exp(-offb*ediff*ediff)  !a*exp(-b*deltaE^2)
      else if (off_basis=="2g") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*exp(-offd*(ediff+offm)&
                 &*(ediff+offm))
      else if (off_basis=="3g") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*exp(-offd*(ediff+offm)&
                 &*(ediff+offm))+offe*exp(-offf*(ediff+offn)*(ediff+offn))
      else if (off_basis=="sp") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*exp(-offd*ediff*ediff)
      else if (off_basis=="sd") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*ediff*exp(-offd*ediff*ediff)
      else if (off_basis=="sd2") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*ediff*&
                 &exp(-offd*ediff*ediff)+offe*ediff*ediff*exp(-offf*ediff*ediff)
      else if (off_basis=="sp2d") then
         offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*&
                 &exp(-offd*ediff*ediff)+offe*ediff*exp(-offf*ediff*ediff)+&
                 &ediff*ediff*offg*exp(-offh*ediff*ediff)
      end if
!
!     Now calculate the analytical evb-qmdff-gradient!
!
      if (num_grad .eqv. .true.) then
          call gradnum(xyz2,g_evb,delta_Q,delta_Q2)
      else
      off4=4.d0*offdiag*offdiag
      root2=sqrt(ediff*ediff+off4)
      e_evb=0.5d0*(e1_shifted+e2_shifted-root2)
      do i=1,n_one
         do j=1,3
            if (off_basis=="const") then
               deldiscr=ediff*(g_one(j,i)-g_two(j,i))
            else if (off_basis=="1g") then
               deldiscr=ediff*(g_one(j,i)-g_two(j,i))*(1.d0-2.d0*off4*offb)
            else if (off_basis=="2g") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)-2d0*&
                  &offc*offd*(offm+ediff)*exp(-offd*(ediff+offm)*(ediff+offm))))  
            else if (off_basis=="3g") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)-2d0*&
                  &offc*offd*(offm+ediff)*exp(-offd*(ediff+offm)*(ediff+offm))&
                  &-2d0*offe*offf*(offn+ediff)*exp(-offf*(ediff+offn)*(ediff+offn))))  
            else if (off_basis=="sp") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)+offc*&
                  &(1d0-2d0*ediff*ediff*offd)*exp(-offd*ediff*ediff)))      
            else if (off_basis=="sd") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)+offc*&
                  &(2d0*ediff-2d0*ediff*ediff*ediff*offd)*exp(-offd*ediff*ediff)))
            else if (off_basis=="sd2") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)+offc*&
                  &(2d0*ediff-2d0*ediff*ediff*ediff*offd)*exp(-offd*ediff*ediff)+&
                  &offe*(2d0*ediff-2d0*ediff*ediff*ediff*offf)*exp(-offf*ediff*ediff)))
            else if (off_basis=="sp2d") then
               deldiscr=(g_one(j,i)-g_two(j,i))*(1d0*ediff+4d0*offdiag*&
                  &(-2d0*offa*offb*ediff*exp(-offb*ediff*ediff)+offc*&
                  &(1d0-2d0*ediff*ediff*offd)*exp(-offd*ediff*ediff)+&
                  &offe*(1d0-2d0*ediff*ediff*offf)*exp(-offf*ediff*ediff)+&
                  &offg*(2d0*ediff-2d0*ediff*ediff*ediff*offh)*exp(-offh*ediff*ediff)))
            end if
            delsqrt=deldiscr/root2
            g_evb(j,i)=0.5d0*(g_one(j,i)+g_two(j,i)-delsqrt)
!
!    Convert from hartree/bohr to kcal/mol*angstöm
!    ---> not in the case of the evb_qmdff program!!
!
            g_evb(j,i)=g_evb(j,i)!*hartree*bohr
         end do
      end do
   end if
   end if
!
!    In the case of 3 QMDFF´s, no analytical results are availiable:
!    instead matrix-calculations are needed to determine energy and
!    gradients
!
else if (nqmdff.eq.3) then
!    The second QMDFF
   call ff_eg_two(n_one,at,xyz2,e_two,g_two)
   call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,r094_mod,sr42,&
                &c6xy_two,e_two,g_two)
   call ff_hb_two(n_one,at,xyz2,e_two,g_two)

!    The third QMDFF
   call ff_eg_three(n_one,at,xyz2,e_three,g_three)
   call ff_nonb_three(n_one,at,xyz2,q_three,r0ab,zab,r094_mod,sr42,&
                &c6xy_three,e_three,g_three)
   call ff_hb_three(n_one,at,xyz2,e_three,g_three)

   gnorm=sqrt(sum(g_one**2))
   gnorm_two=sqrt(sum(g_two**2))
   gnorm_three=sqrt(sum(g_three**2))

!
!   Calculate shifted energies of the QMDFFs (corrected with ab-initio-energies)
!

   e1_shifted=e+E_zero1
   e2_shifted=e_two+E_zero2
   e3_shifted=e_three+E_zero3
   e_qmdff1=e1_shifted
   e_qmdff2=e2_shifted
   e_qmdff3=e3_shifted
!
!     for the dQ-couplingterms
!
   if (use_dq .eqv. .true.) then
      delta_Q=0
      delta_Q2=0
      xyz2=xyz2*bohr
      call xyz_2int(xyz2,int_coord,nat)
      xyz2=xyz2/bohr
!     The geometrical coordinate dq=q-q_ts in internal coordinates
      int_ts2=ts_coordinates
      do i=1,3*nat-6
         dq(i)=int_coord(i)-int_ts2(i)
         if (i .gt. 2*nat-3) then
            dq(i)=sqrt((int_coord(i)-int_ts2(i))*(int_coord(i)-int_ts2(i))+&
                  &(int_coord(i+nat-3)-int_ts2(i+nat-3))*&
                  &(int_coord(i+nat-3)-int_ts2(i+nat-3)))
         end if
      end do
      int_ts2=ts_coordinates2
      do i=1,3*nat-6
         dq2(i)=int_coord(i)-int_ts2(i)
         if (i .gt. 2*nat-3) then
            dq2(i)=sqrt((int_coord(i)-int_ts2(i))*(int_coord(i)-int_ts2(i))+&
                  &(int_coord(i+nat-3)-int_ts2(i+nat-3))*&
                  &(int_coord(i+nat-3)-int_ts2(i+nat-3)))
         end if
      end do

      do i=1,3*natoms
         delta_Q=delta_Q+(dq(i)*dq(i))
         delta_Q2=delta_Q2+(dq2(i)*dq2(i))
      end do
      delta_Q=dsqrt(delta_Q)/50
      delta_Q2=dsqrt(delta_Q2)/50
!
!     if desired,calculate the numerical gradient
!
      JOBZ='V' !eigenvalues and eigenvectors(U)
      UPLO='U' !upper triangle of a
      Nd=3
      LDA=Nd
      INFO=0
      LWORK=Nd*Nd-1
      allocate(A(Nd,Nd),U(Nd,Nd),Ut(Nd,Nd),UtA(nD,Nd))
      allocate(W(Nd))
      allocate(WORK(LWORK))
      if (off_basis=="1g") then
         evb_mat(1,1)=e_qmdff1
         evb_mat(2,2)=e_qmdff2
         evb_mat(3,3)=e_qmdff3
         evb_mat(1,2)=offa*exp(-offb*delta_Q*delta_Q)
         evb_mat(2,1)=evb_mat(1,2)
         evb_mat(2,3)=offc*exp(-offd*delta_Q2*delta_Q2)
         evb_mat(3,2)=evb_mat(2,3)
         evb_mat(1,3)=0
         evb_mat(3,1)=0
         A=evb_mat
         call DSYEV(JOBZ,UPLO,Nd,A,LDA,W,WORK,LWORK,INFO)
         E_evb=W(1)
      else if (off_basis=="sd2") then
         evb_mat(1,1)=e_qmdff1
         evb_mat(2,2)=e_qmdff2
         evb_mat(3,3)=e_qmdff3
         evb_mat(1,2)=offa*exp(-offb*delta_Q*delta_Q)+offc*delta_Q*delta_Q*&
                &exp(-offd*delta_Q*delta_Q)+offe*delta_Q*delta_Q*exp(-offf*delta_Q*delta_Q)
         evb_mat(2,1)=evb_mat(1,2)
         evb_mat(2,3)=offg*exp(-offh*delta_Q2*delta_Q2)+offi*delta_Q2*delta_Q2*&
                &exp(-offj*delta_Q2*delta_Q2)+offk*delta_Q2*delta_Q2*exp(-offl*&
                &delta_Q2*delta_Q2)
         evb_mat(3,2)=evb_mat(2,3)
         evb_mat(1,3)=0
         evb_mat(3,1)=0
         A=evb_mat
         call DSYEV(JOBZ,UPLO,Nd,A,LDA,W,WORK,LWORK,INFO)
         E_evb=W(1)
      end if
!
!     Now calculate the analytical evb-qmdff-gradient!
!     This is done by applying the Hellman-Feynman-theoreme
!     To every component of the gradient-vector, the 
!     formula dD/dq=U^T*dH/dq*U is evaluated
!
      U=A
      do i=1,Nd
         do j=1,Nd
            Ut(i,j)=U(j,i)
         end do
      end do
!
!     if desired,calculate the numerical gradient
!
      if (num_grad .eqv. .true.) then
          call gradnum(xyz2,g_evb,delta_Q,delta_Q2)
      else
         do i=1,n_one
            do j=1,3
               evb_mat(1,1)=g_one(j,i)
               evb_mat(2,2)=g_two(j,i)
               evb_mat(3,3)=g_three(j,i)
               evb_mat(1,2)=0
               evb_mat(2,1)=0
               evb_mat(2,3)=0
               evb_mat(3,2)=0
               evb_mat(1,3)=0
               evb_mat(3,1)=0
               UtA=matmul(Ut,evb_mat)
               g_mat=matmul(UtA,U)
               g_evb(j,i)=g_mat(1,1)
               g_evb(j,i)=g_evb(j,i)*hartree*bohr
            end do
         end do
         gnorm_evb=sqrt(sum(g_evb**2))
      end if
   else 
!
!     for the dE-couplingterms
!
         JOBZ='V' !eigenvalues and eigenvectors(U)
         UPLO='U' !upper triangle of a
         Nd=3
      LDA=Nd
      INFO=0
      LWORK=Nd*Nd-1
      allocate(A(Nd,Nd),U(Nd,Nd),Ut(Nd,Nd),UtA(nD,Nd))
      allocate(W(Nd))
      allocate(WORK(LWORK))
      if (off_basis=="const") then
         write(*,*) "Sorry, this couplingterm isn´t implemented yet!"
         call fatal
      else if (off_basis=="1g") then
         dE12=abs(e_qmdff1-e_qmdff2)
         dE13=abs(e_qmdff1-e_qmdff3)
         dE23=abs(e_qmdff2-e_qmdff3)
         evb_mat(1,1)=e_qmdff1
         evb_mat(2,2)=e_qmdff2
         evb_mat(3,3)=e_qmdff3
         evb_mat(1,2)=offa*exp(-offb*dE12*dE12)
         evb_mat(2,1)=evb_mat(1,2)
         evb_mat(2,3)=offc*exp(-offd*dE23*dE23)
         evb_mat(3,2)=evb_mat(2,3)
         if (full .eqv. .true.) then
            evb_mat(1,3)=offe*exp(-offf*dE13*dE13)
            evb_mat(3,1)=evb_mat(1,3)
         else 
         evb_mat(1,3)=0
         evb_mat(3,1)=0
         end if
         A=evb_mat
         call DSYEV(JOBZ,UPLO,Nd,A,LDA,W,WORK,LWORK,INFO)
         E_evb=W(1) 
!
!     Now calculate the analytical evb-qmdff-gradient!
!     This is done by applying the Hellman-Feynman-theoreme
!     To every component of the gradient-vector, the 
!     formula dD/dq=U^T*dH/dq*U is evaluated
!
         U=A
         do i=1,Nd
            do j=1,Nd
               Ut(i,j)=U(j,i)
            end do
         end do
!
!     if desired,calculate the numerical gradient
!
         if (num_grad .eqv. .true.) then
             call gradnum(xyz2,g_evb,delta_Q,delta_Q2)
         else
             do i=1,n_one
                do j=1,3
                   evb_mat(1,1)=g_one(j,i)
                   evb_mat(2,2)=g_two(j,i)
                   evb_mat(3,3)=g_three(j,i)
                   evb_mat(1,2)=-2d0*offa*offb*(e_qmdff1-e_qmdff2)*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offb*dE12*dE12)
                   evb_mat(2,1)=evb_mat(1,2)
                   evb_mat(2,3)=-2d0*offc*offd*(e_qmdff2-e_qmdff3)*(g_two(j,i)-g_three(j,i))*&
                                &exp(-offd*dE23*dE23)
                   evb_mat(3,2)=evb_mat(2,3)
                   if (full .eqv. .true.) then
                      evb_mat(1,3)=-2d0*offe*offf*(e_qmdff1-e_qmdff3)*(g_one(j,i)-&
                                  &g_three(j,i))*exp(-offf*dE13*dE13)
                      evb_mat(3,1)=evb_mat(1,3)
                   else
                      evb_mat(1,3)=0
                      evb_mat(3,1)=0
                   end if
                   UtA=matmul(Ut,evb_mat)
                   g_mat=matmul(UtA,U)
                   g_evb(j,i)=g_mat(1,1)
                   g_evb(j,i)=g_evb(j,i)*hartree*bohr
                end do
             end do
             gnorm_evb=sqrt(sum(g_evb**2))
         end if
      else if (off_basis=="sd2") then
         dE12=abs(e_qmdff1-e_qmdff2)
         dE13=abs(e_qmdff1-e_qmdff3)
         dE23=abs(e_qmdff2-e_qmdff3)
         evb_mat(1,1)=e_qmdff1
         evb_mat(2,2)=e_qmdff2
         evb_mat(3,3)=e_qmdff3
         evb_mat(1,2)=offa*exp(-offb*dE12*dE12)+offc*dE12*dE12*&
                     &exp(-offd*dE12*dE12)+offe*dE12*dE12*exp(-offf*dE12*dE12)
         evb_mat(2,1)=evb_mat(1,2)
         evb_mat(2,3)=offg*exp(-offh*dE23*dE23)+offi*dE23*dE23*&
                     &exp(-offj*dE12*dE23)+offk*dE23*dE23*exp(-offl*dE23*dE23)
         evb_mat(3,2)=evb_mat(2,3)
         if (full .eqv. .true.) then
            evb_mat(1,3)=offm*exp(-offn*dE13*dE13)
            evb_mat(3,1)=evb_mat(1,3)
         else 
            evb_mat(1,3)=0
            evb_mat(3,1)=0
         end if
         A=evb_mat
         call DSYEV(JOBZ,UPLO,Nd,A,LDA,W,WORK,LWORK,INFO)
         E_evb=W(1) 
!
!     Now calculate the analytical evb-qmdff-gradient!
!     This is done by applying the Hellman-Feynman-theoreme
!     To every component of the gradient-vector, the 
!     formula dD/dq=U^T*dH/dq*U is evaluated
!
         U=A
         do i=1,Nd
            do j=1,Nd
               Ut(i,j)=U(j,i)
            end do
         end do
!
!     if desired,calculate the numerical gradient
!
         if (num_grad .eqv. .true.) then
             call gradnum(xyz2,g_evb,delta_Q,delta_Q2)
         else
             do i=1,n_one
                do j=1,3
                   evb_mat(1,1)=g_one(j,i)
                   evb_mat(2,2)=g_two(j,i)
                   evb_mat(3,3)=g_three(j,i)
                   evb_mat(1,2)=-2d0*offa*offb*(e_qmdff1-e_qmdff2)*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offb*dE12*dE12)+2d0*offc*(e_qmdff1-e_qmdff2)*&
                                &exp(-offd*(dE12*dE12))*(g_one(j,i)-g_two(j,i))-&
                                &2d0*offc*(e_qmdff1-e_qmdff2)**3*offd*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offd*(dE12*dE12))+2d0*offe*(e_qmdff1-e_qmdff2)*&
                                &exp(-offf*(dE12*dE12))*(g_one(j,i)-g_two(j,i))-&
                                &2d0*offe*(e_qmdff1-e_qmdff2)**3*offf*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offf*(dE12*dE12))
                   evb_mat(2,1)=evb_mat(1,2)
                   evb_mat(2,3)=-2d0*offg*offh*(e_qmdff1-e_qmdff2)*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offh*dE12*dE12)+2d0*offi*(e_qmdff1-e_qmdff2)*&
                                &exp(-offj*(dE12*dE12))*(g_one(j,i)-g_two(j,i))-&
                                &2d0*offi*(e_qmdff1-e_qmdff2)**3*offj*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offj*(dE12*dE12))+2d0*offk*(e_qmdff1-e_qmdff2)*&
                                &exp(-offl*(dE12*dE12))*(g_one(j,i)-g_two(j,i))-&
                                &2d0*offk*(e_qmdff1-e_qmdff2)**3*offl*(g_one(j,i)-g_two(j,i))*&
                                &exp(-offl*(dE12*dE12))
                   evb_mat(3,2)=evb_mat(2,3)
                   if (full .eqv. .true.) then
                      evb_mat(1,3)=-2d0*offm*offn*(e_qmdff1-e_qmdff3)*(g_one(j,i)-g_three(j,i))*&
                                  &exp(-offn*dE13*dE13)
                      evb_mat(3,1)=evb_mat(1,3)
                   else
                      evb_mat(1,3)=0
                      evb_mat(3,1)=0
                   end if
                   UtA=matmul(Ut,evb_mat)
                   g_mat=matmul(UtA,U)
                   g_evb(j,i)=g_mat(1,1)
                   g_evb(j,i)=g_evb(j,i)*hartree*bohr
                end do
             end do
             gnorm_evb=sqrt(sum(g_evb**2))
         end if
      else 
      write(*,*) "not implemented yet!"
      end if
      deallocate(A,W,WORK)
   end if
end if
return
end subroutine gradient
