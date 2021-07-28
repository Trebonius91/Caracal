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
!      subroutine egrad_treq: calculate energy and gradient of the 
!        transition region corrected(TRC)-reaction path (RP)-EVB-QMDFF (TREQ)
!
!      part of EVB
!
subroutine egrad_treq(xyz2,e_evb,g_evb,act_bead)
use general
use evb_mod 
implicit none
!      the actual structure  
real(kind=8), intent(inout) :: xyz2(3,natoms)
!      the calculated DG-EVB energy and gradient 
real(kind=8), intent(out) :: e_evb,g_evb(3,natoms)
!      the actual internal coordinate set
real(kind=8)::int_path(nat6)
!      loop indices 
integer::i,j,k
!      the actual topology parameters (path progress s and distance z)
real(kind=8)::s_act,z_act
!      switches for relevant path parts 
logical::calc_ens,calc_evb,rp_transition
!      general EVB calculation parameters 
real(kind=8)::offdiag,off4,root2,deldiscr,delsqrt,root
!      QMDFF energy parameters 
real(kind=8)::e,e_two,e_qmdff1,e_qmdff2,ediff
!      the EVB coupling strength
real(kind=8)::V12
!      the EVB coupling damping function at the path vicinities 
real(kind=8)::zeta
!      number of the actual RPMD bead 
integer::act_bead
!      the projected internal coordinates 
real(kind=8)::int_min(nat6)
!      difference vector between actual and projected structures 
real(kind=8)::deltq(nat6)
!      borders of different spline segments 
integer::klo,khi
!      energy of the reference method along the path
real(kind=8)::Vref
!      stored coupling for analytical gradient 
real(kind=8)::V12_const
!      spline segment borders for gradient and hessian interpolation
integer::igh_left,igh_right
!      weights of gradient/Hessian components interpolated
real(kind=8)::wgh_left,wgh_right
!      the interpolated coupling gradient 
real(kind=8)::dv12_act(nat6)
!      the interpolated coupling Hessian
real(kind=8)::d2v12_act(nat6,nat6)
!      intermediary matrix vector product 
real(kind=8)::dq_dv12_act(nat6)
!      squared difference vector deltq
real(kind=8)::q_square
!      the interpolated referene gradient 
real(kind=8)::grad_act(nat6)
!      the interpolated reference Hessian
real(kind=8)::hess_act(nat6,nat6)
!      intermediary matrix vector product 
real(kind=8)::dq_hess_act(nat6)
!      weighting factor for linear interpolation
real(kind=8)::x_val,smoothstep
!      interpolated EVB/reference energies 
real(kind=8)::en_evb,en_ref
!      parameters for gradient of spline functions 
real(kind=8)::h,apar,b
!      spline derivative of coupling 
real(kind=8)::dvds
!      spline derivatives of internal coordinates 
real(kind=8)::dqds(nat6)
!      absolute value of coordinate derivative 
real(kind=8)::abs_dqds
!      internal and cartesian coupling gradient 
real(kind=8)::g_V12_int(nat6),g_V12(3,natoms)
!      Taylor series part of coupling gradient 
real(kind=8)::dq_taylor(nat6)
!      projected part of Taylor expansion gradient 
real(kind=8)::taylor_proj(nat6)
!      Taylor series gradient for direct interpolation
real(kind=8)::grad_taylor(nat6)
!      the Pi
real(kind=8)::pi
parameter (pi=3.141592653589793238d0)
!      cartesian gradient of energy interpolation
real(kind=8)::grad_v_xyz(3,natoms)  
!      interpolated strcutures/energies along the path
real(kind=8),dimension(nat6+2)::path_s  
!      calculated norm of gradient vector 
real(kind=8)::gnorm
!      if the coupling shall be calculated 
logical::unset_V12
!      internal gradient vector 
real(kind=8)::grad_int(nat6)
!
character(len=80)::adum      
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
            write(adum,*) 66.0!s_t
   smoothstep=0.d0
   x_val=0.d0

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
   smoothstep=x_val
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
           if ((s_act .le. iq_l_hi) .or. (s_act .ge. iq_r_lo)) then
             write(269,*) s_act,z_act,e_evb*2600,smoothstep!+1230000
           else 
             write(269,*) s_act,z_act,e_evb*2600
           end if
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
end subroutine egrad_treq
