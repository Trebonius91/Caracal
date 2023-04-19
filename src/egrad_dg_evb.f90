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
!      subroutine egrad_dg_evb: calculate energy and gradient of the Distributed 
!        Gaussian (DG)-EVB coupling term 
!
!      part of EVB
!
subroutine egrad_dg_evb(xyz2,e_qmdff1,e_qmdff2,ediff,e_evb,g_evb)
use general
use evb_mod


implicit none
!      the actual structure  
real(kind=8), intent(inout) :: xyz2(3,natoms)
!      energies of both QMDFFs and their difference  
real(kind=8), intent(inout) :: e_qmdff1,e_qmdff2,ediff
!      the calculated DG-EVB energy and gradient 
real(kind=8), intent(out) :: e_evb,g_evb(3,natoms)
!      loop indices 
integer::i,j,k
!      the DG-EVB mode 
integer::dg_evb_mode
!      size of the DG-EVB coefficient matrix 
integer::mat_size 
!      the current internal coordinates 
real(kind=8),dimension(nat6)::int_path
!      general EVB calculation parameters 
real(kind=8)::offdiag,off4,root2,deldiscr,delsqrt,root 
!      the EVB coupling strength and its internal and cartesian derivatives
real(kind=8)::V12,g_V12_int(nat6),g_V12(3,natoms)
!      energies of the QMDFFs without shifting 
real(kind=8)::e,e_two
!      for numrical gradient 
real(kind=8)::e_lower,e_upper,step
!      the internal gradient 
real(kind=8)::g_int(nat6)
!      avoid imaginary coupling values 
logical::unset_V12
!   calculate the gradient (analytical and mumerical):
!

dg_evb_mode=dg_mode
if (dg_evb_mode .eq. 1) then
   mat_size=dg_evb_points
else if (dg_evb_mode .eq. 2) then
   mat_size=dg_evb_points*(1+nat6)
else if (dg_evb_mode .eq. 3) then
   mat_size=dg_evb_points*(1+nat6+(nat6)*(nat6+1)/2)
end if
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

            call xyz_2int(xyz2,int_path,natoms)

            e_qmdff1=e+E_zero1
            e_qmdff2=e_two+E_zero2

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

   call xyz_2int(xyz2,int_path,natoms)      
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

call xyz_2int(xyz2,int_path,natoms)
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

return
end subroutine egrad_dg_evb
