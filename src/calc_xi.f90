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
!     subroutine calc_xi: calculate value as well as first and 
!       second derivatives of the reaction coordinate Xi for 
!       the different types of implemented reactions
!     
!     Disclaimer (taken from with slightly changes):
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
!
!   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
!                         William H. Green (whgreen@mit.edu)
!                         Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
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
subroutine calc_xi(coords,xi_ideal,xi_act,dxi_act,d2xi_act,mode)
use general
use evb_mod
implicit none
integer::i,j,k,l,i1,i2,j1,j2  ! loop indices
integer::mode  ! if usual umbrella sampling or recrossing shall be calculated!
integer::s0_terms  ! number of interaction betweem sum_eds educt particles
real(kind=8)::coords(3,natoms) ! the actual coordinates
real(kind=8)::Red(3,sum_eds,sum_eds)  ! all possible educt-educt distances
real(kind=8)::r_eds(sum_eds,sum_eds) ! the actual distane between two educts
real(kind=8)::eds_avg  ! the average distance between two educts
real(kind=8)::R_f(3,form_num),R_b(3,break_num)  ! generalized vectors for all bonds
real(kind=8)::form_act(form_num),break_act(break_num)  ! actual bondlengths
real(kind=8)::xi_act  ! the calculated Xi value
real(kind=8)::dxi_act(3,natoms)  ! the calculated Xi derivative
real(kind=8)::d2xi_act(3,natoms,3,natoms)  ! the calculated Xi hessian 
real(kind=8)::s0,s1   ! dividing surface values 
real(kind=8)::ds0(3,natoms),ds1(3,natoms)  ! the calculated surface derivatives
real(kind=8)::d2s0(3,natoms,3,natoms),d2s1(3,natoms,3,natoms) ! the calculated surface hessians
real(kind=8)::com(sum_eds,3)   ! array with center of mass coordinates
real(kind=8)::Rinv ! the inverse bond length/distance (for derivatives)
real(kind=8)::r_factor(sum_eds,sum_eds)  ! modifying factor in s0 for different distances
real(kind=8)::correct ! local correction for ds0 and d2s0
real(kind=8)::massfactor ! for hessians: factor of two masses
real(kind=8)::dxx,dxy,dxz,dyy,dyz,dzz ! abbreviations for hessian elements
real(kind=8)::xi_ideal  ! the reference xi value (for recrossing)
integer::atom,atom1,atom2 ! the actual atomic indices
!
!    Here methods to calculate the Xi, dXi and d2Xi values/arrays are 
!    listed for all supported types of reactions 
! 
!    Currently these are:
!    1. Bimolecular reaction (one forming, one breaking bond)
!    2. Cycloaddition (two forming bonds)
!    3. Addition (two forming, one breaking bonds)
!    4. Addition3 (three forming, two breaking bonds)
!
!    #############################################################################
!    1 THE BIMOLECULAR REACTION (ONE FORMING, ONE BREAKING)
!
if ((umbr_type .eq. "BIMOLEC") .or. (umbr_type .eq. "CYCLOADD") &
    & .or. (umbr_type .eq. "ADDITION") .or. (umbr_type .eq. "ADDITION3") &
    & .or. (umbr_type .eq. "ADDITION4") .or. (umbr_type .eq. "ADD3_SOLV") &
    & .or. (umbr_type .eq. "ADD4_SOLV") .or. (umbr_type .eq. "MERGING")) then
!
!    STEP1: calculate the bondlengths and then the Xi-value
!
!
!    the forming bonds (in total form_num bonds)
!
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      R_f(1,i)=coords(1,atom1) - coords(1,atom2)
      R_f(2,i)=coords(2,atom1) - coords(2,atom2)
      R_f(3,i)=coords(3,atom1) - coords(3,atom2)
      form_act(i)=sqrt(R_f(1,i)*R_f(1,i)+R_f(2,i)*R_f(2,i)+R_f(3,i)*R_f(3,i))
   end do
!
!    the breaking bonds (in total break_num bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      R_b(1,i)=coords(1,atom1) - coords(1,atom2)
      R_b(2,i)=coords(2,atom1) - coords(2,atom2)
      R_b(3,i)=coords(3,atom1) - coords(3,atom2)
      break_act(i)=sqrt(R_b(1,i)*R_b(1,i)+R_b(2,i)*R_b(2,i)+R_b(3,i)*R_b(3,i))
   end do

!
!    value of the TS dividing surface funtion: sum up all breaking 
!    bond differences and divide this through the total number of 
!    breaking and sum up all forming bond differences and divide this 
!    through the total numbe of forming
!
   s1=0.d0
   do i=1,break_num
      s1=s1+(break_act(i)-break_ref(i))/real(break_num)
   end do
   do i=1,form_num
      s1=s1-(form_act(i)-form_ref(i))/real(form_num)
   end do
   
!
!    value of the reactant dividing surface function
!    calculate the difference between center of masses of both reactands 
!    and compare it to the predefined R_inf distance (no interaction)
!
!   write(*,*) "coords2",coords  
   call calc_com(coords,com)
!
!    calclate the R_ed matrix will all needed combination 
!    (upper triangular matrix)
!    and apply its elements to the s0 coordinate
!
   r_eds=0.d0
   s0=0.d0 
   do i=1,sum_eds
      do j=i+1,sum_eds
         Red(1,i,j)=com(j,1)-com(i,1)
         Red(2,i,j)=com(j,2)-com(i,2)
         Red(3,i,j)=com(j,3)-com(i,3)
         r_eds(i,j)=sqrt(Red(1,i,j)*Red(1,i,j) + Red(2,i,j)*Red(2,i,j) &
            & + Red(3,i,j)*Red(3,i,j)) 
         s0=s0+(R_inf-r_eds(i,j))       
      end do
   end do
!   write(*,*) "s0 calculation",s0
!
!    divide the surface function through the number of interactions 
!    (n_ed^2-n_ed)/2 
!
   s0_terms=(sum_eds*sum_eds-sum_eds)/2
   s0=s0/real(s0_terms)
!
!     Try to enforce similar values for all possible educt-educt distances 
!     calculate average and the deviation to it to correct the actual values
! 
!   eds_avg=sum(r_eds)/real(s0_terms)
!   do i=1,sum_eds
!      do j=i+1,sum_eds
!          r_factor(i,j)=1/((2*r_eds(i,j)-eds_avg)/r_eds(i,j) )
!          write(*,*) i,j,eds_avg,r_eds(i,j),r_factor(i,j)
!         r_eds(i,j)=eds_avg !r_eds(i,j)+(-eds_avg+r_eds(i,j))      
!      end do
!   end do
!   write(*,*) s0_terms,sum(r_eds) 
!   write(*,*) "r12,r13,r23",r_eds(1,2),r_eds(1,3),r_eds(2,3),eds_avg
!   write(*,*) "s0,s1",s0,s1   
!
!    calculate Xi value with the RPMDrate formula
!
!    for usual umbrella sampling
!
   if (mode .eq. 1) then
      xi_act=s0/(s0-s1)
!      write(*,*) "xi_act",xi_ideal,xi_act,s0,s1
!
!    for recrossing calculation
!
   else if (mode .eq. 2) then
      xi_act=xi_ideal*s1+(1-xi_ideal)*s0
!      write(*,*) "xi_act",xi_act,xi_ideal,s0,s1
   end if
!
!     STEP2: calculate the gradient of reaction coordinate
!
!
!     The transition state dividing surface
!
   ds1=0.d0
!
!     the forming bonds (in total, n_form forming bonds)
!
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      Rinv=1.d0/form_act(i)
      ds1(1,atom1) = ds1(1,atom1) - R_f(1,i) * Rinv/real(form_num)
      ds1(2,atom1) = ds1(2,atom1) - R_f(2,i) * Rinv/real(form_num)
      ds1(3,atom1) = ds1(3,atom1) - R_f(3,i) * Rinv/real(form_num)
      ds1(1,atom2) = ds1(1,atom2) + R_f(1,i) * Rinv/real(form_num)
      ds1(2,atom2) = ds1(2,atom2) + R_f(2,i) * Rinv/real(form_num)
      ds1(3,atom2) = ds1(3,atom2) + R_f(3,i) * Rinv/real(form_num)
   end do
!
!     the breaking bonds (in total, n_break breaking bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      Rinv=1.d0/break_act(i)
      ds1(1,atom1) = ds1(1,atom1) + R_b(1,i) * Rinv/real(break_num)
      ds1(2,atom1) = ds1(2,atom1) + R_b(2,i) * Rinv/real(break_num)
      ds1(3,atom1) = ds1(3,atom1) + R_b(3,i) * Rinv/real(break_num)
      ds1(1,atom2) = ds1(1,atom2) - R_b(1,i) * Rinv/real(break_num)
      ds1(2,atom2) = ds1(2,atom2) - R_b(2,i) * Rinv/real(break_num)
      ds1(3,atom2) = ds1(3,atom2) - R_b(3,i) * Rinv/real(break_num)
   end do

! 
!     The educts dividing surface
!
   
   ds0=0.d0
   do i=1,sum_eds
      do j=i+1,sum_eds
         Rinv=1.d0/r_eds(i,j)
         do k=1,n_ed(i)
            atom=at_ed(i,k)
            ds0(1,atom)=ds0(1,atom)+Red(1,i,j)*Rinv*mass(atom)/mass_ed(i)/real(s0_terms)
            ds0(2,atom)=ds0(2,atom)+Red(2,i,j)*Rinv*mass(atom)/mass_ed(i)/real(s0_terms)
            ds0(3,atom)=ds0(3,atom)+Red(3,i,j)*Rinv*mass(atom)/mass_ed(i)/real(s0_terms)
         end do
         do k=1,n_ed(j)
            atom=at_ed(j,k)
            ds0(1,atom)=ds0(1,atom)-Red(1,i,j)*Rinv*mass(atom)/mass_ed(j)/real(s0_terms)
            ds0(2,atom)=ds0(2,atom)-Red(2,i,j)*Rinv*mass(atom)/mass_ed(j)/real(s0_terms)
            ds0(3,atom)=ds0(3,atom)-Red(3,i,j)*Rinv*mass(atom)/mass_ed(j)/real(s0_terms)
         end do
      end do
   end do
!
!     For usual umbrella sampling
!
   if (mode .eq. 1) then
      dxi_act=(s0*ds1-s1*ds0)/((s0-s1)*(s0-s1))
!
!     For recrossing calculation
!
   else if (mode .eq. 2) then
      dxi_act=xi_ideal*ds1+(1-xi_ideal)*ds0
   end if

!
!     STEP3: calculate the hessian of the reaction coordinate
!      (I don´t know why this is needed!)
!
!
!     The transition state dividing surface
!
   d2s1=0.d0
!
!     the forming bonds (in total form_num forming bonds)
!
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      Rinv=1.d0/form_act(i)

      dxx = -(R_f(2,i) * R_f(2,i) + R_f(3,i) * R_f(3,i)) * (Rinv * Rinv * Rinv)
      dyy = -(R_f(3,i) * R_f(3,i) + R_f(1,i) * R_f(1,i)) * (Rinv * Rinv * Rinv)
      dzz = -(R_f(1,i) * R_f(1,i) + R_f(2,i) * R_f(2,i)) * (Rinv * Rinv * Rinv)
      dxy = R_f(1,i) * R_f(2,i) * (Rinv * Rinv * Rinv)
      dxz = R_f(1,i) * R_f(3,i) * (Rinv * Rinv * Rinv)
      dyz = R_f(2,i) * R_f(3,i) * (Rinv * Rinv * Rinv)

      d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx/real(form_num)
      d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy/real(form_num)
      d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz/real(form_num)
      d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy/real(form_num)
      d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy/real(form_num)
      d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz/real(form_num)
      d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz/real(form_num)
      d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz/real(form_num)
      d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz/real(form_num)

      d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx/real(form_num)
      d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy/real(form_num)
      d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz/real(form_num)
      d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy/real(form_num)
      d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy/real(form_num)
      d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz/real(form_num)
      d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz/real(form_num)
      d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz/real(form_num)
      d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz/real(form_num)

      d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx/real(form_num)
      d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy/real(form_num)
      d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz/real(form_num)
      d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy/real(form_num)
      d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy/real(form_num)
      d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz/real(form_num)
      d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz/real(form_num)
      d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz/real(form_num)
      d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz/real(form_num)

      d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx/real(form_num)
      d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy/real(form_num)
      d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz/real(form_num)
      d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy/real(form_num)
      d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy/real(form_num)
      d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz/real(form_num)
      d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz/real(form_num)
      d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz/real(form_num)
      d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz/real(form_num)
   end do
!
!    the breaking bonds (in total break_num breaking bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      Rinv=1.d0/break_act(i)

      dxx = (R_b(2,i) * R_b(2,i) + R_b(3,i) * R_b(3,i)) * (Rinv * Rinv * Rinv)
      dyy = (R_b(3,i) * R_b(3,i) + R_b(1,i) * R_b(1,i)) * (Rinv * Rinv * Rinv)
      dzz = (R_b(1,i) * R_b(1,i) + R_b(2,i) * R_b(2,i)) * (Rinv * Rinv * Rinv)
      dxy = -R_b(1,i) * R_b(2,i) * (Rinv * Rinv * Rinv)
      dxz = -R_b(1,i) * R_b(3,i) * (Rinv * Rinv * Rinv)
      dyz = -R_b(2,i) * R_b(3,i) * (Rinv * Rinv * Rinv)

      d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx/real(break_num)
      d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy/real(break_num)
      d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz/real(break_num)
      d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy/real(break_num)
      d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy/real(break_num)
      d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz/real(break_num)
      d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz/real(break_num)
      d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz/real(break_num)
      d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz/real(break_num)

      d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx/real(break_num)
      d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy/real(break_num)
      d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz/real(break_num)
      d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy/real(break_num)
      d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy/real(break_num)
      d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz/real(break_num)
      d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz/real(break_num)
      d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz/real(break_num)
      d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz/real(break_num)

      d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx/real(break_num)
      d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy/real(break_num)
      d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz/real(break_num)
      d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy/real(break_num)
      d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy/real(break_num)
      d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz/real(break_num)
      d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz/real(break_num)
      d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz/real(break_num)
      d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz/real(break_num)

      d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx/real(break_num)
      d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy/real(break_num)
      d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz/real(break_num)
      d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy/real(break_num)
      d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy/real(break_num)
      d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz/real(break_num)
      d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz/real(break_num)
      d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz/real(break_num)
      d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz/real(break_num)
   end do
!
!     The educts dividing surface
!

   d2s0=0.d0
   do i=1,sum_eds
      do j=i+1,sum_eds   
         Rinv=1.d0/r_eds(i,j)
         dxx = -(Red(2,i,j) * Red(2,i,j) + Red(3,i,j) * Red(3,i,j)) * (Rinv * Rinv * Rinv)
         dyy = -(Red(3,i,j) * Red(3,i,j) + Red(1,i,j) * Red(1,i,j)) * (Rinv * Rinv * Rinv)
         dzz = -(Red(1,i,j) * Red(1,i,j) + Red(2,i,j) * Red(2,i,j)) * (Rinv * Rinv * Rinv)
         dxy = Red(1,i,j) * Red(2,i,j) * (Rinv * Rinv * Rinv)
         dxz = Red(1,i,j) * Red(3,i,j) * (Rinv * Rinv * Rinv)
         dyz = Red(2,i,j) * Red(3,i,j) * (Rinv * Rinv * Rinv)
   
         do k=1,n_ed(i)
            atom1=at_ed(i,k)
            do l=1,n_ed(i)
               atom2=at_ed(i,l)
               massfactor=mass(atom1)/mass_ed(i)*mass(atom2)/mass_ed(i)
               d2s0(1,atom1,1,atom2) = d2s0(1,atom1,1,atom2)+dxx * massfactor/real(s0_terms)
               d2s0(1,atom1,2,atom2) = d2s0(1,atom1,2,atom2)+dxy * massfactor/real(s0_terms)
               d2s0(1,atom1,3,atom2) = d2s0(1,atom1,3,atom2)+dxz * massfactor/real(s0_terms)
               d2s0(2,atom1,1,atom2) = d2s0(2,atom1,1,atom2)+dxy * massfactor/real(s0_terms)
               d2s0(2,atom1,2,atom2) = d2s0(2,atom1,2,atom2)+dyy * massfactor/real(s0_terms)
               d2s0(2,atom1,3,atom2) = d2s0(2,atom1,3,atom2)+dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,1,atom2) = d2s0(3,atom1,1,atom2)+dxz * massfactor/real(s0_terms)
               d2s0(3,atom1,2,atom2) = d2s0(3,atom1,2,atom2)+dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,3,atom2) = d2s0(3,atom1,3,atom2)+dzz * massfactor/real(s0_terms)
            end do
            do l=1,n_ed(j)
               atom2=at_ed(j,l)
               massfactor=mass(atom1)/mass_ed(i)*mass(atom2)/mass_ed(j)
               d2s0(1,atom1,1,atom2) = d2s0(1,atom1,1,atom2)-dxx * massfactor/real(s0_terms)
               d2s0(1,atom1,2,atom2) = d2s0(1,atom1,2,atom2)-dxy * massfactor/real(s0_terms)
               d2s0(1,atom1,3,atom2) = d2s0(1,atom1,3,atom2)-dxz * massfactor/real(s0_terms)
               d2s0(2,atom1,1,atom2) = d2s0(2,atom1,1,atom2)-dxy * massfactor/real(s0_terms)
               d2s0(2,atom1,2,atom2) = d2s0(2,atom1,2,atom2)-dyy * massfactor/real(s0_terms)
               d2s0(2,atom1,3,atom2) = d2s0(2,atom1,3,atom2)-dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,1,atom2) = d2s0(3,atom1,1,atom2)-dxz * massfactor/real(s0_terms)
               d2s0(3,atom1,2,atom2) = d2s0(3,atom1,2,atom2)-dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,3,atom2) = d2s0(3,atom1,3,atom2)-dzz * massfactor/real(s0_terms)
            end do
         end do
         do k=1,n_ed(j)
            atom1=at_ed(j,k)
            do l=1,n_ed(i)
               atom2=at_ed(i,l)
               massfactor=mass(atom1)/mass_ed(j)*mass(atom2)/mass_ed(i)
               d2s0(1,atom1,1,atom2) = d2s0(1,atom1,1,atom2)-dxx * massfactor/real(s0_terms)
               d2s0(1,atom1,2,atom2) = d2s0(1,atom1,2,atom2)-dxy * massfactor/real(s0_terms)
               d2s0(1,atom1,3,atom2) = d2s0(1,atom1,3,atom2)-dxz * massfactor/real(s0_terms)
               d2s0(2,atom1,1,atom2) = d2s0(2,atom1,1,atom2)-dxy * massfactor/real(s0_terms)
               d2s0(2,atom1,2,atom2) = d2s0(2,atom1,2,atom2)-dyy * massfactor/real(s0_terms)
               d2s0(2,atom1,3,atom2) = d2s0(2,atom1,3,atom2)-dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,1,atom2) = d2s0(3,atom1,1,atom2)-dxz * massfactor/real(s0_terms)
               d2s0(3,atom1,2,atom2) = d2s0(3,atom1,2,atom2)-dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,3,atom2) = d2s0(3,atom1,3,atom2)-dzz * massfactor/real(s0_terms)
            end do
            do l=1,n_ed(j)
               atom2=at_ed(j,l)
               massfactor=mass(atom1)/mass_ed(j)*mass(atom2)/mass_ed(j)
               d2s0(1,atom1,1,atom2) = d2s0(1,atom1,1,atom2)+dxx * massfactor/real(s0_terms)
               d2s0(1,atom1,2,atom2) = d2s0(1,atom1,2,atom2)+dxy * massfactor/real(s0_terms)
               d2s0(1,atom1,3,atom2) = d2s0(1,atom1,3,atom2)+dxz * massfactor/real(s0_terms)
               d2s0(2,atom1,1,atom2) = d2s0(2,atom1,1,atom2)+dxy * massfactor/real(s0_terms)
               d2s0(2,atom1,2,atom2) = d2s0(2,atom1,2,atom2)+dyy * massfactor/real(s0_terms)
               d2s0(2,atom1,3,atom2) = d2s0(2,atom1,3,atom2)+dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,1,atom2) = d2s0(3,atom1,1,atom2)+dxz * massfactor/real(s0_terms)
               d2s0(3,atom1,2,atom2) = d2s0(3,atom1,2,atom2)+dyz * massfactor/real(s0_terms)
               d2s0(3,atom1,3,atom2) = d2s0(3,atom1,3,atom2)+dzz * massfactor/real(s0_terms)
            end do
         end do
      end do
   end do
!
!     Determine Hessian of the reaction coordinate from dividing surface 
!     hessians
!
!     for usual umbrella sampling
!
   if (mode .eq. 1) then
      do i1 = 1, 3
         do j1 = 1, Natoms
            do i2 = 1, 3
               do j2 = 1, Natoms
                   d2xi_act(i1,j1,i2,j2) = ((s0 * d2s1(i1,j1,i2,j2) + ds0(i2,j2) * ds1(i1,j1) &
                     &  - ds1(i2,j2) * ds0(i1,j1) - s1 * d2s0(i1,j1,i2,j2)) * (s0 - s1) &
                     &  - 2.0d0 * (s0 * ds1(i1,j1) - s1 * ds0(i1,j1)) &
                     &  * (ds0(i2,j2) - ds1(i2,j2))) &
                     &  / ((s0 - s1) * (s0 - s1) * (s0 - s1))
               end do
            end do
         end do
      end do
!
!     For recrossing factor calculation
!
   else if (mode .eq. 2) then
      d2xi_act=xi_ideal*d2s1+(1-xi_ideal)*d2s0
   end if
!   write(*,*) "Xiii values!",mode
!   write(*,*) "xi_act",xi_act,xi_ideal,s0,s1
!   write(*,*) xi_act
!   write(*,*) "dxi_act"
!   write(*,*) dxi_act
!            write(*,*) "d2s1!"
!            write(*,*) d2s1
!            write(*,*) "d2s0!"
!            write(*,*) d2s0
!
!   write(*,*) "d2xi_act"
!   write(*,*) d2xi_act
!   stop "gizfitcditucD"   
!
!    #############################################################################
!    2 THE ATOM SHIFT (ONE COORDINATE RESTRAINT ON ONE ATOM)
!
!    These coordinate is rather simple: for s0, the expression is x_real-x_ideal,
!    the first derivative is 1 for the relevent coordinate, the second is zero
!
else if (umbr_type .eq. "ATOM_SHIFT") then
!
!    value of the TS dividing surface function: reference to shift_hi  
!
   s1=coords(shift_coord,shift_atom)-shift_hi
!
!    value of the educts dividing surface function: reference to shift_lo
!
   s0=coords(shift_coord,shift_atom)-shift_lo
!
!    Calculate the actual Xi value:
!    for usual umbrella sampling
!
   if (mode .eq. 1) then
      xi_act=s0/(s0-s1)
!
!    for recrossing calculation
!
   else if (mode .eq. 2) then
      xi_act=xi_ideal*s1+(1-xi_ideal)*s0
!      write(*,*) "xi_act",xi_act,xi_ideal,s0,s1
   end if
!
!    gradient of the TS dividing surface function:
!
   ds1=0d0
   ds1(shift_coord,shift_atom)=1d0
!
!    gradient of the educts dividing surface function:
!   
   ds0=0d0
   ds0(shift_coord,shift_atom)=1d0
!
!    calculate the actual Xi derivative:
!
!
!     For usual umbrella sampling
!
   if (mode .eq. 1) then
      dxi_act=(s0*ds1-s1*ds0)/((s0-s1)*(s0-s1))
!
!     For recrossing calculation
!
   else if (mode .eq. 2) then
      dxi_act=xi_ideal*ds1+(1-xi_ideal)*ds0
   end if
!
!     hession of reaction coordinate: both d2s0 and d2s1 are zero!
!
   d2s0=0.d0
   d2s1=0.d0
!
!     Determine Hessian of the reaction coordinate from dividing surface 
!     hessians
!
!     for usual umbrella sampling
!
   if (mode .eq. 1) then
      do i1 = 1, 3
         do j1 = 1, Natoms
            do i2 = 1, 3
               do j2 = 1, Natoms
                   d2xi_act(i1,j1,i2,j2) = ((s0 * d2s1(i1,j1,i2,j2) + ds0(i2,j2) * ds1(i1,j1) &
                     &  - ds1(i2,j2) * ds0(i1,j1) - s1 * d2s0(i1,j1,i2,j2)) * (s0 - s1) &
                     &  - 2.0d0 * (s0 * ds1(i1,j1) - s1 * ds0(i1,j1)) &
                     &  * (ds0(i2,j2) - ds1(i2,j2))) &
                     &  / ((s0 - s1) * (s0 - s1) * (s0 - s1))
               end do
            end do
         end do
      end do
!
!     For recrossing factor calculation
!
   else if (mode .eq. 2) then
      d2xi_act=xi_ideal*d2s1+(1-xi_ideal)*d2s0
   end if
 
!    #############################################################################
!    3 UNIMOLECULAR REACTIONS (GENERAL)
!
!    Here, the s0 dividing surface is defined as the educt structure. Therefore, 
!    both dividing surfaces are defined as a bunch of reference values for 
!    the forming and breaking bonds of interest
!
!

else if ((umbr_type .eq. "CYCLOREVER") .or. (umbr_type .eq. "REARRANGE") .or. &
                  & (umbr_type .eq. "DECOM_1BOND") .or. (umbr_type .eq. "ELIMINATION")) then
!
!    STEP1: calculate the bondlengths and then the Xi-value
!
!    a) the s1 surface: reference to the TS
!
!    the forming bonds (in total form_num bonds)
!
!   write(*,*) "form",form_num,break_num
!   stop "hgpgu"
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      R_f(1,i)=coords(1,atom1) - coords(1,atom2)
      R_f(2,i)=coords(2,atom1) - coords(2,atom2)
      R_f(3,i)=coords(3,atom1) - coords(3,atom2)
      form_act(i)=sqrt(R_f(1,i)*R_f(1,i)+R_f(2,i)*R_f(2,i)+R_f(3,i)*R_f(3,i))
   end do
!
!    the breaking bonds (in total break_num bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      R_b(1,i)=coords(1,atom1) - coords(1,atom2)
      R_b(2,i)=coords(2,atom1) - coords(2,atom2)
      R_b(3,i)=coords(3,atom1) - coords(3,atom2)
      break_act(i)=sqrt(R_b(1,i)*R_b(1,i)+R_b(2,i)*R_b(2,i)+R_b(3,i)*R_b(3,i))
   end do
!
!    value of the TS dividing surface funtion: sum up all breaking 
!    bond differences and divide this through the total number of 
!    breaking and sum up all forming bond differences and divide this 
!    through the total numbe of forming
!
   s1=0.d0
   do i=1,break_num
      s1=s1+(break_act(i)-break_ref(i))/real(break_num)
   end do
   do i=1,form_num
      s1=s1-(form_act(i)-form_ref(i))/real(form_num)
   end do
!     
!    value of the educts dividing surface: analogous to TS dividing surface 
!
   s0=0.d0
   do i=1,break_num
      s0=s0+(break_act(i)-break_ed(i))/real(break_num)
   end do
   do i=1,form_num
      s0=s0-(form_act(i)-form_ed(i))/real(form_num)
   end do
!
!    calculate Xi value with the RPMDrate formula
!
!    for usual umbrella sampling
!
   if (mode .eq. 1) then
      xi_act=s0/(s0-s1)
!
!    for recrossing calculation
!
   else if (mode .eq. 2) then
      xi_act=xi_ideal*s1+(1-xi_ideal)*s0
   end if
!
!     STEP2: calculate the gradient of reaction coordinate
!
!
!     The transition state dividing surface
!
   ds1=0.d0
!
!     the forming bonds (in total, n_form forming bonds)
!
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      Rinv=1.d0/form_act(i)
      ds1(1,atom1) = ds1(1,atom1) - R_f(1,i) * Rinv/real(form_num)
      ds1(2,atom1) = ds1(2,atom1) - R_f(2,i) * Rinv/real(form_num)
      ds1(3,atom1) = ds1(3,atom1) - R_f(3,i) * Rinv/real(form_num)
      ds1(1,atom2) = ds1(1,atom2) + R_f(1,i) * Rinv/real(form_num)
      ds1(2,atom2) = ds1(2,atom2) + R_f(2,i) * Rinv/real(form_num)
      ds1(3,atom2) = ds1(3,atom2) + R_f(3,i) * Rinv/real(form_num)
   end do
!
!     the breaking bonds (in total, n_break breaking bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      Rinv=1.d0/break_act(i)
      ds1(1,atom1) = ds1(1,atom1) + R_b(1,i) * Rinv/real(break_num)
      ds1(2,atom1) = ds1(2,atom1) + R_b(2,i) * Rinv/real(break_num)
      ds1(3,atom1) = ds1(3,atom1) + R_b(3,i) * Rinv/real(break_num)
      ds1(1,atom2) = ds1(1,atom2) - R_b(1,i) * Rinv/real(break_num)
      ds1(2,atom2) = ds1(2,atom2) - R_b(2,i) * Rinv/real(break_num)
      ds1(3,atom2) = ds1(3,atom2) - R_b(3,i) * Rinv/real(break_num)
   end do
!
!     educts dividing surface: derivatives are equal!
!
   ds0=ds1
!
!     Calculate the derivative of the Xi coordinate 
!
!     For usual umbrella sampling
!
   if (mode .eq. 1) then
      dxi_act=(s0*ds1-s1*ds0)/((s0-s1)*(s0-s1))
!
!     For recrossing calculation
!
   else if (mode .eq. 2) then
      dxi_act=xi_ideal*ds1+(1-xi_ideal)*ds0
   end if
!
!     STEP3: calculate the hessian of the reaction coordinate
!      (I don´t know why this is needed!)
!
!
!     The transition state dividing surface
!
   d2s1=0.d0
!
!     the forming bonds (in total form_num forming bonds)
!
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      Rinv=1.d0/form_act(i)

      dxx = -(R_f(2,i) * R_f(2,i) + R_f(3,i) * R_f(3,i)) * (Rinv * Rinv * Rinv)
      dyy = -(R_f(3,i) * R_f(3,i) + R_f(1,i) * R_f(1,i)) * (Rinv * Rinv * Rinv)
      dzz = -(R_f(1,i) * R_f(1,i) + R_f(2,i) * R_f(2,i)) * (Rinv * Rinv * Rinv)
      dxy = R_f(1,i) * R_f(2,i) * (Rinv * Rinv * Rinv)
      dxz = R_f(1,i) * R_f(3,i) * (Rinv * Rinv * Rinv)
      dyz = R_f(2,i) * R_f(3,i) * (Rinv * Rinv * Rinv)

      d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx/real(form_num)
      d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy/real(form_num)
      d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz/real(form_num)
      d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy/real(form_num)
      d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy/real(form_num)
      d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz/real(form_num)
      d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz/real(form_num)
      d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz/real(form_num)
      d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz/real(form_num)

      d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx/real(form_num)
      d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy/real(form_num)
      d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz/real(form_num)
      d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy/real(form_num)
      d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy/real(form_num)
      d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz/real(form_num)
      d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz/real(form_num)
      d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz/real(form_num)
      d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz/real(form_num)

      d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx/real(form_num)
      d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy/real(form_num)
      d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz/real(form_num)
      d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy/real(form_num)
      d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy/real(form_num)
      d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz/real(form_num)
      d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz/real(form_num)
      d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz/real(form_num)
      d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz/real(form_num)

      d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx/real(form_num)
      d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy/real(form_num)
      d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz/real(form_num)
      d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy/real(form_num)
      d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy/real(form_num)
      d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz/real(form_num)
      d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz/real(form_num)
      d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz/real(form_num)
      d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz/real(form_num)
   end do
!
!    the breaking bonds (in total break_num breaking bonds)
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      Rinv=1.d0/break_act(i)

      dxx = (R_b(2,i) * R_b(2,i) + R_b(3,i) * R_b(3,i)) * (Rinv * Rinv * Rinv)
      dyy = (R_b(3,i) * R_b(3,i) + R_b(1,i) * R_b(1,i)) * (Rinv * Rinv * Rinv)
      dzz = (R_b(1,i) * R_b(1,i) + R_b(2,i) * R_b(2,i)) * (Rinv * Rinv * Rinv)
      dxy = -R_b(1,i) * R_b(2,i) * (Rinv * Rinv * Rinv)
      dxz = -R_b(1,i) * R_b(3,i) * (Rinv * Rinv * Rinv)
      dyz = -R_b(2,i) * R_b(3,i) * (Rinv * Rinv * Rinv)

      d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx/real(break_num)
      d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy/real(break_num)
      d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz/real(break_num)
      d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy/real(break_num)
      d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy/real(break_num)
      d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz/real(break_num)
      d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz/real(break_num)
      d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz/real(break_num)
      d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz/real(break_num)

      d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx/real(break_num)
      d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy/real(break_num)
      d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz/real(break_num)
      d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy/real(break_num)
      d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy/real(break_num)
      d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz/real(break_num)
      d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz/real(break_num)
      d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz/real(break_num)
      d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz/real(break_num)

      d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx/real(break_num)
      d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy/real(break_num)
      d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz/real(break_num)
      d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy/real(break_num)
      d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy/real(break_num)
      d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz/real(break_num)
      d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz/real(break_num)
      d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz/real(break_num)
      d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz/real(break_num)

      d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx/real(break_num)
      d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy/real(break_num)
      d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz/real(break_num)
      d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy/real(break_num)
      d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy/real(break_num)
      d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz/real(break_num)
      d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz/real(break_num)
      d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz/real(break_num)
      d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz/real(break_num)
   end do
!
!     The Hessian of the educts dividing surface: equal to the TS one!
!
   d2s0=d2s1
!
!     Determine Hessian of the reaction coordinate from dividing surface 
!     hessians
!
!     for usual umbrella sampling
!
   if (mode .eq. 1) then
      do i1 = 1, 3
         do j1 = 1, Natoms
            do i2 = 1, 3
               do j2 = 1, Natoms
                   d2xi_act(i1,j1,i2,j2) = ((s0 * d2s1(i1,j1,i2,j2) + ds0(i2,j2) * ds1(i1,j1) &
                     &  - ds1(i2,j2) * ds0(i1,j1) - s1 * d2s0(i1,j1,i2,j2)) * (s0 - s1) &
                     &  - 2.0d0 * (s0 * ds1(i1,j1) - s1 * ds0(i1,j1)) &
                     &  * (ds0(i2,j2) - ds1(i2,j2))) &
                     &  / ((s0 - s1) * (s0 - s1) * (s0 - s1))
               end do
            end do
         end do
      end do
!
!     For recrossing factor calculation
!
   else if (mode .eq. 2) then
      d2xi_act=xi_ideal*d2s1+(1-xi_ideal)*d2s0
   end if

end if
return
end subroutine calc_xi
