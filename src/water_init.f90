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
!     Subroutine water_init: Initializes a water model (currently only flexible 
!       SPC available) and checks if the structure includes indeed only 
!       water molecules
!
!     part of EVB
!

subroutine water_init
use evb_mod
use general
use qmdff
implicit none
!     Loop indices
integer::i,j
!
!     Check the atom types: if correctly ordered water molecules are 
!     included in the structure
!
if (mod(natoms,3) .ne. 0) then
   write(*,*) "The number of atoms in the system is no multiple of three!"
   write(*,*) "It seems that no pure water is there.."
   call fatal
end if
nwater=int(natoms/3)
allocate(water_act(natoms))
!
!     Check if the atoms are in the correct order (O, H, H, ...)
!
do i=1,nwater
   if (name((i-1)*3+1) .ne. "O") then
      write(*,*) "The atom ",(i-1)*3+1," should be oxygen!"
   else if (name((i-1)*3+2) .ne. "H") then
      write(*,*) "The atom ",(i-1)*3+2," should be hydrogen!"
   else if (name((i-1)*3+2) .ne. "H") then
      write(*,*) "The atom ",(i-1)*3+3," should be hydrogen!"
   end if
   water_act((i-1)*3+1)=i
   water_act((i-1)*3+2)=i
   water_act((i-1)*3+3)=i
end do

!
!     allocate the water parameters
!
allocate(water_pars(11))
water_pars(1)=1.0d0  ! r_0 (Angstrom)
water_pars(2)=1.633d0 ! r_0HH (Angstrom)
!water_pars(2)=109.47d0 ! theta_0 (degrees)
water_pars(3)=101.9188d0  ! D_e (kcal/mol)
water_pars(4)=2.567d0 ! a (1/Angstrom)
water_pars(5)=328.645606d0 ! k_theta (kcal/mol/Angstrom^2)
water_pars(6)=-211.4672d0 ! k_{r,theta} (kcal/mol/Angstrom^2)
water_pars(7)=111.70765 ! k_rr (kcal/mol/Angstrom^2)
water_pars(8)=0.41d0 ! e_H (e_0)
water_pars(9)=-0.82d0 ! e_O (e_0)
water_pars(10)=3.166d0 ! sigma_OO (Angstrom)
water_pars(11)=0.1554 ! epsilon_OO (kcal/mol)
!
!     Convert the water parameters to atomic units 
!
water_pars(1)=water_pars(1)/bohr
water_pars(2)=water_pars(2)/bohr
!water_pars(2)=water_pars(2)*2d0*pi/180.d0
water_pars(3)=water_pars(3)/hartree
water_pars(4)=water_pars(4)*bohr
water_pars(5)=water_pars(5)/hartree*bohr*bohr
water_pars(6)=water_pars(6)/hartree*bohr*bohr
water_pars(7)=water_pars(7)/hartree*bohr*bohr
water_pars(8)=water_pars(8)
water_pars(9)=water_pars(9)
water_pars(10)=water_pars(10)/bohr
water_pars(11)=water_pars(11)/hartree
!
!     Store the charges of all atoms 
!
allocate(q(natoms))
do i=1,nwater
   q((i-1)*3+1)=water_pars(9)
   q((i-1)*3+2)=water_pars(8)
   q((i-1)*3+3)=water_pars(8)
end do



return
end subroutine water_init
