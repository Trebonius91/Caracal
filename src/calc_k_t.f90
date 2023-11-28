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
!     Subroutine calc_k_t: calculate the rate constant for the 
!       reaction using the results from the rpmd calculation
!       QTST is used, the formula varies with the type of the 
!       calculation (unimolecular, bimolecular etc..)
!  
!     part of EVB
!
subroutine calc_k_t(pmf_max,pmf_min,xi_max,xi_min,kappa,pmf_int)
use general
use evb_mod 
implicit none
real(kind=8)::xi_max,xi_min  ! positions of minimum/maximum of the PMF
real(kind=8)::pmf_max,pmf_min  ! the minimum and maximum of the PMF 
real(kind=8)::pmf_int  ! integral from xi_0 to xi_TS (only unimolecular)
real(kind=8)::kappa    ! the calculated recrossing factor
real(kind=8)::my_R  ! the reduced mass of the reactants
real(kind=8)::k_t  ! the final rate constant
real(kind=8)::pi  ! pi

pi=3.1415926535897932384d0    ! the pi
!
!     for the bimolecular reactions: use the formula from RPMDrate
!
if ((umbr_type .eq. "BIMOLEC") .or. (umbr_type .eq. "CYCLOADD") .or. &
           & (umbr_type .eq. "BIMOL_EXCH") .or. &
           & (umbr_type .eq. "ADDITION") .or. (umbr_type .eq. "ADDITION3") .or. &
           & (umbr_type .eq. "ADDITION4") .or. (umbr_type .eq. "ADD3_SOLV") .or. &
           & (umbr_type .eq. "ADD4_SOLV") .or. (umbr_type .eq. "MERGING")) then
   if (umbr_type .eq. "BIMOLEC") then
      write(15,*) "Type of reaction: BIMOLECULAR SUBSTITUTION"
      write(*,*) "Type of reaction: BIMOLECULAR SUBSTITUTION"
   else if (umbr_type .eq. "BIMOL_EXCH") then
      write(15,*) "Type of reaction: BIMOLECULAR CYCLIC EXCHANGE"
      write(*,*) "Type of reaction: BIMOLECULAR CYCLIC EXCHANGE"
   else if (umbr_type .eq. "CYCLOADD") then
      write(15,*) "Type of reaction: CYCLOADDITION"
      write(*,*) "Type of reaction: CYCLOADDITION"
   else if (umbr_type .eq. "MERGING") then
      write(15,*) "Type of reaction: BIMOLECULAR MERGING"
      write(*,*) "Type of reaction: BIMOLECULAR MERGING"
   else if (umbr_type .eq. "ADDITION") then
      write(15,*) "Type of reaction: BIMOLECULAR ADDITION"
      write(*,*) "Type of reaction: BIMOLECULAR ADDITION"
   else if (umbr_type .eq. "ADDITION3") then
      write(15,*) "Type of reaction: TRIOLECULAR ADDITION"
      write(*,*) "Type of reaction: TRIMOLECULAR ADDITION"
   else if (umbr_type .eq. "ADDITION4") then
      write(15,*) "Type of reaction: TETRAMOLECULAR ADDITION"
      write(*,*) "Type of reaction: TETRAMOLECULAR ADDITION"
   else if (umbr_type .eq. "ADD3_SOLV") then
      write(15,*) "Type of reaction: TRIMOLECULAR ADDITION/PRE-COMPLEX"
      write(*,*) "Type of reaction: TRIMOLECULAR ADDITION/PRE-COMPLEX"
   else if (umbr_type .eq. "ADD4_SOLV") then
      write(15,*) "Type of reaction: TETRAMOLECULAR ADDITION/PRE-COMPLEX"
      write(*,*) "Type of reaction: TETRAMOLECULAR ADDITION/PRE-COMPLEX"
   end if
   write(15,*)
   write(15,*) "- The rate formula is:"
   write(15,*)
   write(15,*) "k(T)=n_paths*kappa*4*pi*dist_inf**2*sqrt(1/(2*pi*beta*my_R))*exp(-beta*(W(xi_TS)-W(0)))"
   write(15,*) " where my_R=m1*m2/(m1+m2)"
   write(15,*)
   write(*,*)
   write(*,*) "- The rate formula is:" 
   write(*,*)
   write(*,*) "k(T)=n_paths*kappa*4*pi*dist_inf**2*sqrt(1/(2*pi*beta*my_R))*exp(-beta*(W(xi_TS)-W(0)))"
   write(*,*) " where my_R=m1*m2/(m1+m2)"
   write(*,*)
!
!     calculate the reduced mass of the reactants
!
   my_R=mass_reac(1)*mass_reac(2)/(mass_reac(1)+mass_reac(2))
   write(15,*) "- The parameters are defined as (in atomic units):"
   write(15,*)
   write(15,'(a,i3)') "   * n_paths = ",npaths
   write(15,'(a,f12.8)') "   * kappa = ",kappa
   write(15,'(a,f15.8,a)') "   * dist_inf = ",R_inf," bohr"
   write(15,'(a,f15.8,a)') "   * beta = ",beta," 1/hartree"
   write(15,'(a,f15.8,a)') "   * m1 = ",mass_reac(1)," amu"
   write(15,'(a,f15.8,a)') "   * m2 = ",mass_reac(2)," amu"
   write(15,'(a,f15.8,a)') "   * my_R = ",my_R," amu"
   write(15,'(a,f15.8,a,f12.8)') "   * W(xi_TS) = ",pmf_max," hartee at xi= ",xi_max
   write(15,'(a,f15.8,a,f12.8)') "   * w(xi_min) = ",pmf_min, " hartree at xi= ",xi_min
   write(15,'(a,f15.8,a)') "   * deltaW = ",(pmf_max-pmf_min)*hartree*joule, " kJ/mol"
   write(*,*) "- The parameters are defined as (in atomic units):"
   write(*,*)
   write(*,'(a,i3)') "   * n_paths = ",npaths
   write(*,'(a,f12.8)') "   * kappa = ",kappa 
   write(*,'(a,f15.8,a)') "   * dist_inf = ",R_inf," bohr"
   write(*,'(a,f15.8,a)') "   * beta = ",beta," 1/hartree"
   write(*,'(a,f15.8,a)') "   * m1 = ",mass_reac(1)," amu"
   write(*,'(a,f15.8,a)') "   * m2 = ",mass_reac(2)," amu"
   write(*,'(a,f15.8,a)') "   * my_R = ",my_R," amu"
   write(*,'(a,f15.8,a,f12.8)') "   * W(xi_TS) = ",pmf_max," hartee at xi= ",xi_max
   write(*,'(a,f15.8,a,f12.8)') "   * w(0) = ",pmf_min, " hartree at xi= ",xi_min
   write(*,'(a,f15.8,a)') "   * deltaW = ",(pmf_max-pmf_min)*hartree*joule, "kJ/mol"
!
!     calculate the rate constant in atomic units
!
   k_t=npaths*kappa*4.d0*pi*R_inf*R_inf*sqrt(1.d0/(2.d0*pi*beta*my_R))*exp(-beta*(pmf_max-pmf_min))
!
!     convert the rate constant into usual units (cm^3*mol/K)
!
   k_t=k_t*1e6 * ((5.2917721092e-11)**3 / 2.418884326505e-17)*6.02214179E23
   write(15,*)
   write(15,*) "----------------------------------------------"
   write(15,'(a)') "  The value of the reaction rate constant is: "
   write(15,'(a,es16.8,a)') "    ",k_t,"  cm^3/(mol*s)"
   write(15,'(a,es16.8,a)') "    ",k_t/avogadro,"  cm^3/(molec*s)"
   write(15,*) "----------------------------------------------"
   write(*,*) 
   write(*,*) "----------------------------------------------"
   write(*,'(a)') "  The value of the reaction rate constant is: "
   write(*,'(a,es16.8,a)') "    ",k_t,"  cm^3/(mol*s)"
   write(*,'(a,es16.8,a)') "    ",k_t/avogadro,"  cm^3/(molec*s)"
   write(*,*) "----------------------------------------------"

else if (umbr_type .eq. "ATOM_SHIFT") then
   write(*,*) "Type of reaction: Atom shift over a certain distance."
   write(*,*) "This is an UNIMOLECULAR REACTION!"
   write(*,*)
   write(*,*) "- The rate formula is:"
   write(*,*)
   write(*,*) "k(T)=kappa*1/sqrt(2*pi*m*beta)*exp(-beta*W(xi_TS))/(INT_(xi_min)^(xi_TS)) exp(-beta*W(s))ds)"
   write(*,*) "  where m=mass of the translated atom"
   write(*,*) "  and INT... the integrated PMF profile from the reactants to the TS"
   write(*,*) 
   write(15,*)
   write(15,*) "- The rate formula is:"
   write(15,*)
   write(15,*) "k(T)=kappa*1/sqrt(2*pi*m*beta)*exp(-beta*W(xi_TS))/(INT_(xi_min)^(xi_TS)) exp(-beta*W(s))ds)"
   write(15,*) "  where m=mass of the translated atom"
   write(15,*) "  and INT... the integrated PMF profile from the reactants to the TS"
   write(15,*)
   write(*,*) "- The parameters are defined as (in atomic units):"
   write(*,*) 
   write(*,'(a,f12.8)') "   * kappa = ",kappa
   write(*,'(a,f15.8,a)') "   * beta = ",beta," 1/hartree"
   write(*,'(a,f15.8,a)') "   * m = ",mass(shift_atom)," amu"   
   write(*,'(a,f15.8,a,f12.8)') "   * W(xi_TS) = ",pmf_max," hartee at xi= ",xi_max
   write(*,'(a,f15.8,a,f12.8)') "   * W(xi_min) = ",pmf_min, " hartree at xi= ",xi_min   
   write(*,'(a,f15.8)') "   * value of the integral from xi_min to xi_TS: ",pmf_int
!
!     calculate the rate constant in atomic units 
!
   k_t=kappa/sqrt((2*pi*beta*mass(shift_atom)))*exp(-beta*pmf_max)/pmf_int
!
!     convert the rate constant into usual units (s^(-1))
!
   k_t=k_t*6.02214179E23*2.418884326505e-17/5.2917721092e-11/pi
!   k_t=k_t*1e3 * ((5.2917721092e-11)**3 / 2.418884326505e-17)*6.02214179E23
   write(15,*)
   write(15,*) "----------------------------------------------"
   write(15,'(a)') "  The value of the reaction rate constant is: "
   write(15,'(a,es16.8,a)') "    ",k_t,"  s^(-1)"
   write(15,*) "----------------------------------------------"
   write(*,*)
   write(*,*) "----------------------------------------------"
   write(*,'(a)') "  The value of the reaction rate constant is: "
   write(*,'(a,es16.8,a)') "    ",k_t,"  s^(-1)"
   write(*,*) "----------------------------------------------"
else if ((umbr_type .eq. "CYCLOREVER") .or. (umbr_type .eq. "REARRANGE") .or. &
            &  (umbr_type .eq. "DECOM_1BOND") .or. (umbr_type .eq. "ELIMINATION")) then 
   if (umbr_type .eq. "CYCLOREVER") then
      write(*,*) "Type of reaction: Cycloreversion (two bonds broken)."
   else if (umbr_type .eq. "REARRANGE") then
      write(*,*) "Type of reaction: Rearrangement (one bond built and one broken)."
   else if (umbr_type .eq. "DECOM_1BOND") then
      write(*,*) "Type of reaction: simple Decomposition (one bond broken)."
   else if (umbr_type .eq. "ELIMINATION") then
      write(*,*) "Type of reaction: Elimination (two bonds broken and one built)."
   end if
   write(*,*) "This is an UNIMOLECULAR REACTION!"
   write(*,*)
   write(*,*) "So far, no general rate constant formula has been implemented for this case."
   write(*,*) "For the different mechanisms, different prefactors will be developed in the "
   write(*,*) " in the next version of CARACAL."
   write(*,*) "Until then, the simple TST rate formula (with kappa) will be used:"
   write(*,*)
   write(*,*) "k(T)=kappa*k_B*T/h*exp(-(W(xi_TS)-W(xi_min))/(k_B*T)) "
   write(*,*) " where T=temperature, k_B=Boltzmann constant, h=Plancks constant"
   write(*,*) "  W(xi_TS), W(xi_min)=free energies at the TS and the reactants minimum."
   write(15,*) "This is an UNIMOLECULAR REACTION!"
   write(15,*)
   write(15,*) "So far, no general rate constant formula has been implemented for this case."
   write(15,*) "For the different mechanisms, different prefactors will be derived in the "
   write(15,*) " in the next version of CARACAL."
   write(15,*) "Until then, the simple TST rate formula (with kappa) will be used:"
   write(15,*)
   write(15,*) "k(T)=kappa*n_paths*k_B*T/h*exp(-(W(xi_TS)-W(xi_min))/(k_B*T)) "
   write(15,*) " where T=temperature, k_B=Boltzmann constant, h=Plancks constant"
   write(15,*) "  W(xi_TS), W(xi_min)=free energies at the TS and the reactants minimum."

   write(*,*) "- The parameters are defined as (in atomic units):"
   write(*,*)
   write(*,'(a,f12.8)') "   * kappa = ",kappa
   write(*,'(a,i3)') "   * n_paths = ",npaths
   write(*,'(a,f15.8,a)') "   * T = ",kelvin," K"
!
!     convert free energies to kcal/mol
!
   pmf_max=pmf_max*2625.50d0
   pmf_min=pmf_min*2625.50d0
   
   write(*,'(a,f15.8,a,f12.8)') "   * W(xi_TS) = ",pmf_max," kJ/mol at xi= ",xi_max
   write(*,'(a,f15.8,a,f12.8)') "   * W(xi_min) = ",pmf_min, " kJ/mol at xi= ",xi_min
!
!     calculate the rate constant in atomic units 
!
   k_t=kappa*npaths*1.3806485E-23*kelvin/6.62607E-34*exp(-(pmf_max-pmf_min)/(0.00831447*kelvin))
!
!     convert the rate constant into usual units (s^(-1))
!
   write(15,*)
   write(15,*) "----------------------------------------------"
   write(15,'(a)') "  The value of the reaction rate constant is: "
   write(15,'(a,es16.8,a)') "    ",k_t,"  s^(-1)"
   write(15,*) "----------------------------------------------"
   write(*,*)
   write(*,*) "----------------------------------------------"
   write(*,'(a)') "  The value of the reaction rate constant is: "
   write(*,'(a,es16.8,a)') "    ",k_t,"  s^(-1)"
   write(*,*) "----------------------------------------------"
else 


   write(*,*) "For this type of reaction no rate constant calculation"
   write(*,*) "has been implemented so far!"
end if

return
end subroutine calc_k_t

