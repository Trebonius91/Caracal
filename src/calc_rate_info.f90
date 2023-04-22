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
!      subroutine calc_rate_info: write out all basic information and settings 
!        for an RPMD k(T) calculation to be started into the logfile.
!        This is called only from the rpmd.x program itself
!
!      part of EVB
!
subroutine calc_rate_info(rank,dt,dtdump,dt_info,gen_step,equi_step,umbr_step,xi_min,&
          & xi_max,pmf_min,pmf_max,nbins,pmf_minloc,print_poly,ts_file,recross_mpi)
use general
use evb_mod
use debug 
implicit none
!     The current MPI rank
integer, intent(in) :: rank
!     the timestep for calculation and write out 
real(kind=8), intent(in) ::dt,dtdump
!     timestep for information of user 
real(kind=8), intent(in) ::dt_info
!     several step numbers for trajectories 
integer, intent(out) ::equi_step,gen_step,umbr_step
!     borders of XI integration area 
real(kind=8), intent(in) :: xi_min,xi_max
!     borders of umbrella windows along Xi reaction path
real(kind=8), intent(in) :: pmf_min,pmf_max
!     number of PMF integratiom bins 
integer, intent(in) :: nbins
!     method for location of PMF minimum
character(len=10), intent (in) ::pmf_minloc
!     position of printout trajectory 
real(kind=8), intent (in) :: print_poly
!     filename of the ts structure (the starting point!)
character(len=80), intent (in) ::ts_file
!     if recrossing shall be parallelized
logical, intent (in)::recross_mpi
!     characters for written out numbers 
character(len=20)::word(8)
!     loop indices 
integer::i,j
!     the pi
real(kind=8)::pi

pi=3.141592653589793238462643d0

!
!     fill dummy index for counting words
!
word=(/"first  ","second ","third  ","fourth ","fifth  ","sixth  ","seventh","eigth  "/)

if (rank .eq. 0) then
   write(15,*) "------------------SETTINGS------------------"
   write(15,'(a,i5)') " * Number of atoms in the system: ", natoms
   write(15,'(a,i5)') " * Number of beads that the ring polymer shall have: ",nbeads
   write(15,'(a,i5)') " * Number of equivalent/symmetric reaction paths: ",npaths
   if (nbeads .eq. 1) write(15,*) "  --> classical dynamics will be applied!"
   write(15,'(a,a)') " * File where the transition state is stored: ", ts_file
   if (pot_type .eq. "h3") then
      write(15,'(a)') " * Used potential energy surface: H+H2 (analytical)"
   else if (pot_type .eq. "brh2") then
      write(15,'(a)') " * Used potential energy surface: HBr+H (analytical)"
   else if (pot_type .eq. "oh3") then
      write(15,'(a)') " * Used potential energy surface: OH2+H (analytical)"
   else if (pot_type .eq. "clnh3") then
      write(15,'(a)') " * Used potential energy surface: NH2+HCl (analytical)"
   else if (pot_type .eq. "ch4h") then
      write(15,'(a)') " * Used potential energy surface: CH4+H (analytical)"
   else if (pot_type .eq. "ch4oh") then
      write(15,'(a)') " * Used potential energy surface: CH4+OH (analytical)"
   else if (pot_type .eq. "geh4oh") then
      write(15,'(a)') " * Used potential energy surface: GeH4+OH (analytical)"
   else if (pot_type .eq. "c2h7") then
      write(15,'(a)') " * Used potential energy surface: C2H6+H (analytical)"
   else 
      write(15,'(a)') " * Used potential energy surface: EVB-QMDFF"
      if (use_dq) then
         write(15,'(a)') "   - Used coupling type: EVB-dQ (distance)"
         write(15,'(a,a)') "   - Used off diagonal term: ",off_basis
      else if (dg_evb) then
         write(15,'(a)') "   - Used coupling type: DG-EVB (distributed gaussian)"
         if (dg_mode .eq. 1) then
            write(15,'(a)') "   - Used DG-EVB mode: 1 (energies)"
         else if (dg_mode .eq. 2) then
            write(15,'(a)') "   - Used DG-EVB mode: 2 (energies + gradients)"
         else
            write(15,'(a)') "   - Used DG-EVB mode: 3 (energies + gradients + hessians)"
         end if
         write(15,'(a,i5)') "   - Number of DG-EVB points: ", dg_evb_points
         write(15,'(a,i5)') "   - Number of alpha coefficients: ", dg_evb_points+add_alph
         write(15,'(a)') "   - Used internal coordinates:"
         do i=1,nat6
            if (coord_def(i,1) .eq. 1) then
               write(15,'(a,i5,i5)') "       bond: ", coord_def(i,2:3)
            else if (coord_def(i,1) .eq. 2) then
               write(15,'(a,i5,i5,i5)') "       bend angle: ", coord_def(i,2:4)
            else if (coord_def(i,1) .eq. 3) then
               write(15,'(a,i5,i5,i5,i5)') "       dihedral angle: ", coord_def(i,2:5)
            else if (coord_def(i,1) .eq. 4) then
               write(15,'(a,i5,i5,i5,i5)') "       out of plane angle: ", coord_def(i,2:5)
            end if
         end do
      else if (treq) then
         write(15,'(a)') "   - Used coupling type: RP-EVB (reaction path)"
         write(15,'(a,i5)') "   - Number of RP-EVB points: ", rp_evb_points
         write(15,'(a,f10.5)') "   - Position of TS (s): ",s_ts
         write(15,'(a,f7.5,a,f10.5)') "   - Borders of left RP/interpolation transition (s): ",&
                    & trans_l_lo," - ",trans_l_hi
         write(15,'(a,f7.5,a,f10.5)') "   - Borders of right RP/interpolation transition (s): ",&
                    & trans_r_lo," - ",trans_r_hi
         write(15,'(a,f15.7)') "   - Exponential damping coefficient: ",pre_exp
         write(15,'(a,f15.7,a,f15.7)') "   - Borders of RP-EVB/QMDFF (lower/upper): ",&
                    & pi/pareta," - ",1-pi/pareta
         write(15,'(a,i6,a)') "   - Used internal coordinates (number: ",nat6,") :"
         do i=1,nat6
            if (coord_def(i,1) .eq. 1) then
               write(15,'(a,i5,i5)') "       bond: ", coord_def(i,2:3)
            else if (coord_def(i,1) .eq. 2) then
               write(15,'(a,i5,i5,i5)') "       bend angle: ", coord_def(i,2:4)
            else if (coord_def(i,1) .eq. 3) then
               write(15,'(a,i5,i5,i5,i5)') "       dihedral angle: ", coord_def(i,2:5)
            else if (coord_def(i,1) .eq. 4) then
               write(15,'(a,i5,i5,i5,i5)') "       out of plane angle: ", coord_def(i,2:5)
            end if
         end do
      else 
         write(15,'(a)') "   - Used coupling type: EVB-dE (energy)"
         write(15,'(a,a)') "   - Used off diagonal term: ",off_basis
      end if
   end if
   write(15,'(a,i5,a)') " * The FIXED_ATOMS option was activated! In total, ",fix_num,&
            & " atoms will be hold fixed."

   write(15,'(a,f12.6)') " * Timestep (fs): ",dt_info*1000d0
   write(15,'(a,f12.6)') " * Temperature (K): ",kelvin
   if (thermo .eq. "ANDERSEN") then
      write(15,'(a)') " * The Andersen thermostat is used, the velocity distribution will be"
      write(15,'(a,i6,a)') "    reset each ",andersen_step," MD steps during dynamics."
   end if
   if (thermo .eq. "NOSE-HOOVER") then
      write(15,'(a)') " * The Nose-Hoover thermostat is used, the velocity distribution will be"
      write(15,'(a,f14.6,a)') "    coupled to heat variable with ",nose_q," coupling strength."
   end if
   if (umbr_type .eq. "BIMOLEC") then
      write(15,'(a)') " * A bimolecular reaction will be investigated."
   else if (umbr_type .eq. "CYCLOADD") then
      write(15,'(a)') " * A cycloaddition reaction will be investigated."
   else if (umbr_type .eq. "MERGING") then
      write(15,'(a)') " * A bimolecular merging reaction will be investigated."
   else if (umbr_type .eq. "ADDITION") then
      write(15,'(a)') " * An addition reaction with two reactants will be investigated."
   else if (umbr_type .eq. "ADDITION3") then
      write(15,'(a)') " * A cyclic addition with three reactants will be investigated."
   else if (umbr_type .eq. "ADD3_SOLV") then
      write(15,'(a)') " * A cyclic addition with three reactants, two of them preclustered,"
      write(15,'(a)') "     will be investigated."
   else if (umbr_type .eq. "ADDITION4") then
      write(15,'(a)') " * A cyclic addition with four reactants will be investigated."
   else if (umbr_type .eq. "ADD4_SOLV") then
      write(15,'(a)') " * A cyclic addition with four reactants, three of them preclustered,"
      write(15,'(a)') "     will be investigated."
   else if (umbr_type .eq. "ELIMINATION") then
      write(15,'(a)') " * An elimination starting from one reactant will be investigated."
   else if (umbr_type .eq. "REARRANGE") then
      write(15,'(a)') " * A rearrangement of one reactant will be investigated."
   else if (umbr_type .eq. "DECOM_1BOND") then
      write(15,'(a)') " * A unimolecular decomposition with one broken bond will be investigated."
   else if (umbr_type .eq. "ATOM_SHIFT") then
      write(15,'(a)') " * A single atom will be moved a certain distance in a certain " 
      write(15,'(a)') "     direction."
      write(15,'(a,i4,a)') "   - The shifted atom has the number ",shift_atom,"." 
      if (shift_coord .eq. 1) then
         write(15,'(a,f12.6,a,f12.6,a)') "   - Its x coordinate will be moved. Xi=0 is x=", &
               & shift_lo," and Xi=1 is x=",shift_hi,"."
      else if (shift_coord .eq. 2) then
         write(15,'(a,f12.6,a,f12.6,a)') "   - Its y coordinate will be moved. Xi=0 is y=", &
               & shift_lo," and Xi=1 is y=",shift_hi,"."
      else 
         write(15,'(a,f12.6,a,f12.6,a)') "   - Its z coordinate will be moved. Xi=0 is z=", &
               & shift_lo," and Xi=1 is z=",shift_hi,"."
      end if

   end if
   if (umbr_type .ne. "ATOM_SHIFT") then
      do i=1,sum_reacs
         write(15,'(a,a,a,i4,a)') "   - The ",trim(word(i))," reactant has ",n_reac(i)," atoms." 
      end do
      do i=1,form_num
         write(15,'(a,a,a,i4,a,i4)') "   - The ",trim(word(i))," forming bond is: ", &
              & bond_form(i,1),"-",bond_form(i,2)
      end do
      do i=1,break_num
         write(15,'(a,a,a,i4,a,i4)') "   - The ",trim(word(i))," breaking bond is: ", &
               & bond_break(i,1),"-",bond_break(i,2)
      end do
      write(15,'(a,f12.6,a)') "   - The asymptotic reactant distance is: ",R_inf*bohr," Angstrom"
   end if
   write(15,'(a)') " * The umbrella samplings are carried out along the reaction path"
   write(15,'(a,f12.6,a,f12.6)') "      along a range between Xi=",umbr_hi," and Xi=",umbr_lo
   write(15,'(a,i4,a)') " * There are",nint((umbr_hi-umbr_lo)/umbr_dist)," umbrella windows."
   write(15,'(a,f12.6,a)') " * The strength of the umbrella force constant is: ",k_force," a.u."
   write(15,'(a)') "     It was multiplied by temperature as usual for umbrella samplings (k=k*T)"
   write(15,'(a,i9,a)') " * A number of ",gen_step," MD steps is used for start structure generation."
   write(15,'(a,i5,a)') " * Per umbrella window",umbr_traj," umbrella trajectories will be sampled."
   write(15,'(a,i9,a)') " * A number of ",equi_step," MD steps are used for each umbrella equilibation."
   write(15,'(a,i9,a)') " * A number of ",umbr_step," MD steps are sampled for each umbrella trajectory."
   if (pmf_method .eq. 'INTEGRATION') then
      write(15,'(a)') " * The free energy profile (PMF) will be calculated using umbrella integration."
   else if (pmf_method .eq. 'WHAM') then
      write(15,'(a)') " * The free energy profile (PMF) will be calculated using the WHAM method."
   else if (pmf_method .eq. 'ALL') then
      write(15,'(a)') " * The free energy profile (PMF) will be calculated using the WHAM method."
      write(15,'(a)') "    as well as umbrella integration. Both curves will be averaged."
   end if
   write(15,'(a,f12.6,a,f12.6)') " * The PMF is calculed between Xi=",xi_min," and Xi=",xi_max
   write(15,'(a,i7,a)') " * The PMF is calculated on ",nbins," single bins along the path."
   if (pmf_minloc .eq. 'ZERO') then
      write(15,'(a)') "  * The PMF_MINLOC option was set to 'ZERO'! Therefore, the reactant asympotic"
      write(15,'(a)') "      will be taken as the lower energy for delta F."
   else if (pmf_minloc .eq. 'PMF_MIN') then 
      write(15,'(a)') "  * The PMF_MINLOC option was set to 'PMF_MIN'! Therefore, the lowest PMF value"
      write(15,'(a)') "      left the TS will be taken as the lower energy for delta F."   
   end if
   write(15,'(a)') " * For the recrossing calculation, the following settings are used:"
   if (.not. recross_mpi) then
      write(15,'(a)') "   - It will be a serial calculation without MPI"
   else  
      write(15,'(a)') "   - It will be a parallel calculation utilizing MPI"
   end if
   write(15,'(a,i9,a)') "   - At first, the parent trajectory is sampled for ",recr_equi," MD steps."
   write(15,'(a,i9,a)') "   - In total ",child_tot," child trajectories will be started."
   write(15,'(a,i6,a)') "   - Per spawn point, ",child_point," child trajectories are started."
   write(15,'(a,i9,a)') "   - Between two spawn points, the parent is sampled for ",child_interv," MD steps."
   write(15,'(a,i9,a)') "   - Each child trajectory will be evaluated for ",child_evol," MD steps."
   write(15,'(a)') " * Switches for error handling during calculation:"
   if (act_check) then
      write(15,'(a,i8,a)') "   - A maximum number of ",err_max," errors is set before canceling the calculation."
      write(15,'(a,f16.6,a)') "   - If the actual energy is",energy_tol,"kJ/mol above the TS energy, give an error."
      if (.not. recross_check) then
         write(15,'(a)') "   - No error checking will be applied for recrossing calculations."
      else 
         write(15,'(a)') "   - Recrossing trajectories will be checked as well."
      end if
   else 
      write(15,'(a)') "   --> No error checking will be done at all!"
   end if
   if (print_poly .gt. -1.d0) then
      write(15,'(a)') " * The 'PRINT_POLYMER' option was activated. Therefore one single umbrella sampling"

   end if
   if (gen_test) then
      write(15,'(a)') " * The 'GEN_TEST' option was activated! Therefore, only start structure generation will be"
      write(15,'(a,i8,a)') "   done, and each ",gen_pr_frac,"'th structure will be written out."
      write(15,'(a)') "   Structures are written to gen_test_struc.xyz, Energies, s-values and z-values will"
      write(15,'(a)') "   be written to gen_test_ens.dat."
   end if
   write(15,*) "--------------------------------------------"
   write(15,*)
end if


end subroutine calc_rate_info
