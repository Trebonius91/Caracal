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
   write(*,*) "------------------DYNAMICAL SETTINGS------------------"
   write(*,'(a,i5)') " * Number of beads that the ring polymer shall have: ",nbeads
   write(*,'(a,i5)') " * Number of equivalent/symmetric reaction paths: ",npaths
   if (nbeads .eq. 1) write(*,*) "  --> classical dynamics will be applied!"
   write(*,'(a,a)') " * File where the transition state is stored: ", trim(ts_file)
   write(*,'(a,i5,a)') " * The FIXED_ATOMS option was activated! In total, ",fix_num,&
            & " atoms will be hold fixed."

   write(*,'(a,f12.6)') " * Timestep (fs): ",dt_info*1000d0
   write(*,'(a,f12.6)') " * Temperature (K): ",kelvin
   if (thermo .eq. "ANDERSEN") then
      write(*,'(a)') " * The Andersen thermostat is used, the velocity distribution will be"
      write(*,'(a,i6,a)') "    reset each ",andersen_step," MD steps during dynamics."
   end if
   if (thermo .eq. "NOSE-HOOVER") then
      write(*,'(a)') " * The Nose-Hoover thermostat is used, the velocity distribution will be"
      write(*,'(a,f14.6,a)') "    coupled to heat variable with ",nose_q," coupling strength."
   end if
   if (umbr_type .eq. "BIMOLEC") then
      write(*,'(a)') " * A bimolecular substitution will be investigated."
   else if (umbr_type .eq. "BIMOL_EXCH") then
      write(*,'(a)') " * A bimolecular cyclic exchange with 2 bonds broken and two bonds "
      write(*,'(a)') "       formed  will be investigated."
   else if (umbr_type .eq. "CYCLOADD") then
      write(*,'(a)') " * A cycloaddition reaction will be investigated."
   else if (umbr_type .eq. "MERGING") then
      write(*,'(a)') " * A bimolecular merging reaction will be investigated."
   else if (umbr_type .eq. "ADDITION") then
      write(*,'(a)') " * An addition reaction with two reactants will be investigated."
   else if (umbr_type .eq. "ADDITION3") then
      write(*,'(a)') " * A cyclic addition with three reactants will be investigated."
   else if (umbr_type .eq. "ADD3_SOLV") then
      write(*,'(a)') " * A cyclic addition with three reactants, two of them preclustered,"
      write(*,'(a)') "     will be investigated."
   else if (umbr_type .eq. "ADDITION4") then
      write(*,'(a)') " * A cyclic addition with four reactants will be investigated."
   else if (umbr_type .eq. "ADD4_SOLV") then
      write(*,'(a)') " * A cyclic addition with four reactants, three of them preclustered,"
      write(*,'(a)') "     will be investigated."
   else if (umbr_type .eq. "ELIMINATION") then
      write(*,'(a)') " * An elimination starting from one reactant will be investigated."
   else if (umbr_type .eq. "REARRANGE") then
      write(*,'(a)') " * A rearrangement of one reactant will be investigated."
   else if (umbr_type .eq. "DECOM_1BOND") then
      write(*,'(a)') " * A unimolecular decomposition with one broken bond will be investigated."
   else if (umbr_type .eq. "ATOM_SHIFT") then
      write(*,'(a)') " * A single atom will be moved a certain distance in a certain " 
      write(*,'(a)') "     direction."
      write(*,'(a,i4,a)') "   - The shifted atom has the number ",shift_atom,"." 
      if (shift_coord .eq. 1) then
         write(*,'(a,f12.6,a,f12.6,a)') "   - Its x coordinate will be moved. Xi=0 is x=", &
               & shift_lo," and Xi=1 is x=",shift_hi,"."
      else if (shift_coord .eq. 2) then
         write(*,'(a,f12.6,a,f12.6,a)') "   - Its y coordinate will be moved. Xi=0 is y=", &
               & shift_lo," and Xi=1 is y=",shift_hi,"."
      else 
         write(*,'(a,f12.6,a,f12.6,a)') "   - Its z coordinate will be moved. Xi=0 is z=", &
               & shift_lo," and Xi=1 is z=",shift_hi,"."
      end if

   end if
   if (umbr_type .ne. "ATOM_SHIFT") then
      do i=1,sum_reacs
         write(*,'(a,a,a,i4,a)') "   - The ",trim(word(i))," reactant has ",n_reac(i)," atoms." 
      end do
      do i=1,form_num
         write(*,'(a,a,a,i4,a,i4)') "   - The ",trim(word(i))," forming bond is: ", &
              & bond_form(i,1),"-",bond_form(i,2)
      end do
      do i=1,break_num
         write(*,'(a,a,a,i4,a,i4)') "   - The ",trim(word(i))," breaking bond is: ", &
               & bond_break(i,1),"-",bond_break(i,2)
      end do
      write(*,'(a,f12.6,a)') "   - The asymptotic reactant distance is: ",R_inf*bohr," Angstrom"
   end if
   write(*,'(a)') " * The umbrella samplings are carried out along the reaction path"
   write(*,'(a,f12.6,a,f12.6)') "      along a range between Xi=",umbr_hi," and Xi=",umbr_lo
   write(*,'(a,i4,a)') " * There are",nint((umbr_hi-umbr_lo)/umbr_dist)," umbrella windows."
   if (k_force_indi) then
      write(*,'(a,a)') " * The umbrella force constants of the windows are read in from file ",trim(k_force_file)
   else
      write(*,'(a,f12.6,a)') " * The strength of the global umbrella force constant is: ",k_force_all," a.u."
   end if
   write(*,'(a)') "     It was multiplied by temperature as usual for umbrella samplings (k=k*T)"
   write(*,'(a,i9,a)') " * A number of ",gen_step," MD steps is used for start structure generation."
   write(*,'(a,i5,a)') " * Per umbrella window",umbr_traj," umbrella trajectories will be sampled."
   write(*,'(a,i9,a)') " * A number of ",equi_step," MD steps are used for each umbrella equilibation."
   write(*,'(a,i9,a)') " * A number of ",umbr_step," MD steps are sampled for each umbrella trajectory."
   if (pmf_method .eq. 'INTEGRATION') then
      write(*,'(a)') " * The free energy profile (PMF) will be calculated using umbrella integration."
   else if (pmf_method .eq. 'WHAM') then
      write(*,'(a)') " * The free energy profile (PMF) will be calculated using the WHAM method."
   else if (pmf_method .eq. 'ALL') then
      write(*,'(a)') " * The free energy profile (PMF) will be calculated using the WHAM method."
      write(*,'(a)') "    as well as umbrella integration. Both curves will be averaged."
   end if
   write(*,'(a,f12.6,a,f12.6)') " * The PMF is calculed between Xi=",xi_min," and Xi=",xi_max
   write(*,'(a,i7,a)') " * The PMF is calculated on ",nbins," single bins along the path."
   if (pmf_minloc .eq. 'ZERO') then
      write(*,'(a)') " * The PMF_MINLOC option was set to 'ZERO'! Therefore, the reactant asympotic"
      write(*,'(a)') "      will be taken as the lower energy for delta F."
   else if (pmf_minloc .eq. 'PMF_MIN') then 
      write(*,'(a)') " * The PMF_MINLOC option was set to 'PMF_MIN'! Therefore, the lowest PMF value"
      write(*,'(a)') "      left the TS will be taken as the lower energy for delta F."   
   end if
   write(*,'(a)') " * For the recrossing calculation, the following settings are used:"
   if (.not. recross_mpi) then
      write(*,'(a)') "   - It will be a serial calculation without MPI"
   else  
      write(*,'(a)') "   - It will be a parallel calculation utilizing MPI"
   end if   
   write(*,'(a,i9,a)') "   - At first, the parent trajectory is sampled for ",recr_equi," MD steps."
   write(*,'(a,i9,a)') "   - In total ",child_tot," child trajectories will be started."
   write(*,'(a,i6,a)') "   - Per spawn point, ",child_point," child trajectories are started."
   write(*,'(a,i9,a)') "   - Between two spawn points, the parent is sampled for ",child_interv," MD steps."
   write(*,'(a,i9,a)') "   - Each child trajectory will be evaluated for ",child_evol," MD steps."
   write(*,'(a)') " * Switches for error handling during calculation:"
   if (act_check) then
      write(*,'(a,i8,a)') "   - A maximum number of ",err_max," errors is set before canceling the calculation."
      write(*,'(a,f16.6,a)') "   - If the actual energy is",energy_tol," kJ/mol above the TS energy, give an error."
      if (.not. recross_check) then
         write(*,'(a)') "   - No error checking will be applied for recrossing calculations."
      else 
         write(*,'(a)') "   - Recrossing trajectories will be checked as well."
      end if
   else 
      write(*,'(a)') "   --> No error checking will be done at all!"
   end if
   if (mirrors) then
      write(*,*) "* The following mirror planes will be present in the system:"
      do i=1,mirror_num
         if (mirror_dims(i) .eq. 1) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : x at",mirror_pos(i)*bohr," Ang."
         else if (mirror_dims(i) .eq. 2) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : y at",mirror_pos(i)*bohr," Ang."
         else if (mirror_dims(i) .eq. 3) then
            write(*,'(a,i6,a,f12.6,a)') "     atom ",mirror_ats(i),", : z at",mirror_pos(i)*bohr," Ang."
         end if
      end do
   end if
   if (print_poly .gt. -1.d0) then
      write(*,'(a)') " * The 'PRINT_POLYMER' option was activated. Therefore one single umbrella sampling"
   end if
   if (print_gen) then
      write(*,'(a,i6,a)') " * The PRINT_GEN option was activated. Every ",print_gen_freq," MD frame "
      write(*,'(a)') "   during the structure generation part will be written to file."
      write(*,'(a)') "   (XDATCAR_gen for periodic systems, traj_gen.xyz for nonperiodic systems)"
   end if
   if (gen_test) then
      write(*,'(a)') " * The 'GEN_TEST' option was activated! Therefore, only start structure generation will be"
      write(*,'(a,i8,a)') "   done, and each ",gen_pr_frac,"'th structure will be written out."
      write(*,'(a)') "   Structures are written to gen_test_struc.xyz, Energies, s-values and z-values will"
      write(*,'(a)') "   be written to gen_test_ens.dat."
   end if
   write(*,*) "--------------------------------------------"
   write(*,*)
end if


end subroutine calc_rate_info
