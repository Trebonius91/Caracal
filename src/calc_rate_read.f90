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
!     subroutine calc_rate_read: read in several settings from caracal.key file 
!      for the rpmd.x program
!
!     part of EVB
!
subroutine calc_rate_read(rank,dt,dtdump,dt_info,gen_step,equi_step,umbr_step,xi_min,&
          & xi_max,pmf_min,pmf_max,nbins,pmf_minloc,print_poly,ts_file,recross_mpi)
use general
use evb_mod
use debug 
use qmdff
implicit none
!     The current MPI rank
integer, intent(out) :: rank 
!     the timestep for calculation and write out 
real(kind=8), intent(out) ::dt,dtdump
!     timestep for information of user 
real(kind=8), intent(out) ::dt_info
!     several step numbers for trajectories 
integer, intent(out) ::equi_step,gen_step,umbr_step
!     borders of XI integration area 
real(kind=8), intent(out) :: xi_min,xi_max
!     borders of umbrella windows along Xi reaction path
real(kind=8), intent(out) :: pmf_min,pmf_max
!     number of PMF integratiom bins 
integer, intent(out) :: nbins
!     method for location of PMF minimum
character(len=10), intent (out) ::pmf_minloc
!     position of printout trajectory 
real(kind=8), intent (out) :: print_poly
!     filename of the ts structure (the starting point!)
character(len=80), intent (out) ::ts_file
!     filename of the reactants structure (for unimolecular)
character(len=80)::reactants_file
character(len=50)::fix_file ! file with numbers of fixed atoms
!     if recrossing shall be parallelized
logical, intent (out)::recross_mpi
!     Loop indices 
integer::i,j,k,l,p,m
!     Needed characters for read in
character(len=20) keyword
character(len=60)::names,names2
character(len=140) record
character(len=140) string
!     local arrays for breaking/forming bonds 
character(len=20),dimension(:),allocatable::input_form,input_break
!     for read in of reactant atoms
character(len=1)::at_index(200)
!     auxiliary characters for read in of bonds 
character(len=1)::act_char
character(len=20)::act_bond
character(len=20)::act_number
logical::active,exist,section
integer::next
!     status of actual read in
integer::readstat,istat

!     Initialize several parameters 
!
!     number of steps per umbrella structure generation trajectory
gen_step = 0
!     number of steps per umbrella equilibration trajectory
equi_step = 0
!     number of steps per umbrella sampling trajectory
umbr_step = 0
!     number of umbrella trajectories per sampling window
umbr_traj = 0
!     temperature for the whole calculations (Kelvin)
kelvin = 0.0d0
!     the nose-hoover damping
nose_q = 100.0
!     the andersen step, after how many MD steps the velocities are rescaled
andersen_step = 60
!     timestep for the trajectories (fs)
dt = -1.d0
!     if a trajectory shall be printed out, frequency for printout structures 
dtdump = 1.d0
!     number of RPMD beads 
nbeads = 0
!     number of equivalent reaction pahts (k(T) multiplier)
npaths = 1
!     if sampling trajectories shall be checked for errors 
act_check = .true.
!     if recrossing trajectories shall be checked for errors 
recross_check = .true.
!     If the loose check option shall be activated (only development!)
loose_check = .false.
!     The umbrella force constant
k_force = -1.d0
!     The lowest umbrella window (with respect to path progress Xi)
umbr_lo = -5.d0
!     The highest umbrella window (with respect to path progress Xi)
umbr_hi = -5.d0
!     The distance between two umbrella windows 
umbr_dist = -1.d0
!     The reaction class with respect to force constant etc.
umbr_type = "none"
!     The asymptotic reactants separation for umbrella samplings
r_inf = 0.d0
!     The file in which the startign TS trajectory shall be stored 
ts_file = " "
!     The lowest reaction path Xi value for PMF calculation
xi_min = 1.d0
!     The highest reaction path Xi value for PMF calculation
xi_max = -1.d0
!     The number of bins for numerical PMF integration
nbins = 0
!     The method for the calculation of PMF surfaces 
pmf_method = "NONE"
!     If the recrossing calculations shall be paralellized 
recross_mpi = .false.
!     Number of steps for initial recrossing equilibration trajectory
recr_equi = 0
!     Total number of child trajectories 
child_tot = 0
!     timesteps between two child spawning events
child_interv = 0
!     Number of child trajectories per spawning interval
child_point = 0
!     timesteps for a single child trajectory
child_evol = 0
!     Method for location of minimum PMF value on path (asympototic or global min.)
pmf_minloc = "ZERO"
!     If a single umbrella trajectory shall be printed to illustrate RPMD
print_poly = -5.d0
!     The energy tolerance in kJ/mol with respect to the TS energy for errors 
energy_tol = 250d0
!     If only the structure generation shall be done for test reasons
gen_test = .false.
!     If the debugging mode shall be set on (write out every single structure)
do_debug = .false.
!     If atoms of the system shall be fixed during calculation
fix_atoms=.false.
!     No NVE ensemble will be used
nve = .false.
!     No NpT ensemble will be used 
npt = .false.
!     No periodic calculation
periodic = .false.

!    
!     Start big loop over all keywords and read them in
!
write(15,*)
write(15,*) "Read in all keywords for the program..."
write(15,*)

! -------------------------------------------------------------------
!     A) First, read in settings that do not belong to a section
! -------------------------------------------------------------------
write(15,*) "A) General parameters..."
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
!
!     reset the read status to successful (check it for each keyword)
!
   readstat= 0
!
!     If a traditional force field shall be calculated, activate global flag!
!     Here also read in a start structure like in dynamic.x to determine the 
!     trad. FF coordinates etc.
!

   if (keyword(1:20) .eq. 'GENERATE_STEPS ') then
      read(record,*,iostat=readstat) names,gen_step
   else if (keyword(1:20) .eq. 'EQUILIBR_STEPS ') then
   else if (keyword(1:20) .eq. 'TS_STRUC ') then
      read(record,*,iostat=readstat) names,ts_file
      if (readstat .ne. 0) then
         write(*,*) "Correct format: TS_STUC [filename]"
         call fatal
      end if
   else if (keyword(1:20) .eq. 'RPMD_BEADS ') then
      read(record,*,iostat=readstat) names,nbeads
      if (readstat .ne. 0) then
         write(*,*) "Correct format: RPMD_BEADS [Number of RPMD beads]"
         call fatal
      end if
   else if (keyword(1:11) .eq. 'DELTAT ') then
      read(record,*,iostat=readstat) names,dt
      if (readstat .ne. 0) then
         write(*,*) "Correct format: DELTAT [MD time step (fs)]"
         call fatal
      end if
   else if (keyword(1:16) .eq. 'PRINT_POLYMER ') then
      read(record,*,iostat=readstat) names,print_poly
      if (readstat .ne. 0) then
         write(*,*) "Correct format: PRINT_POLYMER [Xi-value]"
         call fatal
      end if
   else if (keyword(1:20) .eq. 'FIX_ATOMS ') then
      fix_atoms=.true.
      read(record,*) names,fix_file
      if (readstat .ne. 0) then
         write(*,*) "Correct format: FIX_ATOMS [filename with indices]"
         call fatal
      end if
   else if (keyword(1:11) .eq. 'TDUMP ') then
      read(record,*) names,dtdump
      if (readstat .ne. 0) then
         write(*,*) "Correct format: TDUMP [time interval (fs)]"
         call fatal
      end if
!  
!     calculate the frequency of print outs
!     
      iwrite = nint(dtdump/dt_info*0.001d0)
   else if (keyword(1:16) .eq. 'RPMD_EN_TOL ') then
      read(record,*,iostat=readstat) names,energy_tol
      if (readstat .ne. 0) then
         write(*,*) "Correct format: RPMD_EN_TOL [energy tolerance (kJ/mol)]"
         call fatal
      end if
   else if (keyword(1:20) .eq. 'GEN_TEST ') then
      gen_test = .true.
   else if (keyword(1:20) .eq. 'MAX_ERROR ') then
      read(record,*,iostat=readstat) names,err_max
      if (readstat .ne. 0) then
         write(*,*) "Correct format: MAX_ERROR [Number of errors]"
         call fatal
      end if
   else if (keyword(1:20) .eq. 'NO_CHECK ') then
      act_check=.false.
   else if (keyword(1:16) .eq. 'DEBUG ') then
      do_debug=.true.
      debug_unit=112
   end if
end do
write(15,*) "...done!"

! -------------------------------------------------------------------
!     B) SECOND, read in the NVT section!
! -------------------------------------------------------------------
section=.false.
write(15,*) "B) The NVT section..."
do i = 1, nkey
   next = 1 
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'NVT {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'NVT{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!    The simulation temperature
         if (keyword(1:11) .eq. 'TEMP ') then
            read(record,*,iostat=readstat) names,kelvin
            if (readstat .ne. 0) then
               write(*,*) "Correct format: TEMP [temperature (K)]"
               call fatal
            end if
!    The thermostat 
         else if (keyword(1:20) .eq. 'THERMOSTAT ') then
            read(record,*,iostat=readstat) names,thermo
            if (readstat .ne. 0) then
               write(*,*) "Correct format: THERMOSTAT [thermostat]"
               call fatal
            end if
            call upcase(thermo)
!    How often the Andersen thermostat velocity rescaling is applied
         else if (keyword(1:16) .eq. 'ANDERSEN_STEP ') then
            read(record,*,iostat=readstat) names,andersen_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: ANDERSEN_STEP [No. of steps]"
               call fatal
            end if
!    The damping factor for the Nose-Hoover thermostat
         else if (keyword(1:13) .eq. 'NOSE_DAMP ') then
            read(record,*,iostat=readstat) names,nose_q
            if (readstat .ne. 0) then
               write(*,*) "Correct format: NOSE_DAMP [integer (>1)]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The NVT section has no second delimiter! (})"
            call fatal
         end if
      end do
      exit
   end if
end do
if (.not. section) then
   write(*,*) "Please add a NVT section with parameters to the keyfile!"
   call fatal
end if
write(15,*) "... done!"
section = .false.

! -------------------------------------------------------------------
!     C) THIRD, read in the MECHA section!
! -------------------------------------------------------------------

write(15,*) "C) The MECHA section..."
section=.false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'MECHA {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'MECHA{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!    First, read in only the TYPE, and decide based on this, which parameters must 
!      be read in else 
         if (keyword(1:20) .eq. 'TYPE ') then
            read(record,*,iostat=readstat) names,umbr_type
            call upcase(umbr_type)
            if (readstat .ne. 0) then
               write(*,*) "Correct format: TYPE [Identifier for mechanism]"
               call fatal
            end if
         end if
      end do
   end if
end do
if (.not. section) then
   write(*,*) "Please add a MECHA section with parameters to the keyfile!"
   call fatal
end if

!
!    Read in the type of reaction and define the number of forming and breaking bonds
!
call upcase(umbr_type)
if (umbr_type .eq. "UNIMOLEC") then
   form_num=1
   break_num=1
   sum_reacs=1
else if (umbr_type .eq. "BIMOLEC") then
   form_num=1
   break_num=1
   sum_reacs=2
else if (umbr_type .eq. "CYCLOADD") then
   form_num=2
   break_num=0
   sum_reacs=2
else if (umbr_type .eq. "MERGING") then
   form_num=1
   break_num=0
   sum_reacs=2
else if (umbr_type .eq. "ADDITION") then
   form_num=2
   break_num=1
   sum_reacs=2
else if (umbr_type .eq. "ADDITION3") then
   form_num=3
   break_num=2
   sum_reacs=3
else if (umbr_type .eq. "ADD3_SOLV") then
   form_num=3
   break_num=2
   sum_reacs=2
else if (umbr_type .eq. "ADDITION4") then
   form_num=4
   break_num=3
   sum_reacs=4
else if (umbr_type .eq. "ADD4_SOLV") then
   form_num=4
   break_num=3
   sum_reacs=2
else if (umbr_type .eq. "ELIMINATION") then
   form_num=1
   break_num=2
   sum_reacs=1
else if (umbr_type .eq. "CYCLOREVER") then
   form_num=0
   break_num=2
   sum_reacs=1
else if (umbr_type .eq. "REARRANGE") then
   form_num=1
   break_num=1
   sum_reacs=1
else if (umbr_type .eq. "DECOM_1BOND") then
   form_num=0
   break_num=1
   sum_reacs=1
else if (umbr_type .eq. "ATOM_SHIFT") then
else
   if (rank .eq. 0) then
      write(*,*) "This umbrella coordinate set is currently not supported!"
      write(*,*) "Use one of the supported sets (look into the manual)"
   end if
   call fatal
end if

!
!    Allocate arrays for reaction type parametrization
!
allocate(at_reac(sum_reacs,200),n_reac(sum_reacs),mass_reac(sum_reacs))
at_reac=0
n_reac=0

!
!     allocate arrays for breaking and forming bonds as well as for its 
!     reading in
!

if (umbr_type .ne. "ATOM_SHIFT") then
   allocate(bond_form(form_num,2))
   allocate(bond_break(break_num,2))
   allocate(input_form(form_num))
   allocate(input_break(break_num))
end if

do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'MECHA {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'MECHA{') then
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!
!     Now read in all other informaton of the MECHA section
!
         if (keyword(1:20) .eq. 'REACTANT1 ') then
            read(record,*) names
            k=1
!
!     Read in the atom indices of each single reactant: quite complicated in fortran..
!
            act_number=" "
            do l=1,120
               at_index(l)=record(l:l)
               if (at_index(l) .ne. " ") then
                  if (len(trim(act_number)) .lt. 1) then
                     act_number=at_index(l)
                  else
                     act_number=trim(act_number) // at_index(l)
                  end if
               else
                  if (len(trim(act_number)) .ge. 1) then
                     read(act_number,*,iostat=istat) at_reac(1,k)
                     if (istat .eq. 0) then
                        k=k+1
                     else
                        at_reac(1,k)=0
                     end if
                  end if
                  act_number=" "
               end if
            end do
         else if (keyword(1:20) .eq. 'REACTANT2 ') then
            read(record,*) names
            k=1
!
!     Do the same for the second reactant..
!
            act_number=" "
            do l=1,120
               at_index(l)=record(l:l)
               if (at_index(l) .ne. " ") then
                  if (len(trim(act_number)) .lt. 1) then
                     act_number=at_index(l)
                  else
                     act_number=trim(act_number) // at_index(l)
                  end if
               else
                  if (len(trim(act_number)) .ge. 1) then
                     read(act_number,*,iostat=istat) at_reac(2,k)
                     if (istat .eq. 0) then
                        k=k+1
                     else
                        at_reac(2,k)=0
                     end if
                  end if
                  act_number=" "
               end if
            end do
         else if (keyword(1:20) .eq. 'REACTANT3 ') then
            if (sum_reacs .gt. 2) then
               read(record,*) names
               k=1
!
!     Do the same for the third reactant..
!
               act_number=" "
               do l=1,120
                  at_index(l)=record(l:l)
                  if (at_index(l) .ne. " ") then
                     if (len(trim(act_number)) .lt. 1) then
                        act_number=at_index(l)
                     else
                        act_number=trim(act_number) // at_index(l)
                     end if
                  else
                     if (len(trim(act_number)) .ge. 1) then
                        read(act_number,*,iostat=istat) at_reac(3,k)
                        if (istat .eq. 0) then
                           k=k+1
                        else
                           at_reac(3,k)=0
                        end if
                     end if
                     act_number=" "
                  end if
               end do
            end if
         else if (keyword(1:20) .eq. 'REACTANT4 ') then
            if (sum_reacs .eq. 4) then
               read(record,*) names
               k=1
!     
!     Do the same for the fourth reactant..
!
               act_number=" "
               do l=1,120
                  at_index(l)=record(l:l)
                  if (at_index(l) .ne. " ") then
                     if (len(trim(act_number)) .lt. 1) then
                        act_number=at_index(l)
                     else
                        act_number=trim(act_number) // at_index(l)
                     end if
                  else
                     if (len(trim(act_number)) .ge. 1) then
                        read(act_number,*,iostat=istat) at_reac(4,k)
                        if (istat .eq. 0) then
                           k=k+1
                        else
                           at_reac(4,k)=0
                        end if
                     end if
                     act_number=" "
                  end if
               end do
            end if
!    The infinite reactant separation for bimolecular reactions
         else if (keyword(1:15) .eq. 'DIST_INF ') then
            read(record,*,iostat=readstat) names,r_inf
            if (readstat .ne. 0) then
               write(*,*) "Correct format: DIST_INF [distance (Ang.)]"
               call fatal
            end if
!    Convert the distance from Angstroms to bohr
            r_inf = r_inf/bohr
!    The number of equivalent reaction paths / symmetry number
         else if (keyword(1:20) .eq. 'N_PATHS ') then
            read(record,*,iostat=readstat) names,npaths
            if (readstat .ne. 0) then
               write(*,*) "Correct format: N_PATHS [number]"
               call fatal
            end if
!    The list of formed bonds
         else if (keyword(1:20) .eq. 'BOND_FORM ') then
!
!     Read all bonds with one command
!
            read(record,*) names,input_form
            do m=1,form_num
               l=1
               act_bond=input_form(m)
!
!     Analyze each substring, extract the numbers of the two participating atoms
!
               act_number=" "
               do k=1,20
                  act_char=act_bond(k:k)
                  if (act_char .ne. "-" .and. act_char .ne. " ") then
                     act_number=trim(act_number) // act_char
                  else
                     if (len(trim(act_number)) .ne. 0) then
                        read(act_number,*) bond_form(m,l)
                        l=l+1
                        act_number=" "
                     end if
                  end if
!
!     if both bonded parters are found, exit the inner loop
!
                  if (l .gt. 2) exit
               end do
            end do
!    The list of broken bonds
         else if (keyword(1:20) .eq. 'BOND_BREAK ') then
            read(record,*) names,input_break
            do m=1,break_num
               l=1
               act_bond=input_break(m)
!
!     Analyze each substring, extract the numbers of the two participating atoms
!
               act_number=" "
               do k=1,20
                  act_char=act_bond(k:k)
                  if (act_char .ne. "-" .and. act_char .ne. " ") then
                     act_number=trim(act_number) // act_char
                  else
                     if (len(trim(act_number)) .ne. 0) then
                        read(act_number,*) bond_break(m,l)
                        l=l+1
                        act_number=" "
                     end if
                  end if
!
!     if both bonded parters are found, exit the inner loop
!
                  if (l .gt. 2) exit
               end do
            end do
!
!     Fpr unimolecular reactions, read in the structure of the reactant as well
!
!
!     For unimoleculae reactions: read in structure of reactant molecule!
!
         else if (keyword(1:20) .eq. 'REACTANTS_STRUC ') then
            if (sum_reacs .eq. 1) then
              ! fix_atoms=.true.
               read(record,*) names,reactants_file

               open(unit=31,file=reactants_file,status="old",iostat=readstat)
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The file ",reactants_file," with the reactants structure &
                                       &  could not been found!"
                     call fatal
                  end if
               end if
               allocate(reac_ref(3,natoms))
               read(31,*) ; read(31,*)
               do l=1,natoms
                  read(31,*) names,reac_ref(:,l)
                  end do
               close(31)
!
!     convert reactants structure to bohr
!
               reac_ref=reac_ref/bohr
            end if

!
!     If the special ATOM_SHIFT type is used, read in all needed info and go 
!     directly to the end of this subsection
!
!     Index of the shifted atom
         else if (keyword(1:20) .eq. 'SHIFT_ATOM') then 
            if (umbr_type .eq. "ATOM_SHIFT") then
               read(record,*) names,shift_atom
               if (shift_atom .lt. 1 .or. shift_atom .gt. natoms) then
                  write(*,*) "The atom listed in SHIFT_ATOM is not in the system!"
               end if
            end if
!     Which cartesian coordinate (X,Y,Z) of the atom shift
         else if (keyword(1:20) .eq. 'SHIFT_COORD') then
            if (umbr_type .eq. "ATOM_SHIFT") then
               read(record,*) names,names2 
               call upcase(names2)
               if (trim(names2) .eq. "X") then
                  shift_coord=1
               else if (trim(names2) .eq. "Y") then
                  shift_coord=2
               else if (trim(names2) .eq. "Z") then
                  shift_coord=3
               else
                  write(*,*) "No valid coordinate (x, y or z) is given for SHIFT_COORD!"
                  call fatal
               end if
            end if
!     The interval of the atom shift
         else if (keyword(1:20) .eq. 'SHIFT_INTERV') then
            if (umbr_type .eq. "ATOM_SHIFT") then
               read(record,*) names,shift_lo,shift_hi
!
!     convert coordinate borders to bohr
!
               shift_lo=shift_lo/bohr
               shift_hi=shift_hi/bohr
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The MECHA section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do
write(15,*) "... done!"

! -------------------------------------------------------------------
!     D) FOURTH, read in the UMBRELLA section!
! -------------------------------------------------------------------

section=.false.
write(15,*) "D) The UMBRELLA section..."
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'UMBRELLA {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'UMBRELLA{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!
!     Now read in all other informaton of the MECHA section
!
!     The strength of the umbrella force constant
         if (keyword(1:20) .eq. 'BIAS ') then
            read(record,*,iostat=readstat) names,k_force
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BIAS [force constant (a.u.)]"
               call fatal
            end if
!     The interval of the reaction coordinate Xi where sampling shall happen
         else if (keyword(1:20) .eq. 'BONDS ') then
            read(record,*,iostat=readstat) names,umbr_lo,umbr_hi
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BONDS [lower Xi-val, upper Xi-val]"
               call fatal
            end if
!     Distance between two umbrella windows
         else if (keyword(1:20) .eq. 'DIST ') then
            read(record,*,iostat=readstat) names,umbr_dist
            if (readstat .ne. 0) then
               write(*,*) "Correct format: DIST [Xi-interval between two windows]"
               call fatal
            end if
!     Number of MD steps per window for structure generation
         else if (keyword(1:20) .eq. 'GEN_STEPS ') then
            read(record,*,iostat=readstat) names,gen_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: GEN_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of MD steps per window for umbrella equilibration
         else if (keyword(1:20) .eq. 'EQUI_STEPS ') then
            read(record,*,iostat=readstat) names,equi_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: EQUI_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of MD steps per window for umbrella sampling
         else if (keyword(1:20) .eq. 'SAMPLE_STEPS ') then
            read(record,*,iostat=readstat) names,umbr_step
            if (readstat .ne. 0) then
               write(*,*) "Correct format: SAMPLE_STEPS [Number of MD steps]"
               call fatal
            end if
!     Number of sampling trajectories per umbrella window
         else if (keyword(1:20) .eq. 'SAMPLE_TRAJS ') then
            read(record,*,iostat=readstat) names,umbr_traj
            if (readstat .ne. 0) then
               write(*,*) "Correct format: SAMPLE_TRAJS [Number of trajectories]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The UMBRELLA section has no second delimiter! (})"
            call fatal
         end if

      end do
   end if
end do
if (.not. section) then
   write(*,*) "Please add a UMBRELLA section with parameters to the keyfile!"
   call fatal
end if

write(15,*) "... done!"

! -------------------------------------------------------------------
!     D) FIFTH, read in the PMF section!
! -------------------------------------------------------------------

section=.false.
write(15,*) "D) The PMF section..."
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'PMF {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'PMF{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!     Range along the reaction path variable for integration
         if (keyword(1:20) .eq. 'XI_RANGE ') then
            read(record,*,iostat=readstat) names,xi_min,xi_max
            if (readstat .ne. 0) then
               write(*,*) "Correct format: XI_RANGE [lower Xi-val, upper Xi-val]"
               call fatal
            end if
!     Number of bins in which the PMF shall be calculated
         else if (keyword(1:20) .eq. 'BINS ') then
            read(record,*,iostat=readstat) names,nbins
            if (readstat .ne. 0) then
               write(*,*) "Correct format: BINS [number of PMF calculation bins]"
               call fatal
            end if
!     Which PMF calculation method shall be used
         else if (keyword(1:20) .eq. 'METHOD ') then
            read(record,*,iostat=readstat) names,pmf_method
            if (readstat .ne. 0) then
               write(*,*) "Correct format: METHOD [method identifiert]"
               call fatal
            end if
!     How the minimum along the PMF surface shall be located
         else if (keyword(1:20) .eq. 'MINLOC ') then
            read(record,*,iostat=readstat) names,pmf_minloc
            if (readstat .ne. 0) then
               write(*,*) "Correct format: MINLOC [method identifiert]"
               call fatal
            end if
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The PMF section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do
if (.not. section) then
   write(*,*) "Please add a PMF section with parameters to the keyfile!"
   call fatal
end if
write(15,*) "... done!"

! -------------------------------------------------------------------
!     E) SIXTH, read in the RECROSS section!
! -------------------------------------------------------------------

section=.false.
skip_recross=.false.
write(15,*) "E) The RECROSS section..."
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:15))) .eq. 'RECROSS {' .or. trim(adjustl(record(1:15))) &
              &  .eq. 'RECROSS{') then
      section=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!     Number of MD steps for equilibration trajectory
         if (keyword(1:20) .eq. 'EQUI_STEPS ') then
            read(record,*,iostat=readstat) names,recr_equi
            if (readstat .ne. 0) then
               write(*,*) "Correct format: EQUI_STEPS [Number of MD steps]"
               call fatal
            end if
!     Total number of child trajectories
         else if (keyword(1:20) .eq. 'CHILD_TOTAL ') then
            read(record,*,iostat=readstat) names,child_tot
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_TOTAL [Number of child trajs.]"
               call fatal
            end if
            if (child_tot .eq. 0) then
               if (rank .eq. 0) then
                  write(*,*) "Zero child trajectories shall be calculated, therefore,"
                  write(*,*) " the recrossing calculation will be skipped and a transmission "
                  write(*,*) " of 1.0 will be assumed."
               end if
               skip_recross=.true.
               exit
            end if
!     Number of parent MD steps between child spawnings
         else if (keyword(1:20) .eq. 'CHILD_INTERVAL ') then
            read(record,*,iostat=readstat) names,child_interv
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_INTERVAL [MD steps between childs]"
               call fatal
            end if
!     Number of children trajectories per spawning
         else if (keyword(1:20) .eq. 'CHILD_PERPOINT ') then
            read(record,*,iostat=readstat) names,child_point
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_PERPOINT [Number of children p.p.]"
               call fatal
            end if
!     Number of MD steps for each child trajectory
         else if (keyword(1:20) .eq. 'CHILD_STEPS ') then
            read(record,*,iostat=readstat) names,child_evol
            if (readstat .ne. 0) then
               write(*,*) "Correct format: CHILD_STEPS [Number of MD steps per child]"
               call fatal
            end if
!     If recrossing shall be calculated with MPI
         else if (keyword(1:20) .eq. 'MPI ') then
            recross_mpi = .true.
!     If recrossing part shall not be checked for errors
         else if (keyword(1:20) .eq. 'NO_CHECK ') then
            recross_check = .false.
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The RECROSS section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do
if (.not. section) then
   write(*,*) "Please add a RECROSS section with parameters to the keyfile!"
   call fatal
end if


write(15,*) "... done!"

! -------------------------------------------------------------------
!     F) SEVENTH, read in the FORCE section! (optional)
! -------------------------------------------------------------------

write(15,*) "F) The FORCE section..."
add_force=.false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   call upcase (record)
   string = record(next:120)
   if (trim(adjustl(record(1:11))) .eq. 'FORCE {' .or. trim(adjustl(record(1:11))) .eq. 'FORCE{') then
      add_force=.true.
      do j=1,nkey-i+1
         next=1
         record = keyline(i+j)
         call gettext (record,keyword,next)
         call upcase (keyword)
!     First atom on which force is applied: index and force vector
         if (keyword(1:11) .eq. 'VEC1 ') then
            read(record,*) names,force1_at,force1_k,force1_v(1),force1_v(2),force1_v(3)
!     Second atom on which force is applied: index and force vector
         else if (keyword(1:11) .eq. 'VEC2 ') then
            read(record,*) names,force2_at,force2_k,force2_v(1),force2_v(2),force2_v(3)
         end if
         if (keyword(1:13) .eq. '}') exit
         if (j .eq. nkey-i) then
            write(*,*) "The FORCE section has no second delimiter! (})"
            call fatal
         end if
      end do
   end if
end do


!
!      Check if all read in parameters are valid.
!      Abort the program if something is not valid
!
if (rank .eq. 0) then
!
!      General sampling trajectory settings 
!
   if (gen_step .eq. 0) then
      write(*,*) "No number of start structure generation MD steps defined!"
      write(*,*) "Add the keyword GENERATE_STEPS!"
      call fatal
   else if (equi_step .eq. 0) then
      write(*,*) "No number of umbrella equilibration MD steps defined!"
      write(*,*) "Add the keyword EQUILIBR_STEPS!"
      call fatal
   else if (umbr_step .eq. 0) then
      write(*,*) "No number of umbrella sampling MD steps defined!"
      write(*,*) "Add the keyword UMBRELLA_STEPS!"
      call fatal
   else if (umbr_traj .eq. 0) then
      write(*,*) "No number of single umbrella sampling trajectories per window defined!"
      write(*,*) "Add the keyword UMBRELLA_TRAJS!"
      call fatal
   end if
!
!      Thermostat settings 
!
   if (dt .eq. -1.0d0) then
      write(*,*) "No length of time step defined! Add the keyword DELTAT!"
      call fatal
   else if (kelvin .le. 0.0d0) then
      write(*,*) "No temperature defined! Add the keyword TEMP!"
      call fatal
   else if ((thermo .ne. "ANDERSEN") .and. (thermo .ne. "NOSE-HOOVER")) then
      write(*,*) "No availiable thermostat was chosen! The Andersen"
      write(*,*) "thermostat as well as the Nose-Hoover thermostat are "
      write(*,*) "available at the moment. Add the keyword THERMOSTAT!"
      call fatal
   end if
!
!      The RPMD beads 
!
   if (nbeads .lt. 1) then
      write(*,*) "No valid bead number given! Add the keyword BEAD_NUMBER!"
      call fatal
   end if
!
!     the number of equivalent reaction paths (multiplicator for k(T) value)
!
   if (npaths .lt. 1) then
      write(*,*) "No valid number of equivalent reaction paths given! At least NPATHS 1!"
      call fatal
   end if
!
!     specific umbrella sampling settings 
!
   if (k_force .eq. -1.0) then
      write(*,*) "No umbrella force constant defined! Add the keyword BIAS (UMBRELLA{)!"
      call fatal
   else if ((umbr_lo .eq. -5.0) .or. (umbr_hi .eq. -5.0)) then
      write(*,*) "No range for umbrella samplings defined! Add the keyword BONDS (UMBRELLA{)!"
      call fatal
   else if (umbr_dist .eq. -1.0) then
      write(*,*) "No distance between umbrella samplings defined! Add the keyword DIST (UMBRELLA{)!"
      call fatal
   else if (umbr_type .eq. "none") then
      write(*,*) "No type of biased coordinate system is given! Add the keyword TYPE (MECHA{)!"
      write(*,*) "Coordinate systems currently supported:"
!      write(*,*) "UNIMOLECULAR: A unimolecular reaction: one bond is built, one broken (one molecule)"
      write(*,*) "BIMOLEC: A bimolecular reaction: one bond is built, the other broken."
      write(*,*) "CYCLOADD: A cycloaddition: two bonds are built simultaneously."
      write(*,*) "ADDITION: A bimolecular addition: two bonds are built, one is broken."
      write(*,*) "ADDITION3: A trimolecular addition: three bonds are built, two are broken."
      write(*,*) "ADD3_SOLV: A trimolecular addition: three bonds are built, two are broken, two"
      write(*,*) "           reactants are clustered before the reaction (bimol. k(T))."
      write(*,*) "ADDITION4: A tetramolecular addition: four bonds are built, three are broken."
      write(*,*) "ADD4_SOLV: A tetramolecular addition: four bonds are built, three are broken, three"
      write(*,*) "           reactants are clustered before the reaction (bimol. k(T))."
      write(*,*) "ATOM_SHIFT: A single atom is transported a certain distance, e.g. on a metal surface."
      write(*,*) "CYCLOREV: A cycloreversion: two bonds are broken simultaneously."
      call fatal
   end if
!
!      prevent usage of unimolecular reactions 
!
   if (umbr_type .eq. "UNIMOLEC") then
      write(*,*) "Currently, this type is not implemented!"
      call fatal
   end if
!
!      The umbrella integration settings 
!
   if (xi_min .eq. 1.0 .or. xi_max .eq. -1.0) then
      write(*,*) "No range for PMF calculation defined! Add the keyword PMF_XI_RANGE!"
      call fatal
   else if (nbins .eq. 0) then
      write(*,*) "No number of windows for PMF calculation defined!"
      write(*,*) "Add the keyword PMF_BINS"
      call fatal
   end if
   call upcase(pmf_minloc)
   if (pmf_minloc .ne. 'ZERO' .and. pmf_minloc .ne. 'PMF_MIN') then
      write(*,*) "The keyword PMF_MINLOC was found but no valid option was"
      write(*,*) "chosen! Take ZERO or PMF_MIN!"
      call fatal
   end if
!
!      The asymptotic reactants distance (not needed for ATOM_SHIFT)
!
   if (R_inf .eq. 0.d0) then
     if ((umbr_type .ne. "ATOM_SHIFT") .and. (umbr_type .ne. "CYCLOREV")) then
        write(*,*) "No reactant max distance is given! Add the keyword DIST_INF!"
        call fatal
     end if
   end if
!
!      If no valid TS structure filename was given
!
   if (ts_file .eq. " ") then
      write(*,*) "No filename with the TS structure was given!"
      write(*,*) "Add the keyword TS_STRUC!"
      call fatal
   end if
!
!     check parameters for recrossing calculation
!
   if (.not. skip_recross) then
      if (recr_equi .eq. 0) then
         write(*,*) "No timestep number for the equilibration of the recrossing parent"
         write(*,*) "trajectory is given! Add the keyword RECROSS_EQUI!"
         call fatal
      else if (child_tot .eq. 0) then
         write(*,*) "No total number of child trajectories for the recrossing calculation"
         write(*,*) "is given! Add the keywod CHILD_TOTAL!"
         call fatal
      else if (child_interv .eq. 0) then
         write(*,*) "No interval between a new couple of child trajectories shall be spawned"
         write(*,*) "is given! Add the keyword CHILD_INTERVAL!"
         call fatal
      else if (child_point .eq. 0) then
         write(*,*) "No number of child trajectories that shall be spawned together at each"
         write(*,*) "spawining time is given! Add the keyword CHILD_PERPOINT!"
         call fatal
      else if (child_evol .eq. 0) then
         write(*,*) "No timestep number for each single child trajectory is given!"
         write(*,*) "Add the keyword CHILD_EVOL!"
      end if
      child_times=child_tot/child_point
      if (child_times .ne. real(child_tot/child_point)) then
         write(*,*) "The total number of child trajectories is no multiplicity of the child"
         write(*,*) "number per sampling point! Alter keywords CHILD_TOTAL and/or CHILD_PERPOINT!"
         call fatal
      end if
   else 
      child_times=100
!      recr_equi=1000
   end if
!
!      Check print polymer setting for printout of one trajectory: needs to be 
!       inside the valid sampling area 
!
   if (print_poly .gt. -5.d0) then
      if (print_poly .lt. xi_min .or. print_poly .gt. xi_max) then
         write(*,*) "The 'PRINT_POLYMER' option was activated such that a single umbrella"
         write(*,*) "shall be written out but the Xi location is not in the sampled area!"
         call fatal
      end if
   end if
end if
!
!      Read in dependent information (detailed subinformation of the one given
!      above) or calculate dependent parameters 
!
!     Define type of thermostat
!
if (thermo .eq. "ANDERSEN") then
   thermostat=1
else if (thermo .eq. "NOSE-HOOVER") then
   thermostat=2
else 
   write(*,*) "Only ANDERSEN or NOSE-HOOVER thermostats can be used!"
   call fatal
end if

!
!     Option for fixing of distinct atoms in order to make e.g. metal slab calculations 
!     possible 
!
if (fix_atoms) then
   inquire(file=fix_file,exist=exist)
   if (.not. exist) then
      if (rank .eq. 0) then
         write(*,*) "ERROR! The file ",fix_file," with the fixed atoms is not present!"
      end if
      call fatal
   end if
   allocate(fix_list(1000))
!
!     Read in list with fixed atoms 
! 
   open(unit=213,file=fix_file,status="old")
   k=1
   do
      read(213,*,iostat=readstat) fix_list(k)
      if (readstat .ne. 0) exit
      k=k+1
   end do
   fix_num=k-1
end if


!
!     Control if additional settings are correct/useful/appropriate
!
do i=1,sum_reacs
   if (at_reac(i,1) .eq. 0) then
      if (rank .eq. 0) then
         write(*,'(a,i1,a,i1,a)') "No atom for reactant ",i," is given! &
              &  Add the keyword REACTANT",i,"!"
         call fatal
      end if
   end if
end do

!
!     Control the elements of the reactants arrays for doublings, emptyness etc..
!     for all reactants in one loop!
!
do l=1,sum_reacs
   active=.true.
   i=0
   do while (active)
      i=i+1
      if (at_reac(l,i) .eq. 0) then
         active=.false.
      else if ((at_reac(l,i) .gt. natoms) .or. (at_reac(l,i) .lt. 0)) then
         if (rank .eq. 0) then
            write(*,'(a,i1,a)') "Reactant ",l," has an atom out of range!"
            call fatal
         end if
      else
         do p=1,l
            do j=1,n_reac(p)
               if (at_reac(l,i) .eq. at_reac(p,j)) then
                  if (rank .eq. 0) then
                    write(*,'(a,i4,a,i1,a,i1,a)') " The atom No. ",at_reac(l,i), &
                        & " is shared between reactant ",l," and ",p,"!"
                    call fatal
                  end if
               end if
            end do
         end do
      end if
   end do
   if (i .eq. 1) then
      if (rank .eq. 0) then
         write(*,*) "One of the reactants has no atoms!"
         call fatal
      end if
   end if
   n_reac(l)=i-1
end do
if (sum(n_reac) .ne. natoms .and. umbr_type .ne. "ATOM_SHIFT") then
   if (rank .eq. 0) then
      write(*,*) "The sum of reactant atoms isn´t equal to the total number of atoms!"
      call fatal
   end if
end if
!
!     allocate arrays for breaking and forming bonds as well as for its 
!     reading in
!

!
!      Continue here for the ATOM_SHIFT coordinate
!
22 continue 
!
!      Check methods for PMF integration. If none is set, set it to integration
!        as default 
!
call upcase(pmf_method)
if ((pmf_method .ne. 'WHAM') .and. (pmf_method .ne. 'INTEGRATION') .and. &
              & (pmf_method .ne. 'ALL')) then
   if (rank .eq. 0) then
      write(15,*) "No (existent) method for PMF calculation is given by the user!"
      write(15,*) "Umbrella integration will be used instead."
   end if
   pmf_method="INTEGRATION"
end if 

!
!      Read in the frequency of structure printouts for debug trajectory
!      If no keyword given, take 1 fs interval as default
!
if (print_poly .ne. -5.d0) then
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'TDUMP ') then
         read(record,*) names,dtdump
!
!     calculate the frequency of print outs
!
         iwrite = nint(dtdump/dt_info*0.001d0)
      end if
   end do
end if

!
!     For MECHANOCHEMISTRY simulations: Add forces to certain atoms!
!
!     Main keyword: ADD_FORCE
!
add_force = .false.
do i = 1, nkey
    next = 1
    record = keyline(i)
    call gettext (record,keyword,next)
    call upcase (keyword)
    string = record(next:120)
    if (keyword(1:11) .eq. 'ADD_FORCE ') then
       add_force=.true.
       if (rank .eq. 0) then
          write(*,*) "The keyword ADD_FORCE was found!"
          write(*,*) "Therefore a mechanochemistry trajectory will be"
          write(*,*) "calculated, until a rupture event is monitored and"
          write(*,*) "the user will be informed."
          write(*,*) "The keywords VEC1 and VEC2 must be present, which the formate:"
          write(*,*) "FORCE1/2 [atom index] [force constant (N)] [vector (x,y,z)]"
       end if
       exit
    end if
end do


!
!     If Mechanochemistry was activated: read in the atoms as well as 
!      vectors and strenghts for the applied forces 
!
if (add_force) then
!
!     normalize the force vectors
!
   force1_v=force1_v/(sqrt(force1_v(1)**2+force1_v(2)**2+force1_v(3)*2))
   force2_v=force2_v/(sqrt(force2_v(1)**2+force2_v(2)**2+force2_v(3)*2))
!
!     transform units of force constants from Newton to atomic units 
!
   if ((force1_at .eq. 0) .or. (force1_at .gt. natoms) .or.&
              &  (force2_at .eq. 0) .or. (force2_at .gt. natoms)) then
      write(*,*) "ERROR! The atoms at which the force shall be applied are"
      write(*,*) " not part of the system!"
      call fatal
   end if
   force1_k=force1_k*newton2au
   force2_k=force2_k*newton2au
!
!     Open file in which distances between the atoms with attached force will 
!     be written  
!
!   write(*,*) "Distances between atoms with attached force will be written to"
!   write(*,*) "file force_distances.dat"
!   open(unit=78,file="force_distances.dat",status="replace")
!   write(78,*) "# distances between the atoms with attached force, in Angstrom:"
end if


end subroutine calc_rate_read
