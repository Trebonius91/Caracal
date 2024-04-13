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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine read_pes: read in information about the potential energy 
!     surface to be used for the program currently in use.
!
!     part of EVB
! 
subroutine read_pes(rank)
use general
use evb_mod
use qmdff
use pbc_mod

implicit none
integer::num_arg,input_unit,i,qmdff_energies_unit,asciinum
character(len=70)::fffile1,fffile2,fffile3
real(kind=8)::energy,e_qmdff1,e_qmdff2,e_evb
integer::qmdffnumber,nat,nat3,l,m
integer::readstat  ! for error handling
integer::state_open  ! if a file was opened successfully
integer::nlines  ! number of structures in the IRC
integer::qmdff_index   ! number of QMDFF for corrections
integer::keylines_backup  ! for pGFN-FF bug
character(len=120)::string
character(len=20)::keyword
character(len=120)::record
character(len=70)::fileinfo,filegeo,ts_file
character(len=70)::init_struc_file
character(len=70)::file_irc_struc
character(len=70)::file_irc_ens
character(len=70)::filets,filets2,names
character(len=80)::coul_method,a80
character(len=1)::qmdffnum
real(kind=8),dimension(:,:),allocatable::coord
real(kind=8),dimension(:,:),allocatable::g_evb
real(kind=8),dimension(:,:),allocatable::xyz2,geo_xyz
real(kind=8),dimension(:),allocatable::int_coord,geo_int,geo_xyz1
real(kind=8),dimension(:,:),allocatable::ts_coordinates_a,ts_coordinates2_a
integer::mode,next,j,k,readstatus,dg_evb_mode,mat_size,num_struc
integer::int_mode ! method of defining internal coordinates
real(kind=8)::s_q1_r  ! temporary variable for QMDFF damping range (RP-EVB)
real(kind=8),allocatable::xyz_init(:,:)  ! geometry for pGFN-FF init
logical::read_init_struc
character(len=2),allocatable::names_init(:)   ! element symbols for pGFN-FF init
logical::path_struc,path_energy,coupl,params
logical::evb1,evb2,evb3,ffname1,ffname2,ffname3,defqmdff
logical::exist,exists,has_next,coupl1
logical::read_name
logical::ts_xyz
integer::rank ! the current MPI rank
character(len=80)::corr_name,method
character(len=100)::orca_filename
real(kind=8)::epsilon_val

evb1=.false.
evb2=.false.
evb3=.false.
ffname1=.false.
ffname2=.false.
ffname3=.false.
path_struc=.false.
coupl=.false.
coupl1=.false.
defqmdff=.false.
params=.false.
filegeo="dummy" ! is not used in this context
orca=.false.
water_spc=.false.
pes_topol=.false.
!
!     mainly for test reasons: activate numerical gradient
!
num_grad=.false.
full=.false.
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'NUM_GRAD ') then
      num_grad=.true.
!      if (rank .eq. 0) then
!         write(*,*) "The numerical gradient will be calculated!"
!      end if
      exit
   end if
end do
!
!     if numgrad is activated, ask for the step of elongations
!
if (num_grad) then
   num_grad_step=0.0001d0 ! default value
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:20) .eq. 'NUM_GRAD_STEP ') then
         num_grad=.true.
         read(record,*) names,num_grad_step
         if (rank .eq. 0) then
            write(*,*) "The numerical gradient step will be:", num_grad_step
         end if
      end if
   end do
end if
!
!     If the system shall be calculated with periodic boundary 
!     conditions 
!
periodic=.false.
box_walls=.false.
boxlen_x=0.d0
boxlen_y=0.d0
boxlen_z=0.d0

do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:13) .eq. 'PERIODIC ') then
      read(record,*,iostat=readstat) names,boxlen_x,boxlen_y,boxlen_z
      periodic=.true.
      if (readstat .ne. 0) then
         write(*,*) "Correct format: PERIODIC [xlen,ylen,zlen (A)]"
         call fatal
      end if
   else if (keyword(1:11) .eq. 'BOX_WALLS ') then
      read(record,*,iostat=readstat) names,boxlen_x,boxlen_y,boxlen_z
      box_walls=.true.
      if (readstat .ne. 0) then
         write(*,*) "Correct format: BOX_WALLS [xlen,ylen,zlen (A)]"
         call fatal
      end if
   end if
end do

!
!     The general Method keyword
!
method=""
do i = 1, nkey_lines
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'PES ') then
      read(record,*) names,method
      call upcase(method)
      exit
   end if
end do

!
!      Choose the EVB coupling term to be optimized
!
evb_dq=.false.
dg_evb=.false.
treq=.false.
if (method .eq. "QMDFF")  then
   if (rank .eq. 0) then
      qmdffnumber=1
   end if
else if (method .eq. "DE_EVB")  then
   evb_de=.true.
   if (rank .eq. 0) then
      write(15,*) "The dE-EVB (energy) coupling term will be used!"
      write(15,*)
      write(*,*) "The dE-EVB (energy) coupling term will be used!"
      write(*,*)
   end if
else if (method .eq. "DQ_EVB")  then
   evb_dq=.true.
   if (rank .eq. 0) then
      write(15,*) "The dQ-EVB (coordinate) coupling term will be used!"
      write(15,*)
      write(*,*) "The dQ-EVB (coordinate) coupling term will be used!"
      write(*,*)
   end if
else if (method .eq. "DG_EVB") then
   dg_evb=.true.
   if (rank .eq. 0) then
      write(15,*) "The DG-EVB (distributed Gaussian) coupling term will be used!"
      write(15,*)
      write(*,*) "The DG-EVB (distributed Gaussian) coupling term will be used!"
      write(*,*)
   end if
else if (method .eq. "TREQ") then
   treq=.true.
   if (rank .eq. 0) then
      write(15,*) "The TREQ potential energy function will be used!"
      write(15,*)
      write(*,*) "The TREQ potential energy function will be used!"
      write(*,*)
   end if
!
!     If the SPC water model shall be used (only for water, of course...)
!
else if (method .eq. "WATER_SPC") then
   water_spc=.true.
   if (rank .eq. 0) then
   end if

else if (method .eq. "GFN-XTB") then
   gfn_xtb = .true.
   if (rank .eq. 0) then
      write(*,*) "The included tblite program will be called for "
      write(*,*) " calculation of the GFN-xTB energy and gradient."
      write(*,*) "Details of the calculations are written to 'gfn_xtb.log'."
   end if
!
!     Counter for better listing of single SCFs in gfn_xtb.log
!   
   xtb_calc_num = 1
   open(unit=84,file="gfn_xtb.log",status="replace")
!
!     If the gradient shall be calculated on the fly with orca
!
else if (method .eq. "PGFN-FF") then
   pgfn_ff = .true.
   if (rank .eq. 0) then
   !   write(*,*) "The GULP program library will be called for"
   !   write(*,*) " calculation of the pGFN-FF energy and gradient."
   end if

else if (method .eq. "ORCA") then
   orca=.true.
   if (rank .eq. 0) then
   end if
!
!     If the energy/gradient shall be calculated on the fly with an arbitrary 
!      external program
!
else if (method .eq. "EXTERNAL") then
   call_ext=.true.
!
!     If the program calc_rate is used, the number of atoms must be read in from the 
!     start structure (TS) since the number of atoms is not clear from the beginning!
!
   if (use_calc_rate) then
      do i = 1, nkey_lines
         next = 1
         record = keyline(i) 
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
!  
!     reset the read status to successful (check it for each keyword)
!
         readstat= 0
         if (keyword(1:20) .eq. 'TS_STRUC ') then
            read(record,*,iostat=readstat) names,ts_file
            if (readstat .ne. 0) then
               write(*,*) "Correct format: TS_STUC [filename]"
               call fatal
            end if
            exit
         end if
      end do
!
!     If the TS file has been found, read in the number of atoms from its first line
!
      open(unit=37,file=ts_file,status="old")
      read(37,*) natoms
      close(37)
   end if

   if (rank .eq. 0) then
      write(15,*) "An external potential energy function will be used!"
      write(15,*)
      write(*,*) "An external potential energy function will be used!"
      write(*,*)
   end if

!
!     If a new PES has been included by a user, this is the preferred way to connect
!     it to the code
!
else if (method .eq. "CUSTOM") then
!
!     A general flag is activated, indicating the usage of the imported method
!   
   call_cust=.true.

!
!     If the program calc_rate is used, the number of atoms must be read in from the 
!     start structure (TS) since the number of atoms is not clear from the beginning!
!
   if (use_calc_rate) then
      do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
!  
!     reset the read status to successful (check it for each keyword)
!
         readstat= 0
         if (keyword(1:20) .eq. 'TS_STRUC ') then
            read(record,*,iostat=readstat) names,ts_file
            if (readstat .ne. 0) then
               write(*,*) "Correct format: TS_STUC [filename]"
               call fatal
            end if
            exit
         end if
      end do
!
!     If the TS file has been found, read in the number of atoms from its first line
!
      open(unit=37,file=ts_file,status="old")
      read(37,*) natoms
      close(37)
   end if

!
!     The analytical potential energy surfaces! Mainly for benchmark reasons(?)
!
else if (method .eq. "ANA_H3") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The H+H2 potential energy surface is used."
      write(15,*) "References:  D. G. Truhlar, C. J. Horowitz J. Chem. Phys.,"
      write(15,*) "             Vol. 68, p. 2466, 1978."
      write(15,*)
      write(*,*) "The H+H2 potential energy surface is used."
      write(*,*) "References:  D. G. Truhlar, C. J. Horowitz J. Chem. Phys.,"
      write(*,*) "             Vol. 68, p. 2466, 1978."
      write(*,*)
!     (no initialization needed)
   end if
   natoms=3  ! store number of atoms..
   pot_ana = .true.
   pot_type = "h3"

else if (method .eq. "ANA_BRH2") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The HBr+H potential energy surface is used"
      write(15,*) "References:  D. C. Clary, Chem. Phys. 71, 117 (1982)"
      write(15,*) "             I. Last and M. Baer, in Potential Energy Surfaces "
      write(15,*) "             and Dynamics, edited by D. G. Truhlar, p. 519 "
      write(15,*)
      write(*,*) "The HBr+H potential energy surface is used"
      write(*,*) "References:  D. C. Clary, Chem. Phys. 71, 117 (1982)"
      write(*,*) "             I. Last and M. Baer, in Potential Energy Surfaces "
      write(*,*) "             and Dynamics, edited by D. G. Truhlar, p. 519 "
      write(*,*)
   end if
   call initialize_brh2
   pot_type="brh2"
   natoms=3  ! store number of atoms..

else if (method .eq. "ANA_O3") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The O2+O potential energy surface is used"
      write(15,*) "References:  Z. Varga, Y. Paukku, and D. G. Truhlar"
      write(15,*) "             J. Chem. Phys. 147, 154312/1-17 (2017) "
      write(15,*)
      write(15,*) "The O2+O potential energy surface is used"
      write(15,*) "References:  Z. Varga, Y. Paukku, and D. G. Truhlar"
      write(15,*) "             J. Chem. Phys. 147, 154312/1-17 (2017) "
      write(15,*)
   end if
   pot_type="o3"
   natoms=3
else if (method .eq. "ANA_OH3") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The OH2+H potential energy surface is used"
      write(15,*) "References:  G. C. Schatz and H. Elgersma"
      write(15,*) "             Chem. Phys. Lett. 73, 21 (1980) "
      write(15,*)
      write(*,*) "The OH2+H potential energy surface is used"
      write(*,*) "References:  G. C. Schatz and H. Elgersma"
      write(*,*) "             Chem. Phys. Lett. 73, 21 (1980) "
      write(*,*)
   end if
   pot_type="oh3"
   call initialize_oh3
   natoms=4  ! store number of atoms..
   pot_ana=.true.
else if (method .eq. "ANA_H2CO") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The H2CO unimolecular potential energy surface is used"
      write(15,*) "References:  X. Wang, P. L. Houston and J. M. Bowman"
      write(15,*) "             Phil. Trans. R. Soc. A 375, 20160194 (2016)"
      write(15,*)
      write(*,*) "The H2CO unimolecular potential energy surface is used"
      write(*,*) "References:  X. Wang, P. L. Houston and J. M. Bowman"
      write(*,*) "             Phil. Trans. R. Soc. A 375, 20160194 (2016) "
      write(*,*)
   end if
   pot_type="h2co"
   call initialize_h2co
   natoms=4  ! store number of atoms..
   pot_ana=.true.
else if (method .eq. "ANA_CLNH3") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The NH2+HCl potential energy surface is used"
      write(15,*) "References:  M. Monge-Palacios, C. Rangel, J.C. Corchado and J. Espinosa-Garcia"
      write(15,*) "             Int.J.Quantum.Chem. 112, 1887 (2012) "
      write(15,*)
      write(*,*) "The NH2+HCl potential energy surface is used"
      write(*,*) "References:   M. Monge-Palacios, C. Rangel, J.C. Corchado and J. Espinosa-Garcia"
      write(*,*) "              Int.J.Quantum.Chem. 112, 1887 (2012) "
      write(*,*)
   end if
   pot_type="clnh3"
   call initialize_clnh3
   natoms=5  ! store number of atoms..
   pot_ana=.true.
else if (method .eq. "ANA_CH4H") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The CH4H potential energy surface is used."
      write(15,*) "Reference: J.C. Corchado, J.L. Bravo and J. Espinosa-Garcia "
      write(15,*) "             J.Chem.Phys.130,184314 (2009)."
      write(15,*)
      write(*,*) "The CH4H potential energy surface is used."
      write(*,*) "Reference: J.C. Corchado, J.L. Bravo and J. Espinosa-Garcia "
      write(*,*) "             J.Chem.Phys.130,184314 (2009)."
      write(*,*)
   end if
   pot_type="ch4h"
   call initialize_ch4h
   natoms=6
   pot_ana=.true.
else if (method .eq. "ANA_NH3OH") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The NH3OH potential energy surface is used."
      write(15,*) "References: M. Monge-Palacios, C. Rangel and J. Espinosa-Garcia "
      write(15,*) "           J.Chem.Phys. 138, 084305 (2013)."
      write(15,*)
      write(*,*) "The NH3OH potential energy surface is used."
      write(*,*) "References: M. Monge-Palacios, C. Rangel and J. Espinosa-Garcia "
      write(*,*) "           J.Chem.Phys. 138, 084305 (2013)."
      write(*,*)
   end if
   pot_type="nh3oh"
   call initialize_nh3oh
   natoms=6
   pot_ana=.true.
else if (method .eq. "ANA_CH4OH") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The CH4OH potential energy surface is used."
      write(15,*) "References:  J. Espinosa-Garcia, J. C. Corchado, J. Chem. Phys.,"
      write(15,*) "             Vol. 112, p. 5731, 2000."
      write(15,*)
      write(*,*) "The CH4OH potential energy surface is used."
      write(*,*) "References:  J. Espinosa-Garcia, J. C. Corchado, J. Chem. Phys.,"
      write(*,*) "             Vol. 112, p. 5731, 2000."
      write(*,*)
   end if
   pot_type="ch4oh"
   call initialize_ch4oh
   natoms=7
   pot_ana=.true.
else if (method .eq. "ANA_CH4CN") then
   if (rank .eq. 0) then
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The CH4CN potential energy surface is used."
      write(15,*) "References:  J. Espinosa-Garcia, C. Rangel and Y. V. Suleimanov,"
      write(15,*) "             Phys. Chem. Chem. Phys., 2017, 19, 19341."
      write(15,*)
      write(*,*) "The CH4CN potential energy surface is used."
      write(*,*) "References:  J. Espinosa-Garcia, C. Rangel and Y. V. Suleimanov,"
      write(*,*) "             Phys. Chem. Chem. Phys., 2017, 19, 19341."
      write(*,*)
   end if
   pot_type="ch4cn"
   call initialize_ch4cn
   natoms=7
   pot_ana=.true.
else if (method .eq. "ANA_GEH4OH") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The GeH4+OH potential energy surface is used."
      write(15,*) "References:  J. Espinosa-Garcia, C. Rangel and J.C. Corchado"
      write(15,*) "             Phys. Chem. Chem. Phys. 18, 16941 (2016)"
      write(15,*)
      write(*,*) "The GeH4+OH potential energy surface is used."
      write(*,*) "References:  J. Espinosa-Garcia, C. Rangel and J.C. Corchado"
      write(*,*) "             Phys. Chem. Chem. Phys. 18, 16941 (2016)"
      write(*,*)
   end if
   pot_type="geh4oh"
   call initialize_geh4oh
   natoms=7  ! store number of atoms..
   pot_ana=.true. 
else if (method .eq. "ANA_C2H7") then
   if (rank .eq. 0) then
      open(unit=16,file="pot_info.out",status="unknown")
      write(16,*) "Informational output about the used analytical PES:"
      if (rank .eq. 0) then
         write(*,*) "----------"
      end if
      write(15,*) "The C2H6+H potential energy surface is used."
      write(15,*) "References:  Arindam Chakraborty, Yan Zhao, Hai Lin, and Donald G. Truhlar"
      write(15,*) "             J. Chem. Phys., 124, 044315 (2006)."
      write(15,*)
      write(*,*) "The C2H6+H potential energy surface is used."
      write(*,*) "References:  Arindam Chakraborty, Yan Zhao, Hai Lin, and Donald G. Truhlar"
      write(*,*) "             J. Chem. Phys., 124, 044315 (2006)."
      write(*,*)
   end if
   pot_type="c2h7"
!   (no initialization needed)
   natoms=9  ! store number of atoms..
   pot_ana=.true.
else if (method .eq. "TOPOL") then
   if (rank .eq. 0) then
      if (.not. use_explore) then
         write(*,*) "Currently, the PES TOPOL keyword can only be used with the explore program!"
         call fatal
      end if
      pes_topol = .true.
      write(*,*) "A topology analysis will be done for all structures given."
      write(*,*) " A list of all bonds, angles, dihedrals and out of planes will be"
      write(*,*) " given as well as the Wilson matrix elements and derivatives."
      
   end if
else
   if (rank .eq. 0) then
      write(*,*) "No valid potential energy surface was chosen!"
      write(*,*) "Choose a valid one:"
      write(*,*) "QMDFF, DE_EVB, DQ_EVB, DG_EVB, TREQ, ORCA (QM call),"
      write(*,*) "EXTERNAL (providing your own energy/gradient program),"
      write(*,*) "analytical PES (ANA_H3, ANA_BRH2, ANA_O3, ANA_OH3, "
      write(*,*) " ANA_CH4H, ANA_NH3OH, ANA_CH4OH, ANA_GEH4OH, ANA_C2H7), "
      write(*,*) "TOPOL (coordinate analysis with explore program)"
      call fatal
   end if
end if

if (pot_ana) then
   if (rank .eq. 0) then
      write(*,*) "You have chosen one of the analytical PES functions."
      write(*,*) "The following surfaces are available as well (ascending atom number)"
      write(*,*) "1)   H2 + H    (H3)      (3 atoms: H, H, H)"
      write(*,*) "2)   BrH + H   (BrH2)    (3 atoms: H, Br, H)"
      write(*,*) "3)   O2 + O    (OH3)     (3 atoms: O, O, O)"
      write(*,*) "4)   OH2 + H   (OH3)     (4 atoms: O, H, H, H)"
      write(*,*) "5)   H2CO      (H2CO)    (4 atoms: C, O, H, H)"
      write(*,*) "6)   NH2 + HCl (ClNH2)   (5 atoms: H, N, H, H, Cl)"
      write(*,*) "7)   CH4 + H   (CH4H)    (6 atoms: H, C, H, H, H, H)"
      write(*,*) "8)   NH3 + OH  (NH3OH)   (6 atoms: H, N, H, H, O, H)"
      write(*,*) "9)   CH4 + OH  (CH4OH)   (7 atoms: H, C, H, H, H, O, H)"
      write(*,*) "10)   CH4 + CN  (CH4CN)   (7 atoms: H, C, H, H, H, C, N)"
      write(*,*) "11)   GeH4 + OH (GeH4OH)  (7 atoms: H, Ge, H, H, H, O, H)"
      write(*,*) "12)   C2H6 + H  (C2H7)    (9 atoms: C, H, H, H, C, H, H, H, H)"
      write(*,*)
      write(15,*) "Infos about potential initialization written to file pot_info.out."
      write(15,*)
      write(*,*) "Infos about potential initialization written to file pot_info.out."
      write(*,*)
   end if
   close(16)
   return
end if

!
!     Call the water model initialization routine if needed 
!     For the egrad program, read in the number of atoms from the structure file 
!      (first line in a xyz file..)
!
if (water_spc) then
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:16) .eq. 'COORDS_FILE ') then
         read(record,*,iostat=readstat) names,filegeo
         open(unit=38,file=filegeo,status="old") 
         read(38,*) natoms
         read(38,*)
         do j=1,natoms
            read(38,*) name(j)
         end do
         close(38) 
      end if
   end do


   call water_init
   goto 678
end if

!
!     For GFN-xTB calculations: initialize the respective object
!
if (gfn_xtb) then
!
!      Set default values 
!
   hamil_string = "xxx"
   xtb_charge = 0
   xtb_accuracy = 1.0d-6
   xtb_maxiter = 200
   xtb_el_temp = 0.0d0
   solv_string = "none"
   solv_spec = "none"
   epsilon_val = -1000d0
   exist_spec = .false.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'GFN-XTB {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'GFN-XTB{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
!
!     The GFN-xTB Hamiltonian, available are: GFN1-xTB, GFN2-xTB and IPEA1-xTB
!
            if (keyword(1:18) .eq. 'HAMILTONIAN ') then
               read(record,*,iostat=readstat) names,hamil_string
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The HAMILTONIAN keyword seems to be corrupted!"
                  end if
                  call fatal
               end if                
               call upcase (hamil_string) 
               if (trim(adjustl(hamil_string)) .eq. "GFN1-XTB") then
                  if (rank .eq. 0) then
         !            write(*,*) "The GFN1-xTB Hamiltonian will be used."
                  end if
               else if (trim(adjustl(hamil_string)) .eq. "GFN2-XTB") then   
                  if (rank .eq. 0) then
         !            write(*,*) "The GFN2-xTB Hamiltonian will be used."
                  end if
               else if (trim(adjustl(hamil_string)) .eq. "IPEA1-XTB") then
                  if (rank .eq. 0) then
         !            write(*,*) "The IPEA1-xTB Hamiltonian will be used."
                  end if
               else 
                  if (rank .eq. 0) then
                     write(*,*) "No valid GFN-xTB Hamiltonian has been chosen!"
                     write(*,*) "Either GFN1-xTB, GFN2-xTB or IPEA1-xTB can be used."
                     call fatal
                  end if
               end if
!
!     The total charge of the system
!
            else if (keyword(1:18) .eq. 'CHARGE ') then
               read(record,*,iostat=readstat) names,xtb_charge
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The CHARGE keyword seems to be corrupted!"
                  end if
                  call fatal
               end if

!
!     The SCF convergence criterion: maximum energy change
!
            else if (keyword(1:18) .eq. 'DELTAE_CONV ') then
               read(record,*,iostat=readstat) names,xtb_accuracy
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The DELTAE_CONV keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
!
!     The maximum number of steps for the SCF cycle
!
            else if (keyword(1:18) .eq. 'SCF_MAXITER ') then
               read(record,*,iostat=readstat) names,xtb_maxiter
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The SCF_MAXITER keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
!
!     The electronic temperature for the Fermi occupation smearing
!
            else if (keyword(1:18) .eq. 'EL_TEMP ') then
               read(record,*,iostat=readstat) names,xtb_el_temp
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The EL_TEMP keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
!
!     The model for implicit solvation (if activated at all)
!
            else if (keyword(1:18) .eq. 'SOLV_MODEL ') then
               read(record,*,iostat=readstat) names,solv_string
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The SOLV_MODEL keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
               call upcase (solv_string)
               if (trim(adjustl(solv_string)) .eq. "ALPB") then
                  if (rank .eq. 0) then
         !            write(*,*) "The ALPB solvation model will be used."
                  end if
               else if (trim(adjustl(solv_string)) .eq. "CPCM") then
                  if (rank .eq. 0) then
          !           write(*,*) "The CPCM solvation model will be used."
                  end if
               else 
                  if (rank .eq. 0) then
                     write(*,*) "No valid solvation model chosen! Take either ALPB or CPCM."
                  end if
                  call fatal
               end if
!
!     One of the predefined solvent species
!
            else if (keyword(1:16) .eq. 'SOLV_SPEC ') then
               read(record,*,iostat=readstat) names,solv_spec
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The SOLV_SPEC keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
               call lowcase(solv_spec)
               exist_spec = .true.
!
!     An individual epsilon keyword for the solvation
!
            else if (keyword(1:17) .eq. 'SOLV_EPSILON ') then
               read(record,*,iostat=readstat) names,epsilon_val
               write(solv_epsilon,'(f20.12)') epsilon_val
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The SOLV_EPSILON keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The GFN-XTB section has no second delimiter! (})"
               call fatal
            end if

         end do
      end if
   end do
!
!     Check if any Hamiltonian has been chosen
!
   if (trim(adjustl(hamil_string)) .eq. "xxx") then
      write(*,*) "No valid GFN-xTB Hamiltonian has been chosen!"
      write(*,*) "Please add the HAMILTONIAN keyword!"
      call fatal
   end if
!
!     Check the validity of implicit solvation input
!
   if (trim(adjustl(solv_string)) .ne. "none") then
      if (exist_spec) then
         if (trim(adjustl(solv_spec)) .eq. "none") then
            write(*,*) "Implicit solvation was activated but no valid solvent"
            write(*,*) " species or epsilon is given!"
            call fatal         
         end if 
      else 
         if (epsilon_val + 1000.d0 .le. 0.0001d0) then
            write(*,*) "Implicit solvation was activated but no valid solvent"
            write(*,*) " species or epsilon is given!"
            call fatal
         end if
      end if
   end if
!
!     Convert the accuracy to internal unit
!
   xtb_accuracy = xtb_accuracy*1.0D6
!
!     Convert the electronic temperature from Kelvin to kT in eV
!
   xtb_el_temp = xtb_el_temp*8.61732476104d-5
   goto 678

end if

!
!     If the topology of a structure or of a number of structures in a 
!     trajectory shall be analyzed with explore.x, read in the method for 
!     bond determination (van der Waals radii or extended Hückel theory)
!
if (pes_topol) then
   topol_vdw_scale = 1.0
   topol_eht_cutoff = 1.0 
   topol_bonds = "EHT"    !  default: extended Hückel calculation
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'TOPOL {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'TOPOL{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
!
!     In which way the bonds of the system shall be calculated: VDW (comparison to
!       pairwise van-der-Waals radii) or EHT (Extended Hückel theory)
!
            if (keyword(1:14) .eq. 'BONDS ') then
               read(record,*,iostat=readstat) names,topol_bonds
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The BONDS keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
               call upcase (topol_bonds)
               if (trim(adjustl(topol_bonds)) .eq. "VDW") then
                  if (rank .eq. 0) then
                  end if
               else if (trim(adjustl(topol_bonds)) .eq. "EHT") then
                  if (rank .eq. 0) then
                  end if
               else
                  if (rank .eq. 0) then
                     write(*,*) "No valid TOPOL BONDS has been chosen!"
                     write(*,*) "Either VDW or EHT can be used."
                     call fatal
                  end if
               end if
!
!     Scaling for VDW bond comparisons
!
            else if (keyword(1:14) .eq. 'VDW_SCALE ') then
               read(record,*,iostat=readstat) names,topol_vdw_scale
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The VDW_SCALE keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
!
!     Cutoff for Wiberg-Mayer bond orders in EHT
!

            else if (keyword(1:14) .eq. 'EHT_CUTOFF ') then
               read(record,*,iostat=readstat) names,topol_eht_cutoff
               if (readstat .ne. 0) then
                  if (rank .eq. 0) then
                     write(*,*) "The EHT_CUTOFF keyword seems to be corrupted!"
                  end if
                  call fatal
               end if
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The TOPOL section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do

   goto 678

end if

!
!     For pGFN-FF calculations with the GULP library: Read in the initial
!      structure but no other keywords so far
!      
!
if (pgfn_ff) then
!
!     If the GULP libraries were not included in the compilation, give an error
!
#ifdef GULP

#else
      write(*,*)
      write(*,*) "Please compile Caracal with GULP to enable pGFN-FF calculations!"
      write(*,*) " Look into the Caracal wiki for details."
      call fatal
#endif
   read_init_struc = .false.
!
!      Set default values 
!
   keylines_backup = nkey_lines

   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'PGFN-FF {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'PGFN-FF{') then
         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
!
!     The initial structure for bond matrix initialization 
!
            if (keyword(1:16) .eq. 'STRUC_INIT ') then
               if (.not. read_init_struc) then
                  read(record,*,iostat=readstat) names,filegeo
                  open(unit=38,file=filegeo,status="old")
                  read(38,*) natoms
                  read(38,*)
                  allocate(xyz_init(3,natoms))
                  allocate(names_init(natoms))
                  do k=1,natoms
                     read(38,*) names_init(k),xyz_init(:,k)
                  end do
                  close(38)
                  init_struc_file=filegeo
                  read_init_struc = .true.
               end if   
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The EXTERNAL section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do
   if (.not. read_init_struc) then
      write(*,*) "Please give the keyword STRUC_INIT in the PGFN-FF section!"
      call fatal
   end if
!
!     Call the actual initialization routine
!
   call pgfn_init(natoms,xyz_init,names_init)
   
   nkey_lines=keylines_backup
   goto 678
end if
!
!     For orca calculations: Read in the orca header line of the orca 
!       input file 
!      It must contain everything above the charge and multiplicity line:
!    ! PALx MP2 ....
!    %scf ... [special settings]
!      end
!
if (orca) then

   if (rank .eq. 0) then
      write(*,*) 
      write(*,*) "Used PES: direct call to orca QM package"
      write(*,*) "Each structure will be submitted to an external orca"
      write(*,*) " calculation for each bead, separately!"
   end if
!
!      Set default values 
!
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'ORCA {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'ORCA{') then

         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
            if (keyword(1:18) .eq. 'HEADER_FILE ') then
               read(record,*) names,orca_filename
            else if (keyword(1:16) .eq. 'CHARGE ') then
               read(record,*) names,orca_charge
            else if (keyword(1:16) .eq. 'MULTI ') then
               read(record,*) names,orca_multi
            else if (keyword(1:16) .eq. 'SYMLINK ') then
               read(record(8:120),'(a)') call_orca
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The ORCA section has no second delimiter! (})"
               call fatal
            end if

         end do
      end if
   end do
   open(unit=78,file=orca_filename,status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",trim(orca_filename)," with the orca commands is not"
      write(*,*) " there or corrupted!"
      call fatal
   end if
   orca_com_num=1
   do 
      read(78,'(a)',iostat=readstat) orca_com(orca_com_num)
      if (readstat .ne. 0) exit
      orca_com_num=orca_com_num+1
   end do
   close(78)

   goto 678
end if
!
!    For calculations with external arbitrary programs: read in the 
!     link/symlink of that program
!

if (call_ext) then

   if (rank .eq. 0) then
      write(*,*)
      write(*,*) "Used PES: direct call to an external program"
      write(*,*) "Each structure will be submitted to an external "
      write(*,*) " calculation for each bead, separately!"
   end if
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:13))) .eq. 'EXTERNAL {' .or. trim(adjustl(record(1:13))) &
              &  .eq. 'EXTERNAL{') then

         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
            if (keyword(1:16) .eq. 'SYMLINK ') then
               read(record(8:120),'(a)') symlink_ext
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The EXTERNAL section has no second delimiter! (})"
               call fatal
            end if

         end do
      end if
   end do
   goto 678
end if

!
!     For calculation with a custom PES routine: If you want to use 
!     several PESs with one Caracal program, simply choose the right one 
!     with the integer
!

if (call_cust) then

   if (rank .eq. 0) then
      write(*,*)
      write(*,*) "Used PES: call to (one of) your custom routines"
      write(*,*) "The right routine will be chosen based on the specifier"
      write(*,*) " given after the keyword PES_NUMBER"
   end if
   cust_number = 0
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:13))) .eq. 'CUSTOM {' .or. trim(adjustl(record(1:13))) &
              &  .eq. 'CUSTOM{') then

         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            record=adjustl(record)
            if (keyword(1:17) .eq. 'PES_NUMBER ') then
               read(record,*) names, cust_number
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The CUSTOM section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do
   if (cust_number .lt. 1) then
      write(*,*) "You have chosen to use a custom PES routine, but did not specify"
      write(*,*) " its index with the keyword PES_NUMBER!"
      call fatal
   else 
!
!     The initialization management routine is called, where global
!     PES parameters are defined or read in if needed
!
      call custom_init

   end if
   goto 678
end if

!
!     If QMDFF or EVB-QMDFF shall be used as method, read in additional parameters
!

if (qmdffnumber .eq. 1 .or. evb_de .or. evb_dq .or. dg_evb .or. treq) then
   evb2=.false.
   evb3=.false.
   read_name=.false.
   exist=.false.
   ewald=.false.
   zahn=.true.
   coul_method="ZAHN"
   coul_cut=10.d0
   vdw_cut=10.d0
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      call upcase (record)
      string = record(next:120)
      if (trim(adjustl(record(1:11))) .eq. 'QMDFF {' .or. trim(adjustl(record(1:11))) &
           &  .eq. 'QMDFF{') then

         do j=1,nkey_lines-i+1
            next=1
            record = keyline(i+j)
            call gettext (record,keyword,next)
            call upcase (keyword)
            string = record(next:120)
!
!     Read in the names of the QMDFF files
!
            if (keyword(1:11) .eq. 'FFNAMES ') then
               if (qmdffnumber .eq. 1) then
                  read(record,*,iostat=readstat) names,fffile1
                  if (readstat .ne. 0) then
                     write(*,*) "Please check the QMDFFNAMES keyword!"
                     call fatal
                  end if
               end if
               read(record,*,iostat=readstat) names,fffile1,fffile2,fffile3
               ffname3=.true.
               defqmdff=.true.
               if (readstat .ne. 0) then
                  read(record,*,iostat=readstat) names,fffile1,fffile2
                  ffname2=.true.
                  defqmdff=.true.
                  if (readstat .ne. 0) then
                     qmdffnumber=1
                  else if (readstat .eq. 0) then
                     evb2=.true.
                     qmdffnumber=2
                  end if
               else if (readstat .eq. 0) then
                  evb3=.true.
                  qmdffnumber=3
               end if
               read_name = .true.
!
!     Read in the relative energy shift of the QMDFFs 
!
            else if (keyword(1:11) .eq. 'ESHIFT ') then
               if (qmdffnumber .eq. 1) then
                  read(record,*,iostat=readstat) names,E_zero1
                  exist=.true.
                  if (readstat .ne. 0) then
                     write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                     call fatal
                  end if
               end if
               if (evb2) then
                  read(record,*,iostat=readstat) names,E_zero1,E_zero2
                  exist=.true.
                  if (readstat .ne. 0) then
                     write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                     call fatal
                  end if
               end if
               if (evb3) then
                  read(record,*,iostat=readstat) names,E_zero1,E_zero2,E_zero3
                  exist=.true.
                  if (readstat .ne. 0) then
                     write(*,*) "The ESHIFT keyword in the QMDFF section seems to be corrupted!"
                     call fatal
                  end if
               end if
!
!     The method used for coulomb interactions 
!
            else if (keyword(1:11) .eq. 'COULOMB ') then
               read(record,*,iostat=readstat) names,a80
               if (readstat .ne. 0) then
                  write(*,*) "Correct format: COULOMB [method]"
                  call fatal
               end if
               call upcase(a80)
               if (a80 .eq. "PPPME") then
                  ewald=.true.
                  write(*,*) "The Particle Particle Particle Mesh Ewald (PPPME) method"
                  write(*,*) " will be used for the handling of Coulomb interactions."
                  write(*,*) "Please choose another method, the implementation is not complete!"
                 call fatal
               end if
               if (a80 .eq. 'ZAHN ') then
                  zahn=.true.
               end if
               coul_method=a80
!    The cutoff distance for Coulomb interactions
            else if (keyword(1:14) .eq. 'COUL_CUTOFF ') then
               read(record,*,iostat=readstat) names,coul_cut
               if (readstat .ne. 0) then
                  write(*,*) "Correct format: COUL_CUTOFF [value (A)]"
                  call fatal
               end if
!    The cutoff distance for VDW interactions
            else if (keyword(1:13) .eq. 'VDW_CUTOFF ') then
               read(record,*,iostat=readstat) names,vdw_cut
               if (readstat .ne. 0) then
                  write(*,*) "Correct format: VDW_CUTOFF [value (A)]"
                  call fatal
               end if
!           
!     Determine if the QMDFF energy shifts shall be corrected automatically
!     to exactly reproduce the first/last energy of the path or not
!    
            else if (keyword(1:16) .eq. 'SHIFT_MANUAL') then
               shift_man=.true.
            end if
            if (keyword(1:11) .eq. '}') exit
            if (j .eq. nkey_lines-i) then
               write(*,*) "The QMDFF section has no second delimiter! (})"
               call fatal
            end if
         end do
      end if
   end do
   if (.not. read_name) then
      write(*,*) "Please check the FFNAMES keyword in the QMDFF section!"
      call fatal
   end if


   if (.not. exist) then
      if (qmdffnumber .gt. 1) then
         write(*,*) "No ESHIFT keyword given in the QMDFF section!"
         write(*,*) " Please give the QMDFF shift energies in order to "
         write(*,*) " ensure a useful calculation!"
         call fatal
      end if
   end if
end if


!
!     initialize the single QMDFF´s
!
call prepare (fffile1,fffile2,fffile3,qmdffnumber)
nqmdff=qmdffnumber
!
!     So far only for one 
!
!
!-------2 QMDFFs-----------------------------------------------------------------
!
if (qmdffnumber.eq.2) then

   if (evb_de) then
!
!      Set default values 
!
      off_basis="1g"
      do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         call upcase (record)
         string = record(next:120)
         if (trim(adjustl(record(1:11))) .eq. 'DE_EVB {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'DE_EVB{') then

            do j=1,nkey_lines-i+1
               next=1
               record = keyline(i+j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               if (keyword(1:11) .eq. 'COUPLING ') then
                  read(record,*) names,off_basis
                  if (off_basis .eq. "1g" .or. off_basis.eq."2g" &
                      & .or. off_basis.eq."3g" &
                      & .or.off_basis.eq."sp" .or. off_basis.eq."sd" &
                      & .or. off_basis.eq."sd2" .or. off_basis.eq."sp2d3" &
                      & .or. off_basis.eq."sp2d"  ) then
                  else
                     if (rank .eq. 0) then
                        write(*,*) "No valid coupling function was defined."
                        write(*,*) "We will use the simple 1g-dE coupling instead."
                     end if
                     off_basis="1g"
                  end if
               end if
               if (keyword(1:11) .eq. '}') exit
               if (j .eq. nkey_lines-i) then
                  write(*,*) "The DE-EVB section has no second delimiter! (})"
                  call fatal
               end if
            end do
         end if
      end do
      goto 289
   end if

!
!      Read in the detailed settings for the dE-EVB coupling term
!      Loop over structure with brackets
!
   if (evb_dq) then
!
!      Set default values 
!
      off_basis="1g"
      do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         call upcase (record)
         string = record(next:120)
         if (trim(adjustl(record(1:11))) .eq. 'DQ_EVB {' .or. trim(adjustl(record(1:11))) &
              &  .eq. 'DQ_EVB{') then

            do j=1,nkey_lines-i+1
               next=1
               record = keyline(i+j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               if (keyword(1:11) .eq. 'COUPLING ') then
                  read(record,*) names,off_basis
                  if (off_basis .eq. "1g" .or. off_basis.eq."2g" &
                      & .or. off_basis.eq."sd2") then
                  else
                     if (rank .eq. 0) then
                        write(*,*) "No valid coupling function was defined."
                        write(*,*) "We will use the simple 1g-dQ coupling instead."
                     end if
                     off_basis="1g"
                  end if
               else if (keyword(1:11) .eq. 'TS_FILE ') then
                  read(record,*) names,filets
                  inquire(file=filets,exist=ts_xyz)
               end if


               if (keyword(1:11) .eq. '}') exit
               if (j .eq. nkey_lines-i) then
                  write(*,*) "The DE-EVB section has no second delimiter! (})"
                  call fatal
               end if
            end do
         end if
      end do
      goto 289
   end if



!
!     If one or both QMDFFs shall be corrected in the nonbonded part,
!     this will be activated and read in here
!
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'CORR_NONB ') then
         read(record,*) names,qmdff_index
         call qmdff_corr(qmdff_index)
      end if
   end do


   if (dg_evb) then
      if (rank .eq. 0) then
         write(*,*) "Used coupling method: DG-EVB"
      end if
!
!      Set default values 
!
      dg_evb_mode=0
      dg_evb_points=0
      dg_ref_file="grad_hess.dat"
      g_thres=1E-10
      do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'DG_EVB { ' .or. keyword(1:11) .eq. 'DG_EVB{ ') then
            do j=1,nkey_lines-i+1
               next=1
               record = keyline(i+j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               if (keyword(1:16) .eq. 'MODE ') then
                  read(record,*) names,dg_evb_mode
               else if (keyword(1:16) .eq. 'POINTS ') then
                  read(record,*) names,dg_evb_points
               else if (keyword(1:16) .eq. 'DOUBLE_ALPHA ') then
                  double_alpha=.true.
               else if (keyword(1:16) .eq. 'READ_COORD ') then
                  read_coord=.true.
               else if (keyword(1:16) .eq. 'POINTS_REF ') then
                  read(record,*) names,dg_ref_file
               else if (keyword(1:16) .eq. 'DG_MAT_NUM ') then
                  dg_mat_num=.true.
               else if (keyword(1:18) .eq. 'GAUSS_THRESHOLD ') then
                  read(record,*) names,g_thres
!
!     If the internal gradients for each structure shall be written in an extra 
!     file, for example to do a vector plot with 3d_evb.pl
!
               else if (keyword(1:16) .eq. 'INT_GRAD_PLOT ') then
                  int_grad_plot=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The INT_GRAD_PLOT option is activated! Internal gradient"
                     write(*,*) "components will be written to file int_grad.out!"
                  end if
!
!     If the internal coordinates for each structure shall be written in an extra 
!     file, for example to do a line plot with 3d_evb.pl
!
               else if  (keyword(1:16) .eq. 'INT_COORD_PLOT ') then
                  int_coord_plot=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The INT_COORD_PLOT option is activated! Internal coordinate"
                     write(*,*) "components will be written to file int_coord.out!"
                  end if

               end if
               if (record .eq. '}') exit
               if (j .eq. nkey_lines-i) then
                  write(*,*) "The DG-EVB section has no second delimiter! (})"
                  call fatal
               end if
            end do
         end if
      end do
      if (double_alpha) then
         add_alph=dg_evb_points
      end if
!
!     Catch missing information
!
      inquire(file=dg_ref_file,exist=exist)
      if (.not. exist) then
         if (rank .eq. 0) then
            dg_evb=.false.
            write(*,*) "The file ", trim(dg_ref_file)," containing the DG-ref. informations"
            write(*,*) "is not availiable!"
            call fatal
         end if
      end if
      if (dg_evb_mode .eq. 0) then
         write(*,*) "Please add the keyword 'MODE' to the DG_EVB section!"
         call fatal
      end if
      if (dg_evb_points .eq. 0) then
         write(*,*) "Please add the keywod 'POINTS' to the DG-EVB section!"
         call fatal
      end if

!
!     Now define set of internal coordinates: In this case always read in from file!
!
      read_coord=.true. 
      num_struc=1 ! dummy variable for number of structures for init_int; never used there
      int_mode=2
      if (read_coord) int_mode=1
      call init_int(filegeo,num_struc,rank,int_mode)
!
!   allocate needed arrays
!
      nat=natoms
      nat3=3*nat
      if (dg_evb_mode .eq. 1) then
         mat_size = dg_evb_points
         if (rank .eq. 0) then
            write(*,*) "Used DG-EVB mode: 1 - energies (E)"
         end if
         dg_mode=1
      else if (dg_evb_mode .eq. 2) then
         mat_size = dg_evb_points*(1+nat6)
         if (rank .eq. 0) then
            write(*,*) "Used DG-EVB mode: 2 - energies and gradients (E+G)"
         end if
         dg_mode=2
      else if (dg_evb_mode .eq. 3) then
         mat_size=dg_evb_points*(1+nat6+(nat6)*(nat6+1)/2)
         if (rank .eq. 0) then
            write(*,*) "Used DG-EVB mode: 3 - energies, gradients and hessians (E+G+H)"
         end if
         dg_mode=3
      else 
         if (rank .eq. 0) then
            write(*,*) "You have chosen the DG-EVB-Mode",dg_evb_mode
            write(*,*) "Only the modes 1 to 3 are availiable!"
            call fatal
         end if
      end if 
      if (rank .eq. 0) then
         write(*,*) "The dimension of the linear equation F=DB is ",mat_size
      end if
      nat=natoms
      allocate(xyz2(3,natoms))
      allocate(int_coord(nat6))
      allocate(geo_int(nat6))
      allocate(geo_xyz(3,nat))
      allocate(geo_xyz1(3*nat))
      allocate(all_xyz(nat*3,dg_evb_points))
      allocate(all_int(nat6,dg_evb_points))
      allocate(point_int(nat6,dg_evb_points))
      allocate(alph_opt(dg_evb_points))
      allocate(b_opt(mat_size))
      allocate(b_vec(mat_size))
      b_vec=b_opt
! 
!   Read in the needed reference structures
!   Use the same method as used for DG-EVB optimization
!
      allocate(all_ens(dg_evb_points))

      call read2int(1)
      point_int=all_int
!
!   Read in the evb_pars.dat-file and initialize the alpha 
!   and b_vec values
!
      open (unit=443,file="evb_pars.dat",status="old")
      
      read(443,*) names
      do i=1,mat_size
         read(443,*,iostat=readstat) b_vec(i)
         if (readstat .ne. 0) then
            if (rank .eq. 0) then
               write(*,*) "The evb_pars.dat file contains too few lines or has a wrong formate!"
               call fatal
            end if
         end if

      end do
      do i=1,dg_evb_points
         read(443,*,iostat=readstat) alph_opt(i)
         if (readstat .ne. 0) then
            if (rank .eq. 0) then
               write(*,*) "The evb_pars.dat file contains too few lines or has a wrong formate!"
               call fatal
            end if
         end if
      end do
      close(443)
!
!  Go back to main programs after printing PES information
!
      goto 678
   end if
!
!     read in the TREQ coupling term! No optimization was 
!     needed, so only reference data from the system will be read in and the 
!     interpolation will be initialized
!
   if (treq) then
      if (rank .eq. 0) then
         write(*,*) "Used coupling method: TREQ"
         write(*,*)
      end if
!
!   default values
!
      int_grad_plot=.false. 
      diff_ana=.false.
      read_coord=.false.
      dist_matrix=.false.
      interp_tol=0.0000000005d0
      rp_ana_step=0.0000001d0
      par_epsi=0.d0
      no_evb=.false.  ! if no RP-EVB shall be used but only direct interpolation and QMDFFs
      no_evb_xi=0.8d0  ! for NO EVB: maximum xi value that is calculated with pure QMDFF
      iq_l_lo=0.0000001d0
      iq_l_hi=0.0000001d0 
      iq_r_lo=0.9999999d0
      iq_r_hi=0.9999999d0
      rp_mid_tot=0.5d0   ! total region around TS (in s) where reference energies are interpolated
      rp_mid_trans=0.1d0  ! transition region between energies and EVB-QMDFF, inside rp_mid_tot
      pre_exp=4.d0  ! pre factor for exponential damping of RP-EVB
      pareta=20.d0  ! shape parameter for damping function
      path_dist_limit=100.d0
!      corr_max=1.0d0  ! maximal derivative of the QMDFF correction (Emax)
!      deltp=0.40  ! interval in which the QMDFF correction is applied
!      ddeltp=0.10d0  ! how much shorter the interval for correcting QMDFF hessians is 
                     ! this needs to be shorter because else divergence would occur
      rp_evb_mode=3   ! What kind of reference data shall be used (2: E+G, 3: E+G+H)
      s_ts=0.5d0  ! position of the transition state along the path
      s_bord_man=.false. ! if the positions of the different IRC areas shall be determined manually
      trans_l_lo=0.d0
      trans_l_hi=0.d0
      trans_r_lo=0.d0
      trans_r_hi=0.d0
      int_grad_plot=.false.
      int_coord_plot=.false.
      dg_ref_file="grad_hess.dat"
      file_irc_struc=""
      file_irc_ens=""
      irc_local=5
      rp_evb_points=0

      do i = 1, nkey_lines
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         call upcase (record)
         string = record(next:120)
         if (trim(adjustl(record(1:11))) .eq. 'TREQ { ' .or. & 
                     & trim(adjustl(record(1:11))) .eq. 'TREQ{ ') then
            do j=1,nkey_lines-i+1
               next=1
               record = keyline(i+j)
               call gettext (record,keyword,next)
               call upcase (keyword)
               if (keyword(1:16) .eq. 'POINTS ') then
                  read(record,*) names,rp_evb_points
               end if

               if (keyword(1:16) .eq. 'DIFF_ANA ') then
                  diff_ana=.true.
               end if
               if (keyword(1:16) .eq. 'READ_COORD ') then
                  read_coord=.true.
               end if
!
!     The tolerance: How small the interval along the IRC shall be 
!     that will determine the optimal s-value of a structure
!
               if (keyword(1:16) .eq. 'INTERP_TOL ') then
                  read(record,*) names,interp_tol
                  if (rank .eq. 0) then
                     write(*,*) "The maximal uncertainty for numerical calculation"
                     write(*,*) "of actual s values will be:",interp_tol
                  end if
               end if
!
!     Analytical RP_EVB gradient: so far, the "analytical" gradient 
!     needs to be calculated numerically in parts: the slope along 
!     the IRC path cannot be calculated analytically. Therefore it 
!     is done numerically in the space of internal coordinates
!
               if (.not. num_grad) then
                  if (keyword(1:16) .eq. 'RP_ANA_STEP ') then
                     read(record,*) names,rp_ana_step
                     if (rank .eq. 0) then
                        write(*,*) "The stepsize for numerical calculation of gradient"
                        write(*,*) "components along the IRC will be:",rp_ana_step
                     end if
                  end if
               end if
!
!     If the borders for the different areas along the IRC according 
!     to the RP-EVB paper shall be read in manually:
!     s_Q1-R, s_R-RI, s_RI-I, s_I-RI, s_RI-R (s_R_Q2 symmetric to s_Q1_R)
!     Nevertheless they will be corrected in order to coincide with 
!     spline interpolation borders 
!
               if (keyword(1:16) .eq. 'S_BORDERS_MAN ') then
                  read(record,*) names,s_q1_r,trans_l_lo,trans_l_hi,trans_r_lo,trans_r_hi
                  if ((trans_l_lo .eq. 0) .or. (trans_l_hi .eq. 0) .or. & 
                       & (trans_r_lo .eq. 0) .or. (trans_l_hi .eq. 0)) then
                     write(*,*) "ERROR! You activated the S_BORDERS_MAN option but one of "
                     write(*,*) "  the s-borders is still zero!"
                  end if
                  if (rank .eq. 0) then
                     write(*,*) "The S_BORDERS_MAN option was activated! Therefore the borders of"
                     write(*,*) " the different areas along the IRC are read in manually."
                  end if
                  s_bord_man=.true.
               end if
!
!     If no EVB shall be used, instead: only direct interpolation in the 
!     TS region and direct transition to pure QMDFFs at the asympotics
!
               if (keyword(1:16) .eq. 'NO_EVB ') then
                  no_evb=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The NO_EVB option was activated! Therefore no RP-EVB part will be"
                     write(*,*) "calculated: only direct interpolation and QMDFFs will be used, with"
                     write(*,*) "a smooth transition between them."
                  end if
                  read(record,*,iostat=readstat) names,rp_mid_tot,rp_mid_trans
                  if (readstat .ne. 0) then
                     write(*,*) "The keyword NO_EVB has the wrong line formate!"
                     call fatal
                  end if
               end if
!
!     For new method: near the TS the reference energies itself are 
!     interpolated to give a reasonable description of the PES in that 
!     region
!
               if (.not. no_evb) then
                  if (keyword(1:16) .eq. 'RP_EVB_MID ') then
                     read(record,*) names,rp_mid_tot,rp_mid_trans
                  end if
               end if
!
!     Smooth transition from full RP-EVB coupling at the transition path
!     to zero coupling at the Minima to avoid jumps 
!     calculation: s<=pi/pareta, convert it directly
!
               if (keyword(1:16) .eq. 'PAR_ETA ') then
                  read(record,*) names,pareta
          !  pareta=pi/pareta
               end if

!
!     Exponential coefficient of the exponential damping prefactor for 
!     the whole coupling term
!
               if (keyword(1:20) .eq. 'RP_EXP_COEFF') then
                  read(record,*) names,pre_exp
               end if
!
!     Determine if the QMDFF energy shifts shall be corrected automatically
!     to exactly reproduce the first/last energy of the path or not
!
               if (keyword(1:20) .eq. 'SHIFT_MANUAL ') then
                  shift_man=.true.
               end if
!
!     For acceleration of RP-EVB calculations: for s calculation, do not 
!     evaluate all reference IRC structures in each step but calculate 
!     only +/- irc_local steps, in reference to the previous timestep
!
               if (keyword(1:16) .eq. 'IRC_LOCAL ') then
                  read(record,*) names,irc_local
               end if
!
!     What kind of reference data shall be used for the taylor expansion
!     Default: Energies, gradients and hessians (mode 3)
!     Possible: Energies and gradients (mode 2)
!
               if (keyword(1:16) .eq. 'RP_EVB_MODE ') then
                  read(record,*) names,rp_evb_mode
               end if
!
!     If the internal gradients for each structure shall be written in an extra 
!     file, for example to do a vector plot with 3d_evb.pl
!
               if (keyword(1:16) .eq. 'INT_GRAD_PLOT ') then
                  int_grad_plot=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The INT_GRAD_PLOT option is activated! Internal gradient" 
                     write(*,*) "components will be written to file int_grad.out!"
                  end if
               end if
!
!     If the internal coordinates for each structure shall be written in an extra 
!     file, for example to do a line plot with 3d_evb.pl
!
               if (keyword(1:20) .eq. 'INT_COORD_PLOT ') then
                  int_coord_plot=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The INT_COORD_PLOT option is activated! Internal coordinate"
                     write(*,*) "components will be written to file int_coord.out!"
                  end if
               end if
!
!     Check if the file with the IRC structures is there
!
               if (keyword(1:11) .eq. 'IRC_STRUC ') then
                  read(record,*) names,file_irc_struc
                  inquire(file=file_irc_struc,exist=exist)
                  if (.not. exist) then
                     write(*,*) "The file", file_irc_struc, "containing the structures of the "
                     write(*,*) "reference path could not been found! Control the IRC_STRUC keyword!"
                     call fatal
                  end if
               end if
!
!     Check if the file with the IRC energies is there
!
               if (keyword(1:11) .eq. 'IRC_ENS ') then
                  read(record,*) names,file_irc_ens
                  inquire(file=file_irc_ens,exist=exist)
                  if (.not. exist) then
                     write(*,*) "The file", file_irc_ens, "containing the energies of the "
                     write(*,*) "reference path could not been found! Control the IRC_ENS keyword!"
                     call fatal
                  end if
!
!     Determine the number of lines in the energy file to get the number of 
!     structures in the IRC!
!
                  open(unit=56,file=file_irc_ens,status="old")
                  nlines=0
                  do
                     read(56,*,iostat=readstat)
                     if (readstat .ne. 0) exit
                     nlines=nlines+1
                  end do
                  close(56)
               end if
!
!     Maximum allowed z-value of the actual structure from the minimum path, above which it 
!     will be assumed that the finding of a useful limit failed and the structure is 
!     in reality far away from the path. QMDFF1 or 2 will be associated depending on the 
!     relative distances to the endpoints of the path   
!

               if (keyword(1:20) .eq. 'PATH_DIST_LIMIT ') then
                  read(record,*) names,path_dist_limit
                  if (rank .eq. 0) then
                     write(*,*) "The PATH_DIST_LIMIT keyword was found! Therefore a new tolerance value"
                     write(*,*) " for the maximum allowed z-value (internal distance from path) will be "
                     write(*,*) " stated. Above this value the s-value determination will be marked as"
                     write(*,*) " failed and either QMDFF1 or 2 will be used for calculation."
                     write(*,*) " The read in value is",path_dist_limit," (default: 100.0000000)"
                  end if
               end if
!
!     If the set of internal coordinates shall be defined a a simple distance 
!     matrix of all atoms in the system.
!     This should only be done if all other sets fail, because the number of 
!     coordinates scales quadratically with the systems size!
!
               if (keyword(1:20) .eq. 'DIST_MATRIX ') then
                  dist_matrix=.true.
                  if (rank .eq. 0) then
                     write(*,*) "The DIST_MATRIX option is activated! Set set of internal coordinates"
                     write(*,*) "will be defined as a simple distance matrix!"
                     write(*,*) "This should only be done if all other options fail because the"
                     write(*,*) "number of internals scales quadratically with the system size!"
                  end if
               end if

               if (keyword(1:13) .eq. '}') exit
               if (j .eq. nkey_lines-i) then
                  write(*,*) "The TREQ section has no second delimiter! (})"
                  call fatal
               end if
            end do
         end if
      end do

!
!     If manual borders of areas along IRC were read in, determine QMDFF damping area
!      
      if (s_bord_man) then
         pareta=pi/s_q1_r
      end if


!
!     Abort if missing/wrong information is given
!
      if (rp_evb_points .lt. 1) then
         write(*,*) "Too few gradient+Hessian reference points are given!"
         call fatal
      end if
      if (file_irc_struc .eq. "") then
         write(*,*) "Add the keyword IRC_STRUC!"
         call fatal
      end if

      if (file_irc_ens .eq. "") then
         write(*,*) "Add the keyword IRC_ENS!"
         call fatal
      end if

!
!     Call the initialization routine for this coupling type
!
      call rp_evb_init(file_irc_struc,file_irc_ens,nlines,rank)
!
!     go directly to the end of the subroutine
!
      goto 678
   end if
   289 continue
!    
!     Initialize the dE/dQ couplingterms and their parameters
!     Read them from the evb_pars.dat file obtained from EVB optimization
!  

   open (unit=443,file="evb_pars.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file evb_pars.dat containing the parameters for the EVB"
      write(*,*) " coupling term is not there! Please generate it using evbobt.x!"
      call fatal
   end if
   read(443,*,iostat=readstat) names
   read(443,*,iostat=readstat) offa
   read(443,*,iostat=readstat) offb
   if (off_basis .ne. "1g") then
      read(443,*,iostat=readstat) offc
      read(443,*,iostat=readstat) offd
   end if
   if (off_basis .eq. "2g") then
      read(443,*,iostat=readstat) offm
   else if (off_basis .eq. "3g") then
      read(443,*,iostat=readstat) offe
      read(443,*,iostat=readstat) offf
      read(443,*,iostat=readstat) offm
      read(443,*,iostat=readstat) offn
   else if (off_basis .eq. "sd2") then
      read(443,*,iostat=readstat) offe
      read(443,*,iostat=readstat) offf
   else if (off_basis .eq. "sp2d") then
      read(443,*,iostat=readstat) offe
      read(443,*,iostat=readstat) offf
      read(443,*,iostat=readstat) offg
      read(443,*,iostat=readstat) offh
   else if (off_basis .eq. "sp2d3") then
      read(443,*,iostat=readstat) offe
      read(443,*,iostat=readstat) offf
      read(443,*,iostat=readstat) offg
      read(443,*,iostat=readstat) offh
      read(443,*,iostat=readstat) offi
      read(443,*,iostat=readstat) offj
      read(443,*,iostat=readstat) offk
      read(443,*,iostat=readstat) offl
   end if
   close(443) 
   if (readstat .ne. 0) then
      write(*,*) "The evb_pars.dat file contains too few lines or has a wrong formate!"
      call fatal
   end if
!
!   The coupling with respect to the geometrical coordinate (dQ)
!   Read in the structure of a reference structure (mostly the TS)
!
   if (use_dq .eqv. .true.) then
      allocate(ts_coordinates_a(3,natoms))
      allocate(ts_coordinates(3*natoms-6))
      open(unit=33,file=filets,status='old')
      call next_geo(ts_coordinates_a,natoms,33,has_next)
      call xyz_2int(ts_coordinates_a,ts_coordinates,natoms)
      close(33)
   else 
   end if
   if (rank .eq. 0) then
      write(*,*) "The chosen coupling function: ",off_basis
   end if
end if

!
!-------3 QMDFFs-----------------------------------------------------------------
!

!
!     for 3x3-EVB: if numerical gradients should be used
!

if (qmdffnumber.eq.3) then
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'COUPLING ') then
         read(record,*) names,off_basis
         coupl=.true.
         if (off_basis .eq."1g" .or. off_basis .eq."sd2") then
         else
            off_basis="1g"
         end if
      end if
   end do
   if (.not.coupl) then
      off_basis="1g"
   end if
   full=.false.
   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'FULL-MATRIX ') then
         write(*,*) "The 1,3-EVB-matrix-elements are parameterized!"
         full=.true.
      end if
   end do
   if (.not. full) then
      write(*,*) "The 1,3-EVB-matrix-elements are set to zero."
   end if

   do i = 1, nkey_lines
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'EVB_PAR ') then
         if (off_basis.eq."1g") then
            read(record,*) names,offa,offb,offc,offd
            params=.true.
         else if (off_basis.eq."sd2") then
            read(record,*) names,offa,offb,offc,offd,offe,offf,offg,offh,offi,&
                           &offj,offk,offl,offm,offn
            params=.true.
         end if
      end if
   end do
!  
!     The coupling with respect to the geometrical coordinate (dQ)
!
   do i = 1, nkey_lines
       next = 1
       record = keyline(i)
       call gettext (record,keyword,next)
       call upcase (keyword)
       string = record(next:120)
       if (keyword(1:11) .eq. 'EVB_DQ ') then
          read(record,*) names,filets,filets2
          use_dq = .true.
       end if
   end do
   if (use_dq .eqv. .true.) then
      allocate(ts_coordinates(3*natoms-6))
      allocate(ts_coordinates2(3*natoms-6))
      allocate(ts_coordinates_a(3,natoms))
      allocate(ts_coordinates2_a(3,natoms))
      open(unit=33,file=filets,status='old')
      call next_geo(ts_coordinates_a,natoms,33,has_next)
      call xyz_2int(ts_coordinates_a,ts_coordinates,nat)
      close(33)
      open(unit=34,file=filets2,status='old')
      call next_geo(ts_coordinates2_a,natoms,34,has_next)
      call xyz_2int(ts_coordinates2_a,ts_coordinates2,nat)
      close(34)
      if (rank .eq. 0) then
         write(*,*) "The dQ-couplingterm is used. Basis:",off_basis
      end if
   end if
end if

!
!     From now on, perform general final setup
!
678 continue



!     
!    Set the periodicitiy internally if needed
!
!     For periodic systems: convert the box lengths to bohr as well!
!
if (periodic) then
   boxlen_x=boxlen_x/bohr
   boxlen_y=boxlen_y/bohr
   boxlen_z=boxlen_z/bohr

   call set_periodic(periodic,boxlen_x,boxlen_y,boxlen_z,ewald,zahn,ewald_brute,coul_cut,vdw_cut)
     ! set global variables in qmdff module..
end if
!
!     For systems with hard box walls: convert the box length to bohr
if (box_walls) then
   walldim_x=boxlen_x/bohr
   walldim_y=boxlen_y/bohr
   walldim_z=boxlen_z/bohr
end if

!
!     Finally, print out info concering the actual PES is usage 
!


if (rank .eq. 0) then
   write(*,*) "----------------------------------------------------------"
!
!     A single QMDFF
!
   if (qmdffnumber .eq. 1) then
      write(*,*) "PES description: A single QMDFF"
      write(*,*) " - QMDFF file: ",trim(fffile1)
      write(*,'(a,i8)') "  - Number of atoms: ",natoms
      if (periodic) then
         write(*,'(a)') "  - The system is simulated in a cubix periodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "       x=", boxlen_x*bohr," Ang.  ,y=",boxlen_y*bohr, &
                     & " Ang., z=",boxlen_z*bohr,"Ang."
         write(*,*) " - Long range Coulomb interactions treated with: ",trim(coul_method)
         write(*,'(a,f11.3,a)') "  - (Short) Coulomb interactions have a cutoff of ",coul_cut, " A"
         write(*,'(a,f11.3,a)') "  - VDW interactions have a cutoff of ",vdw_cut, " A"
      else if (box_walls) then
         write(*,'(a)') "  - The system is simulated in a hard-walls nonperiodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "       x=", boxlen_x," Ang.  ,y=",boxlen_y, &
                     & " Ang., z=",boxlen_z,"Ang."
      else if (coord_vasp) then
         write(*,'(a)') "  - The system is simulated in a periodic box given by POSCAR: "
         write(*,'(a,3f13.7)') "     a = ",vasp_a_vec(:)
         write(*,'(a,3f13.7)') "     b = ",vasp_b_vec(:)
         write(*,'(a,3f13.7)') "     c = ",vasp_c_vec(:)

      else 
         write(*,*) " - The simulated system has no periodicity."
      end if
      

   end if
!
!     One of the GFN-xTB methods 
!
   if (gfn_xtb) then
      write(*,*) "PES description: GFN-xTB semiempirics" 
      write(*,'(a,i8)') "  - Number of atoms: ",natoms
      if (periodic) then
         write(*,'(a)') "  - The system is simulated in a cubix periodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "       x=", boxlen_x*bohr," Ang.  ,y=",boxlen_y*bohr, &
                     & " Ang., z=",boxlen_z*bohr,"Ang."
      else if (box_walls) then
         write(*,'(a)') "  - The system is simulated in a hard-walls nonperiodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "       x=", boxlen_x," Ang.  ,y=",boxlen_y, &
                     & " Ang., z=",boxlen_z,"Ang."
      else if (coord_vasp) then
         write(*,'(a)') "  - The system is simulated in a periodic box given by POSCAR: "
         write(*,'(a,3f13.7)') "     a = ",vasp_a_vec(:)
         write(*,'(a,3f13.7)') "     b = ",vasp_b_vec(:)
         write(*,'(a,3f13.7)') "     c = ",vasp_c_vec(:)
      else
         write(*,*) " - The simulated system has no periodicity."
      end if
      write(*,'(a,f12.6)') "  - Total charge of the system: ",xtb_charge
      write(*,*) " - Used Hamiltonian: ",trim(hamil_string)

      write(*,'(a,e10.4)') "  - SCF energy convergence criterion (Hartrees): ", & 
            & xtb_accuracy*1D-6
      write(*,'(a,i6)') "  - Maximum number of SCF iterations: ",xtb_maxiter
      write(*,'(a,f16.6)') "  - Electronic temperature (Kelvin): ",xtb_el_temp/8.61732476104d-5
      if (solv_string .eq. "none") then  
         write(*,*) " - No implicit solvation will be used."
      else 
         write(*,*) " - Used solvation model: ",trim(solv_string)
         if (exist_spec) then
            write(*,*) " - Used solvent species: ",trim(solv_spec) 
         else 
            write(*,*) " - Dielectric constant used for solvation:",trim(solv_epsilon)
         end if
      end if
   end if
!
!     The pGFN-FF method by external  GULP call
!
   if (pgfn_ff) then
      write(*,*) "PES description: pGFN-FF force field (from GULP call)"
      write(*,'(a,i8)') "  - Number of atoms: ",natoms
      write(*,'(a,a)') "  - Initial setup from structure in: ",init_struc_file
      if (periodic) then
         write(*,'(a)') "  - The system is simulated in a cubix periodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "        x=", boxlen_x*bohr," Ang.  ,y=",boxlen_y*bohr, &
                     & " Ang., z=",boxlen_z*bohr,"Ang."
      else if (box_walls) then
         write(*,'(a)') "  - The system is simulated in a hard-walls nonperiodic box: "
         write(*,'(a,f13.7,a,f13.7,a,f13.7,a)') "       x=", boxlen_x," Ang.  ,y=",boxlen_y, &
                     & " Ang., z=",boxlen_z,"Ang."
      else
         write(*,*) " - The simulated system has no periodicity."
      end if

   end if

  
!
!     The SPC water model
!
   if (water_spc) then
      write(*,*) "PES desciption: Simple point charge (SPC) water model"
      write(*,'(a,i8)') "  - Number of atoms: ",natoms
      write(*,'(a,i8)') "  - Number of water molecules: ",natoms/3
   end if   

!
!     Inform the user if the numerical gradient is calculated 
!
   if (num_grad) then
      write(*,*) " - The gradient will be calculated numerically."
   end if
end if
if (rank .eq. 0) then
   write(*,*)
end if
return
end subroutine read_pes
