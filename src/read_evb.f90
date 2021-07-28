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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine read_evb: reads in QMDFFs generated from qmdffgen.x and EVB
!     coupling terms generated from evbopt.x and prepares upfollowing 
!     calculations done via egrad.x, dynamic.x, evb_qmdff.x (or rpmd.x)
!
!     part of EVB
! 
subroutine read_evb(rank)
use general
use evb_mod

implicit none
integer::num_arg,input_unit,i,qmdff_energies_unit,asciinum
character(len=70)::fffile1,fffile2,fffile3
real(kind=8)::energy,e_qmdff1,e_qmdff2,e_evb
integer::qmdffnumber,nat,nat3,l,m
integer::readstat  ! for error handling
integer::state_open  ! if a file was opened successfully
integer::nlines  ! number of structures in the IRC
integer::qmdff_index   ! number of QMDFF for corrections
character(len=120)::string
character(len=20)::keyword
character(len=120)::record
character(len=70)::fileinfo,filegeo
character(len=70)::file_irc_struc
character(len=70)::file_irc_ens
character(len=70)::filets,filets2,names
character(len=1)::qmdffnum
real(kind=8),dimension(:,:),allocatable::coord
real(kind=8),dimension(:,:),allocatable::g_evb
real(kind=8),dimension(:,:),allocatable::xyz2,geo_xyz
real(kind=8),dimension(:),allocatable::int_coord,geo_int,geo_xyz1
real(kind=8),dimension(:,:),allocatable::ts_coordinates_a,ts_coordinates2_a
integer::mode,next,j,k,readstatus,dg_evb_mode,mat_size,num_struc
integer::int_mode ! method of defining internal coordinates
real(kind=8)::s_q1_r  ! temporary variable for QMDFF damping range (RP-EVB)
logical::path_struc,path_energy,coupl,params
logical::evb1,evb2,evb3,ffname1,ffname2,ffname3,defqmdff
logical::exist,exists,has_next,coupl1
integer::rank ! the current MPI rank
real(kind=8)::pi

pi=3.141592653589793238462643d0

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
mueller_brown=.false.

!
!     First, check if an analytical potential energy surface is desired
!
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'POT_ANA ') then
      read(record,*) names,pot_type
      pot_ana = .true.
!
!    Print general info message about availiable potentials
!
      if (rank .eq. 0) then
         write(*,*) "You have invoked the POT_ANA option for analytical potentials!"
         write(*,*) "The following surfaces are availiable (ascending atom number)"
         write(*,*) "0)   Müller-Brown        (1 atom, 2-dimensional)"
         write(*,*) "1)   H2 + H    (H3)      (3 atoms: H, H, H)"
         write(*,*) "2)   BrH + H   (BrH2)    (3 atoms: H, Br, H)"
         write(*,*) "3)   O2 + O    (OH3)     (3 atoms: O, O, O)"
         write(*,*) "4)   OH2 + H   (OH3)     (4 atoms: O, H, H, H)"
         write(*,*) "5)   NH2 + HCl (ClNH2)   (5 atoms: H, N, H, H, Cl)"
         write(*,*) "6)   CH4 + H   (CH4H)    (6 atoms: H, C, H, H, H, H)"
         write(*,*) "7)   NH3 + OH  (NH3OH)   (6 atoms: H, N, H, H, O, H)"
         write(*,*) "8)   CH4 + OH  (CH4OH)   (7 atoms: H, C, H, H, H, O, H)"
         write(*,*) "9)   CH4 + CN  (CH4CN)   (7 atoms: H, C, H, H, H, C, N)"
         write(*,*) "10)   GeH4 + OH (GeH4OH)  (7 atoms: H, Ge, H, H, H, O, H)"
         write(*,*) "11)   C2H6 + H  (C2H7)    (9 atoms: C, H, H, H, C, H, H, H, H)"
        write(*,*)
      end if
      exit
   end if
end do

!
!     Check if ab-initio MD with orca shall be activated!
!
do i = 1, nkey
    next = 1
    record = keyline(i)
    call gettext (record,keyword,next)
    call upcase (keyword)
    string = record(next:120)
    if (keyword(1:11) .eq. 'ORCA ') then
       orca=.true.
       write(*,*) "The keyword orca was found!"
       write(*,*) "Ab-initio MD will be conducted!"
       write(*,*) "Input commands will be read in from 'orca_com.dat'"
       exit
    end if
end do
!
!     Read in orca command line!
!
if (orca) then
   open(unit=261,file="orca_com.dat",status="old")
   read(261,'(a)') orca_com
   close(261)
   return
end if


!
!     Initialize the analytical PES if requested, write additional PES infos 
!     to pot_info.out file
!
if (pot_ana) then
   open(unit=16,file="pot_info.out",status="unknown")
   write(16,*) "Informational output about the used analytical PES:"
   if (rank .eq. 0) then 
      write(*,*) "---------->"
   end if
   if (pot_type .eq. "mueller_brown") then
      if (rank .eq. 0) then
         write(15,*) "The Müller-Brown model surface for 2D systems will be used!"
         write(15,*) 
         write(*,*) "The Müller-Brown model surface for 2D systems will be used!"
         write(*,*)
      end if
      natoms=1
      nat6=2
      mueller_brown=.true.
      return
   else if (pot_type .eq. "h3") then
      if (rank .eq. 0) then
         write(15,*) "The H+H2 potential energy surface is used."
         write(15,*) "References:  D. G. Truhlar, C. J. Horowitz J. Chem. Phys.,"
         write(15,*) "             Vol. 68, p. 2466, 1978."
         write(15,*)
         write(*,*) "The H+H2 potential energy surface is used."
         write(*,*) "References:  D. G. Truhlar, C. J. Horowitz J. Chem. Phys.,"
         write(*,*) "             Vol. 68, p. 2466, 1978."
         write(*,*)
      end if
!     (no initialization needed)
      natoms=3  ! store number of atoms..
      return
   else if (pot_type .eq. "brh2") then
      if (rank .eq. 0) then
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
      natoms=3  ! store number of atoms..
!      stop "hhff"
      return 
   else if (pot_type .eq. "o3") then
      if (rank .eq. 0) then
         write(15,*) "The O2+O potential energy surface is used"
         write(15,*) "References:  Z. Varga, Y. Paukku, and D. G. Truhlar"
         write(15,*) "             J. Chem. Phys. 147, 154312/1-17 (2017) "
         write(15,*)
         write(15,*) "The O2+O potential energy surface is used"
         write(15,*) "References:  Z. Varga, Y. Paukku, and D. G. Truhlar"
         write(15,*) "             J. Chem. Phys. 147, 154312/1-17 (2017) "
         write(15,*)
      end if
      natoms=3    
      return 
   else if (pot_type .eq. "oh3") then
      if (rank .eq. 0) then
         write(15,*) "The OH2+H potential energy surface is used"
         write(15,*) "References:  G. C. Schatz and H. Elgersma"
         write(15,*) "             Chem. Phys. Lett. 73, 21 (1980) "
         write(15,*)
         write(*,*) "The OH2+H potential energy surface is used"
         write(*,*) "References:  G. C. Schatz and H. Elgersma"
         write(*,*) "             Chem. Phys. Lett. 73, 21 (1980) "
         write(*,*)
      end if
      call initialize_oh3
      natoms=4  ! store number of atoms..
!      stop "hhff" 
      return
   else if (pot_type .eq. "clnh3") then
      if (rank .eq. 0) then
         write(15,*) "The NH2+HCl potential energy surface is used"
         write(15,*) "References:  M. Monge-Palacios, C. Rangel, J.C. Corchado and J. Espinosa-Garcia"
         write(15,*) "             Int.J.Quantum.Chem. 112, 1887 (2012) "
         write(15,*)
         write(*,*) "The NH2+HCl potential energy surface is used"
         write(*,*) "References:   M. Monge-Palacios, C. Rangel, J.C. Corchado and J. Espinosa-Garcia"
         write(*,*) "              Int.J.Quantum.Chem. 112, 1887 (2012) "
         write(*,*)
      end if
      call initialize_clnh3
      natoms=5  ! store number of atoms..
!      stop "hhff" 
      return     
   else if (pot_type .eq. "ch4h") then
      if (rank .eq. 0) then
         write(15,*) "The CH4H potential energy surface is used."
         write(15,*) "Reference: J.C. Corchado, J.L. Bravo and J. Espinosa-Garcia "
         write(15,*) "             J.Chem.Phys.130,184314 (2009)."
         write(15,*)
         write(*,*) "The CH4H potential energy surface is used."
         write(*,*) "Reference: J.C. Corchado, J.L. Bravo and J. Espinosa-Garcia " 
         write(*,*) "             J.Chem.Phys.130,184314 (2009)."
         write(*,*)
      end if
      call initialize_ch4h  
      natoms=6   
      return    ! go back to main program directly
   else if (pot_type .eq. "nh3oh") then
      if (rank .eq. 0) then
         write(15,*) "The NH3OH potential energy surface is used."
         write(15,*) "References: M. Monge-Palacios, C. Rangel and J. Espinosa-Garcia "
         write(15,*) "           J.Chem.Phys. 138, 084305 (2013)."
         write(15,*)
         write(*,*) "The NH3OH potential energy surface is used."
         write(*,*) "References: M. Monge-Palacios, C. Rangel and J. Espinosa-Garcia "
         write(*,*) "           J.Chem.Phys. 138, 084305 (2013)."
         write(*,*)
      end if
      call initialize_nh3oh
      natoms=6
      return
   else if (pot_type .eq. "ch4oh") then
      if (rank .eq. 0) then
         write(15,*) "The CH4OH potential energy surface is used."
         write(15,*) "References:  J. Espinosa-Garcia, J. C. Corchado, J. Chem. Phys.,"
         write(15,*) "             Vol. 112, p. 5731, 2000."
         write(15,*)
         write(*,*) "The CH4OH potential energy surface is used."
         write(*,*) "References:  J. Espinosa-Garcia, J. C. Corchado, J. Chem. Phys.,"
         write(*,*) "             Vol. 112, p. 5731, 2000."
         write(*,*)
      end if
      call initialize_ch4oh
      natoms=7
      return    ! go back to main program directly
   else if (pot_type .eq. "ch4cn") then
      if (rank .eq. 0) then
         write(15,*) "The CH4CN potential energy surface is used."
         write(15,*) "References:  J. Espinosa-Garcia, C. Rangel and Y. V. Suleimanov,"
         write(15,*) "             Phys. Chem. Chem. Phys., 2017, 19, 19341."
         write(15,*)
         write(*,*) "The CH4CN potential energy surface is used."
         write(*,*) "References:  J. Espinosa-Garcia, C. Rangel and Y. V. Suleimanov,"
         write(*,*) "             Phys. Chem. Chem. Phys., 2017, 19, 19341."
         write(*,*)
      end if
      call initialize_ch4cn
      natoms=7
      return    ! go back to main program directly
   else if (pot_type .eq. "geh4oh") then
      if (rank .eq. 0) then
         write(15,*) "The GeH4+OH potential energy surface is used."
         write(15,*) "References:  J. Espinosa-Garcia, C. Rangel and J.C. Corchado"
         write(15,*) "             Phys. Chem. Chem. Phys. 18, 16941 (2016)"
         write(15,*) 
         write(*,*) "The GeH4+OH potential energy surface is used."
         write(*,*) "References:  J. Espinosa-Garcia, C. Rangel and J.C. Corchado"
         write(*,*) "             Phys. Chem. Chem. Phys. 18, 16941 (2016)"
         write(*,*)
      end if
      call initialize_geh4oh
      natoms=7  ! store number of atoms..
      return      
    else if (pot_type .eq. "c2h7") then
      if (rank .eq. 0) then
         write(15,*) "The C2H6+H potential energy surface is used."
         write(15,*) "References:  Arindam Chakraborty, Yan Zhao, Hai Lin, and Donald G. Truhlar"
         write(15,*) "             J. Chem. Phys., 124, 044315 (2006)."
         write(15,*)
         write(*,*) "The C2H6+H potential energy surface is used."
         write(*,*) "References:  Arindam Chakraborty, Yan Zhao, Hai Lin, and Donald G. Truhlar"
         write(*,*) "             J. Chem. Phys., 124, 044315 (2006)."
         write(*,*)
      end if
!     (no initialization needed)
      natoms=9  ! store number of atoms..
      return
   else 
      if (rank .eq. 0) then  
         write(*,*) "No valid potential energy surface was chosen!"
         write(*,*) "Choose a valid one (H3,BrH2,O3,OH3,CH4H,NH3OH,CH4OH,GeH4OH) or remove keyword"
         call fatal
      end if
   end if
   if (rank .eq. 0) then
      write(15,*) "Infos about potential initialization written to file pot_info.out."
      write(15,*)
      write(*,*) "Infos about potential initialization written to file pot_info.out."
      write(*,*)
   end if
   close(16)
end if
!
!     Read in the Force field Parameters
!
E_zero1=0d0; E_zero2=0d0; E_zero3=0d0
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'EQMDFF ') then
      read(record,*) names,E_zero1
      evb1=.true.
   else if (keyword(1:11) .eq. '2EVB ') then
      read(record,*) names,E_zero1,E_zero2
      evb2=.true.
   else if (keyword(1:11) .eq. '3EVB ') then
      read(record,*) names,E_zero1,E_zero2,E_zero3
      evb3=.true.
   end if
end do
qmdffnumber=0
fffile1=""; fffile2=""; fffile3=""
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'FFNAME ') then
      if (evb1) then
         read(record,*) names,fffile1
         ffname1=.true.
         defqmdff=.true.
         qmdffnumber=1
         nqmdff=1
      else if (evb2) then
         read(record,*) names,fffile1,fffile2
         ffname2=.true.
         defqmdff=.true.
         qmdffnumber=2
         nqmdff=2
      else if (evb3) then
         read(record,*) names,fffile1,fffile2,fffile3
         ffname3=.true.
         defqmdff=.true.
         qmdffnumber=3
         nqmdff=3
      end if
   end if
end do

! 
!     In the case that no QMDFF´s are defined in tinker.key,
!     read in dem manually 
!

if (.not. ffname1 .and. .not.  ffname2 .and. .not. ffname3) then
   exist=.false.
   do while (.not. exist)
      if (use_mpi) then
         if (rank .eq. 0) then
            write(*,*) "No valid keyfile with QMDFFs given!"
            call fatal
         end if
      else 
         write(iout,'(/," Number of QMDFF´s: ",$)')
      end if
      read (*,'(a1)')  qmdffnum
      asciinum = ICHAR(qmdffnum)
!
!     Convert vom ASCII to the real number (1=49,2=50,3=51)
!
      select case (asciinum)
      case (49)
         qmdffnumber=1
         nqmdff=1
      case (50)
         qmdffnumber=2
         nqmdff=2
      case (51)
         qmdffnumber=3
         nqmdff=3
      end select
      if (qmdffnumber.eq.1 .or. qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
         exist=.true.
         defqmdff=.true.
      end if
   end do
end if

if (.not.defqmdff) then
   exist=.false.
   do while (.not. exist)
      write(iout,'(/," Name of the first force field:  ",$)')
      read (*,'(a120)')  fffile1
      inquire(file=fffile1,exist=exist)
   end do
   if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,'(/," Name of the second force field:  ",$)')
         read (*,'(a120)')  fffile2
         inquire(file=fffile2,exist=exist)
      end do
   end if
   if (qmdffnumber.eq.3) then
      exist=.false.
      do while (.not. exist)
         write(iout,'(/," Name of the third force field:  ",$)')
         read (*,'(a120)')  fffile1
         inquire(file=fffile1,exist=exist)
      end do
   end if
end if
!
!     Read in the QMDFF-energys if not defined via keyfile
!
if ((.not. evb1) .and. (.not. evb2) .and. (.not. evb3)) then
   if (use_mpi) then
      if (rank .eq. 0) then
         write(*,*) "No QMDFF energies given in keyfile!"
         call fatal
      end if
   else 
      write(iout,'(/," QMDFF-Energy of the first QMDFF :  ",$)')
      read (*,*) E_zero1
      if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
         write(iout,'(/," QMDFF-Energy of the second QMDFF :  ",$)')
         read (*,*) E_zero2
      end if
      if (qmdffnumber.eq.3) then
         write(iout,'(/," QMDFF-Energy of the third QMDFF :  ",$)')
         read (*,*) E_zero3
      end if
   end if
end if
!
!     initialize the single QMDFF´s
!
call prepare (fffile1,fffile2,fffile3,qmdffnumber)
!
!     mainly for test reasons: activate numerical gradient
!
num_grad=.false.
full=.false.
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
   if (keyword(1:11) .eq. 'NUM_GRAD ') then
      num_grad=.true.
      if (rank .eq. 0) then
         write(*,*) "The numerical gradient will be calculated!"
      end if
   end if
end do
!
!     if numgrad is activated, ask for the step of elongations
!
if (num_grad) then
   num_grad_step=0.0001d0 ! default value
   do i = 1, nkey
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
!-------2 QMDFFs-----------------------------------------------------------------
!
if (qmdffnumber.eq.2) then
!
!     If one or both QMDFFs shall be corrected in the nonbonded part,
!     this will be activated and read in here
!
   do i = 1, nkey
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

!
!     This part for the Distributed Gaussian couplingterm
!     Jump to the gradient, if this keyword is noticed
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'DG_EVB ') then
         read(record,*) names,dg_evb_points
         dg_evb=.true.

         inquire(file="ref.input",exist=exist)
         if (.not. exist) then
            dg_evb=.false. 
            if (rank .eq. 0) then
               write(*,*) "You have ordered a Distributed Gaussian(DG)-EVB-calculation,"
               write(*,*) "but the files with the reference structures ref.input"
               write(*,*) "is not availiable!"
               write(*,*) "Look into the manual for further informations."
               call fatal
            end if
         end if
      end if
   end do
   
    

   if (dg_evb) then
      if (rank .eq. 0) then
         write(*,*) "Used coupling method: DG-EVB"
      end if
!
!   default values
!
      dg_evb_opt_multi=.false.
      maxstart=100
      diff_ana=.false.
      lower_bond=1
      upper_bond=30
      int_grad_plot=.false.
      int_coord_plot=.false.
      g_thres=1E-10
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:16) .eq. 'DG_EVB_OPT_MULTI ') then
            dg_evb_opt_multi=.true.
         end if
         if (keyword(1:16) .eq. 'START_POINTS ') then
            read(record,*) names,maxstart
         end if
         if (keyword(1:16) .eq. 'DIFF_ANA ') then
            diff_ana=.true.
         end if
         if (keyword(1:16) .eq. 'RANDOM_BONDS ') then
            read(record,*) names,lower_bond,upper_bond
         end if
         if (keyword(1:16) .eq. 'DG_EVB_MODE ') then
            read(record,*) names,dg_evb_mode
         end if
         if (keyword(1:16) .eq. 'READ_COORD ') then
            read_coord=.true.
         end if
         if (keyword(1:16) .eq. 'DOUBLE_ALPHA ') then
             double_alpha=.true.
             add_alph=dg_evb_points
         end if
         if (keyword(1:16) .eq. 'GAUSS_THRESHOLD ') then
            read(record,*) names,g_thres
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
         if (keyword(1:16) .eq. 'INT_COORD_PLOT ') then
            int_coord_plot=.true.
            if (rank .eq. 0) then
               write(*,*) "The INT_COORD_PLOT option is activated! Internal coordinate"
               write(*,*) "components will be written to file int_coord.out!"
            end if
         end if

      end do
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
!   Read in the evb.info-file and initialize the alpha 
!   and b_vec values
!
      open (unit=443,file="evb.info",status="old")
      
      read(443,*) names
      do i=1,mat_size
         read(443,*,iostat=readstat) b_vec(i)
         if (readstat .ne. 0) then
            if (rank .eq. 0) then
               write(*,*) "The evb.info file contains too few lines or has a wrong formate!"
               call fatal
            end if
         end if

      end do
      do i=1,dg_evb_points
         read(443,*,iostat=readstat) alph_opt(i)
         if (readstat .ne. 0) then
            if (rank .eq. 0) then
               write(*,*) "The evb.info file contains too few lines or has a wrong formate!"
               call fatal
            end if
         end if
      end do
      close(443)
!
!  Go back to main programs
!
      return
   end if
!
!     In case of a reaction path EVB (RP-EVB), check if the method shall be 
!     activated
!
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'RP_EVB ') then
         read(record,*) names,rp_evb_points
         rp_evb=.true.
         names="ref.input"
         inquire(file=names,exist=exist)
         if (.not. exist) then
            rp_evb=.false.
            write(*,*) "You have ordered a Reaction Path EVB calculation (RP-EVB),"
            write(*,*) "but the file ", trim(names)," containing the reference informations"
            write(*,*) "is not availiable!"
            write(*,*) "Look into the manual for further informations."
            call fatal
         end if
      end if
   end do 
!
!     read in the reaction path (RP)-EVB coupling term! No optimization was 
!     needed, so only reference data from the system will be read in and the 
!     interpolation will be initialized
!
   if (rp_evb) then
      if (rank .eq. 0) then
         write(*,*) "Used coupling method: RP-EVB"
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
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
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
         irc_local=5
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
      end do
!
!     If manual borders of areas along IRC were read in, determine QMDFF damping area
!      
      if (s_bord_man) then
         pareta=pi/s_q1_r
      end if

!
!     If the set of internal coordinates shall be defined a a simple distance 
!     matrix of all atoms in the system.
!     This should only be done if all other sets fail, because the number of 
!     coordinates scales quadratically with the systems size!
!
!     Only if read_coord is not activated!
!
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:20) .eq. 'DIST_MATRIX ') then
            if (.not. read_coord) then
               dist_matrix=.true.
               if (rank .eq. 0) then
                  write(*,*) "The DIST_MATRIX option is activated! Set set of internal coordinates"
                  write(*,*) "will be defined as a simple distance matrix!"
                  write(*,*) "This should only be done if all other options fail because the"
                  write(*,*) "number of internals scales quadratically with the system size!"
               end if
            end if
         end if
      end do

!
!     Maximum allowed z-value of the actual structure from the minimum path, above which it 
!     will be assumed that the finding of a useful limit failed and the structure is 
!     in reality far away from the path. QMDFF1 or 2 will be associated depending on the 
!     relative distances to the endpoints of the path   
!
path_dist_limit=100d0
do i = 1, nkey
   next = 1
   record = keyline(i)
   call gettext (record,keyword,next)
   call upcase (keyword)
   string = record(next:120)
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
end do

!
!      TEST: read in maximum xi value for which only QMDFF will be calculated
!
if (no_evb) then
   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:20) .eq. 'NO_EVB_XI ') then
         read(record,*) names,no_evb_xi
         write(*,*) "For the NO_EVB option, the minimum Xi value below it only QMDFFs shall"
         write(*,*) " be considered is set to",no_evb_xi
      end if
   end do
end if

!
!     For modification of the whole coupling: only the linear Taylor expansion term
!
      if (rank .eq. 0) then
         if (rp_evb_mode .eq. 2) then
            write(*,*) "The RP-EVB mode shall be 2! Therefore, only energies"
            write(*,*) "and gradients are used to extrapolate the path at the sides."
         else if (rp_evb_mode .eq. 1) then
            write(*,*) "You have chosen RP-EVB mode=1! This is not implemented!"
            call fatal
         end if
      end if
!
!     Check if the file with the IRC structures is there
!
      file_irc_struc=""
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'IRC_STRUC ') then
            read(record,*) names,file_irc_struc
            inquire(file=file_irc_struc,exist=exist)
            if (.not. exist) then
               write(*,*) "The file", file_irc_struc, "containing the structures of the "
               write(*,*) "reference path could not been found! Control the IRC_STRUC keyword!"
               call fatal
            end if
         end if
      end do
      if (file_irc_struc .eq. "") then
         write(*,*) "Add the keyword IRC_STRUC!"
         call fatal
      end if
!
!     Check if the file with the IRC energies is there
!
      file_irc_ens=""
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
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
      end do
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
      return
   end if
!
!     read in other coupling-terms: if there is no known term,
!     use the 1g-dE-coupling as default.
!

   do i = 1, nkey
      next = 1
      record = keyline(i)
      call gettext (record,keyword,next)
      call upcase (keyword)
      string = record(next:120)
      if (keyword(1:11) .eq. 'EVB_DQ ') then
         read(record,*) names,filets
         use_dq = .true.
      end if
   end do

   if (.not. dg_evb) then
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'COUPLING ') then
            read(record,*) names,off_basis
            coupl=.true.
            if (use_dq) then
                if ((off_basis .eq."1g") .or. (off_basis .eq."3g") .or. &
                    &(off_basis .eq."sd2")) then
                else
                   if (rank .eq. 0) then 
                      write(*,*) "You do not have submitted a valid coupling basis."
                      write(*,*) "Please use evbopt with a valid coupling basis first!"
                      call fatal
                   end if
                end if
            else 

 
            end if
               if ((off_basis.eq."1g") .or. (off_basis.eq."2g") &
                   & .or. (off_basis .eq."3g") &
                   & .or. (off_basis .eq."sp") .or. (off_basis .eq."sd") &
                   & .or. (off_basis .eq."sd2") .or. (off_basis .eq."sp2d") &
                   & .or. (off_basis .eq."sp2d3")) then
               else
                  if (rank .eq. 0) then
                     write(*,*) "You do not have submitted a valid coupling basis."
                     write(*,*) "Please use evbopt with a valid coupling basis first!"
                     call fatal
                  end if
               end if 
                 
         end if
      end do
      if (.not.coupl) then
         off_basis="1g"
      end if
     
!    
!     Initialize the dE/dQ couplingterms and their parameters
!     Read them from the evb.info file obtained from EVB optimization
!  

      open (unit=443,file="evb.info",status="old")
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
         write(*,*) "The evb.info file contains too few lines or has a wrong formate!"
         call fatal
      end if
!  
!     The coupling with respect to the geometrical coordinate (dQ)
!     Read in the structure of a reference structure (mostly the TS)
!
      if (use_dq .eqv. .true.) then
         allocate(ts_coordinates_a(3,natoms))
         allocate(ts_coordinates(3*natoms-6))
         open(unit=33,file=filets,status='old')
         call next_geo(ts_coordinates_a,natoms,33,has_next)
         call xyz_2int(ts_coordinates_a,ts_coordinates,natoms)
         close(33)
         if (rank .eq. 0) then
            write(*,*) "Used coupling method: EVB-dQ"
         end if
      else 
         if (rank .eq. 0) then
            write(*,*) "Used coupling method: EVB-dE"
         end if
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

else if (qmdffnumber.eq.3) then
   do i = 1, nkey
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
   do i = 1, nkey
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

   do i = 1, nkey
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
   do i = 1, nkey
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

return
end subroutine read_evb
