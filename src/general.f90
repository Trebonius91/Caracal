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
!     ##################################################################
!     ##                                                              ##
!     ##  module general  --  general parameters for the programs     ##
!     ##                                                              ##
!     ##################################################################
!
!     Here all needed global variables are listed which aren´t needed from
!     QMDFF or EVB-QMDFF directly but rather by file-IO, dynamics calculations
!     maximal sizes etc..
!
module general
implicit none
!   
!     This integer shows the several subroutines which program is ran in the moment
!     - 1: qmdffgen
!     - 2: evbopt
!     - 3: egrad
!     - 4: dynamic
!     - 5: evb_qmdff
!
integer::whichprog
!
!     from tinker module "files"
!     --> filenames etc.
integer::nprior
integer::ldir,leng
character(len=120)::filename
character(len=120)::outfile
!
!     from tinker module "keys"
!     --> array with command lines
!
integer::maxkey
parameter (maxkey=25000)
integer::nkey_lines
character(len=120)::keyline(maxkey)
!
!     from tinker module "ionuit" 
!     --> read/write units for input/output
!
integer::input
integer::iout
!
!     from tinker module "inform"
!     --> some general status variables..
!
integer::digits,iprint
integer::iwrite,isend
logical::debug,holdup,abort
!
!     from tinker module "argue"
!     --> number of command line arguments
!
integer::maxarg
parameter (maxarg=20)
integer::narg
logical::listarg(0:maxarg)
character(len=120)::arg(0:maxarg)
!
!     from tinker module "ascii"
!     --> ascii-numbers for a number of special characters
!
integer null,tab
integer::linefeed,formfeed,carriage,escape
integer::space,exclamation,quote,pound
integer::dollar,percent,ampersand,apostrophe
integer::asterisk,plus,comma,minus
integer::period,frontslash,colon,semicolon
integer::equal,question,atsign,backslash
integer::caret,underbar,vertical,tilde
parameter (null=0)
parameter (tab=9)
parameter (linefeed=10)
parameter (formfeed=12)
parameter (carriage=13)   !  carriage return
parameter (escape=27)
parameter (space=32)   !  blank space
parameter (exclamation=33)
parameter (quote=34)  !  double quote
parameter (pound=35)  !  pound sign
parameter (dollar=36)
parameter (percent=37)
parameter (ampersand=38)
parameter (apostrophe=39)  !  single quote
parameter (asterisk=42)
parameter (plus=43)
parameter (comma=44)
parameter (minus=45)
parameter (period=46)
parameter (frontslash=47)
parameter (colon=58)
parameter (semicolon=59)
parameter (equal=61)
parameter (question=63)
parameter (atsign=64)
parameter (backslash=92)
parameter (caret=94)
parameter (underbar=95)
parameter (vertical=124)  !  vertical bar
parameter (tilde=126)
!
!     from tinker module "sizes"
!     --> maximum number of different parameters/array sizes
!
integer::maxatm,maxtyp
parameter (maxatm=10000)
parameter (maxtyp=5000)
!
!     from tinker module "atoms"
!     --> coordinates for dynamics and optimizations
!
integer::n
integer::type(maxatm)  ! elements of the atoms
real(kind=8)::x(maxatm)
real(kind=8)::y(maxatm)
real(kind=8)::z(maxatm)
!
!     from tinker module "bath"
!     --> thermostate settings for dynamics!
!     maxnose     maximum length of Nose-Hoover thermostat chain
!     voltrial    mean number of steps between Monte Carlo moves
!     kelvin      target value for the system temperature (K)
!     atmsph      target value for the system pressure (atm)
!     tautemp     time constant for Berendsen thermostat (psec)
!     taupres     time constant for Berendsen barostat (psec)
!     compress    isothermal compressibility of medium (atm-1)
!     collide     collision frequency for Andersen thermostat
!     eta         velocity value for Bussi-Parrinello barostat
!     volmove     maximum volume move for Monte Carlo barostat (Ang**3)
!     vbar        velocity of log volume for Nose-Hoover barostat
!     qbar        mass of the volume for Nose-Hoover barostat
!     gbar        force for the volume for Nose-Hoover barostat
!     vnh         velocity of each chained Nose-Hoover thermostat
!     qnh         mass for each chained Nose-Hoover thermostat
!     gnh         force for each chained Nose-Hoover thermostat
!     isothermal  logical flag governing use of temperature control
!     isobaric    logical flag governing use of pressure control
!     anisotrop   logical flag governing use of anisotropic pressure
!     thermostat  choice of temperature control method to be used
!     barostat    choice of pressure control method to be used
!     volscale    choice of scaling method for Monte Carlo barostat
!
integer::maxnose
parameter (maxnose=4)
integer::voltrial
real(kind=8)::kelvin,atmsph
real(kind=8)::tautemp,taupres
real(kind=8)::compress,collide
real(kind=8)::eta,volmove
real(kind=8)::vbar,qbar,gbar
real(kind=8)::vnh(maxnose)
real(kind=8)::qnh(maxnose)
real(kind=8)::gnh(maxnose)
logical::isothermal
logical::isobaric
logical::anisotrop
character(len=9)::volscale

!
!     from tinker module "units"
!     --> different natural constants and conversation factors
!     literature references:
!
!     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended
!     Values of the Fundamental Physical Constants: 2010", Reviews of
!     Modern Physics, 84, 1527-1605 (2012)
!
!     Most values below are taken from 2010 CODATA reference values;
!     available on the web from the National Institute of Standards
!     and Technology at http://physics.nist.gov/constants/
!
!     The conversion from calorie to Joule is the definition of the
!     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)
!
!     The "coulomb" energy conversion factor is found by dimensional
!     analysis of Coulomb's Law, ie, by dividing the square of the
!     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
!     the permittivity of vacuum (the "electric constant"); note that
!     eps0 is typically given in F/m, equivalent to C**2/(J-m)
!
!     The approximate value used for the Debye, 3.33564 x 10-30 C-m,
!     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)
!
!     The value of "prescon" is based on definition of 1 atmosphere
!     as 101325 Pa set by the 10th Conference Generale des Poids et
!     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3
!
!     avogadro    Avogadro's number (N) in particles/mole
!     lightspd    speed of light in vacuum (c) in cm/ps
!     boltzmann   Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
!     gasconst    ideal gas constant (R) in kcal/mole/K
!     emass       mass of an electron in atomic mass units
!     planck      Planck's constant (h) in J-s
!     joule       conversion from calories to joules
!     convert     conversion from kcal to g*Ang**2/ps**2
!     bohr        conversion from Bohrs to Angstroms
!     hartree     conversion from Hartree to kcal/mole
!     evolt       conversion from Hartree to electron-volts
!     efreq       conversion from Hartree to cm-1
!     coulomb     conversion from electron**2/Ang to kcal/mole
!     debye       conversion from electron-Ang to Debyes
!     prescon     conversion from Hartree/bohr**3 to Atm
!
real(kind=8)::avogadro
real(kind=8)::lightspd
real(kind=8)::boltzmann
real(kind=8)::gasconst
real(kind=8)::emass,planck
real(kind=8)::joule,convert
real(kind=8)::bohr,hartree
real(kind=8)::evolt,efreq
real(kind=8)::coulomb,debye
real(kind=8)::prescon
real(kind=8)::newton2au
parameter (avogadro=6.02214179d+23)
parameter (lightspd=2.99792458d-2)
parameter (boltzmann=0.831446215d0)
parameter (gasconst=1.98720415d-3)
parameter (emass=5.485799095d-4)
parameter (planck=6.62606896d-34)
parameter (joule=4.1840d0)
parameter (convert=4.1840d+2)
parameter (bohr=0.52917721092d0)
parameter (hartree=627.5094743d0)
parameter (evolt=27.21138503d0)
parameter (efreq=2.194746313708d+5)
parameter (coulomb=332.063714d0)
parameter (debye=4.80321d0)
! parameter (prescon=737.395359668)
parameter (prescon=290371576.61)
parameter (newton2au=12137804.11081)
!
!     THERMOSTATS: Andersen or Nose-Hoover?
!
integer::thermostat  ! which thermostat (0: Andersen, 1: Nose-Hoover)
real(kind=8)::nose_zeta  ! the Nose-Hoover friction variable zeta
real(kind=8)::nose_q  ! the Nose-Hoover damping factor
real(kind=8),allocatable::nhc_zeta(:,:)  ! for NHC: array with thermostat DOFs
integer::nhc_length  ! length of the NHC thermostat chain
!
!     BAROSTAT
!
real(kind=8)::nose_tau  ! the Nose-Hoover barostat damping factor
!
!     from tinker module "atomid"
!     --> elements and masses of atoms in system
!
character(len=3)::name(maxatm)
!     Element indices of the atoms
integer,allocatable::elem_index(:)
real(kind=8)::mass(maxatm)
!
!     from tinker module "moldyn"
!     --> MD trajectory velocity and accelerations
!

real(kind=8),allocatable::v(:,:)
real(kind=8),allocatable::a(:,:)
real(kind=8),allocatable::aalt(:,:)

!
!     from tinker module "mdstuf"
!     --> MD trajectory control options
!

integer::irest
integer::nfree
logical::dorest
logical::velsave
logical::frcsave
logical::uindsave
!
!     For manual setting of random seed for andersen thermostat if desired
!
logical::seed_manual
integer::seed_value
!
!     For evaluation of coordinates during dynamic.x calculation
!
logical::eval_coord  ! if evaluation takes place
integer::eval_step  ! after how many MD steps a new evaluation shall be done 
integer::eval_number   ! number of evaluated coordinates 
integer,allocatable::eval_inds(:,:)   ! list of evaluated coordinates (indices)
!
!     from tinker module "usage"
!     --> used atoms for energy calculation: OBSOLETE
!
integer::nuse
integer,allocatable::iuse(:)
logical,allocatable::use(:)
!
!     from tinker module "molcul"
!     --> individual molecules in current system    
!
integer::nmol
integer,allocatable::imol(:,:)
integer,allocatable::kmol(:)
integer,allocatable::molcule(:)
real(kind=8),allocatable::molmass(:)
!
!     from subroutine getkey: to define it globally
!
character(len=120)::keyfile
!
!     If an analytical potential energy surface from the 
!     literature shall be used 
!
logical::pot_ana
character(len=80)::pot_type
!
!     For printout of GFN-xTB calculation details: number of calculations
!
integer::xtb_calc_num
!
!     For certain structure analysis tools: if the program explore is used
!
logical::use_explore
!
!     If PES topology analysis shall be done 
!
logical::pes_topol
character(len=3)::topol_bonds
real(kind=8)::topol_vdw_scale,topol_eht_cutoff
!
!     for MPI parallelization: tells, when it is used
!
logical::use_mpi
end module general


