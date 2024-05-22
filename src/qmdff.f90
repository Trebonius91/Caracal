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
!     ##  module qmdff  --  global variables used in QMDFF program    ##
!     ##                                                              ##
!     ##################################################################
!
!     All variables that were stated in common blocks in the f77 version
!     QMDFF are now defined in this module! 
!

module qmdff
implicit none
!     Pi and several exponents of it
real(kind=8)::pi,pi2,pi6,spi
parameter(pi=3.1415926535897932384626433832795029d0)
parameter(pi2=6.28318530717958623199592693708837d0)
parameter(pi6=961.389193575304212170602974626298d0)
parameter(spi=1.77245385090551599275151910313925d0)
!
!     variables from former ffparam-file (QMDFF1):
!
integer::ndim,ntterm
parameter(ndim=500000)
parameter(ntterm=4)
integer::bond(2,ndim),angl(3,ndim),tors(6,ndim),hb(3,ndim),nci(3,20*ndim)
real(kind=8)::vbond(3,ndim),vangl(2,ndim),vtors(3*ntterm+2,ndim)
real(kind=8)::vhb(2,ndim)
integer::nbond,nbond12,nangl,ntors,nhb,nnci
integer::nmols ! number of molecules in the QMDFF
real(kind=8)::molnum(ndim)
real(kind=8),allocatable::q_glob(:)  ! global charge array for solvent QMDFF
real(kind=8)::scalehb_glob(94),scalexb_glob(94)
!
!     variables from former ffparam_two file (QMDFF2):
!
integer::bond_two(2,ndim),angl_two(3,ndim),tors_two(6,ndim),hb_two(3,ndim)
integer::nci_two(3,20*ndim)
real(kind=8)::vbond_two(3,ndim),vangl_two(2,ndim),vtors_two(3*ntterm+2,ndim)
real(kind=8)::vhb_two(2,ndim)
integer::nbond_two,nbond12_two,nangl_two,ntors_two,nhb_two,nnci_two
!
!     variables from former ffparam_three file (QMDFF3):
!
integer::bond_three(2,ndim),angl_three(3,ndim),tors_three(6,ndim),hb_three(3,ndim)
integer::nci_three(3,20*ndim)
real(kind=8)::vbond_three(3,ndim),vangl_three(2,ndim),vtors_three(3*ntterm+2,ndim)
real(kind=8)::vhb_three(2,ndim)
integer::nbond_three,nbond12_three,nangl_three,ntors_three,nhb_three,nnci_three

!
!     for qmdffgen input: used software package for reference
!
logical::read_software ! if the information is given into the keyfile
!  kind of used reference software (O(orca), G(gaussian), T(turbomole),C (CP2K))
character(len=1)::software

!     variables from former cutoffs-file
!     from common cutoffs
real(kind=8)::athr  !180 - athr or athr = linear or planar inversion

!     from common dftd3
real(kind=8)::c6ab(94,94,5,5,3),rcov(94)
real(kind=8)::r2r4(94),r094(94,94),a1,a2,s8
integer::maxci(94)
!     from common restrainpot
integer::nrs,rest(2,10000)
real(kind=8)::vrest(2,10000)

!     from common rotplot
!     --> set values > 0 (atom numbers) to get the rot plot
integer::nrot1,nrot2

!     from common hscal
!     --> scale factor 
real(kind=8)::scalh
!     from common screenes
!     -->
real(kind=8)::eps1(6)
!     from common screenrd
!     -->
real(kind=8)::eps2(6)
!     from common screenrd
!     -->
real(kind=8)::eps3(6)
!     from common atomradii
!     --> atomic radii of all elements!!
real(kind=8)::rad(94)
!
!     COVALENT RADII of all elements
!     based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!     in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!     edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!     corrected Nov. 17, 2010 for the 92nd edition.
!
!data  rad / &
! & 0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0, &
! & 0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0, &
! & 2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0, &
! & 1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0, &
! & 2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0, &
! & 1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0, &
! & 2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0, &
! & 1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0, &
! & 1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0, &
! & 1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0, &
! & 1.84D0,1.83D0,1.80D0,1.80D0 /

data rad/0.699,0.605,2.532,1.701,1.549,1.455,1.417,1.379,1.342, &
     &       1.304,2.910,2.457,2.230,2.098,2.003,1.928,1.871, &
     &       1.833,3.704,3.288,2.721,2.570,2.362,2.499,2.627, &
     &       2.362,2.381,2.286,2.607,2.475,2.381,2.305,2.248, &
     &       2.192,2.154,2.079,3.987,3.628,3.061,2.797,2.589, &
     &       2.740,2.948,2.381,2.551,2.476,2.891,2.797,2.721, &
     &       2.665,2.608,2.551,2.513,2.456,4.251,3.742,3.194, &
     &       0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000, &
     &       0.000,0.000,0.000,0.000,0.000,3.024,2.835,2.608, &
     &       2.759,3.005,2.419,2.590,2.419,2.721,2.815,2.797, & 
     &       2.778,2.759,0.000,0.000,2.740,0.000,0.000,0.000, &
     &       0.000,0,000,0.000,0.000/ 
!     from common valelec
!     --> number of valence electrons for all elements
integer::z(94)
!     from common zeatb
!     --> Hermans Slater exponents for EHT calculations
real(kind=8)::zet(86,3)
!real(kind=8)::parz(86,3) ! for subroutine basis etc.
!     from common voiptb
!     --> ionization potentials (?)
real(kind=8)::ip(86,3)
real(kind=8)::parip(86,3)  ! for subroutine basis etc.
!     from common ehtcommon in ehtcommon.f
!     --> parameter for Extended Hückel calculations
integer::maxao
parameter (maxao=20000)
integer::lao(maxao),nprim(maxao),aoat(maxao),fila(2,maxao/10)
integer::nshell(maxao/10),shell(5,maxao/10)
real(kind=8)::alp(maxao*3),cont(maxao*3),hdiag(maxao)
!     from common metalatoms
!     --> which atoms are metals?
integer::metal(94)
!     from common atomEN
!     --> electronegativity values for all elements
real(kind=8)::en(94)
data en/ 2.300,4.160,0.912,1.576,2.051,2.544,3.066,3.610,4.193  &
     &        ,4.789,0.869,1.293,1.613,1.916,2.253,2.589,2.869,3.242 &
     &        ,0.734,1.034,1.19, 1.38 ,1.53 ,1.65 ,1.75 ,1.80 ,1.84 &
     &        ,1.88 ,1.85 ,1.59, 1.756,1.994,2.211,2.434,2.685,2.966 &
     &        ,0.706,0.963,1.12, 1.32 ,1.41 ,1.47 ,1.51 ,1.54 ,1.56 &
     &        ,1.59 ,1.87 ,1.52, 1.656,1.824,1.984,2.158,2.359,2.582 &
     &        ,0.659,0.881,1.09, 1.16 ,1.34 ,1.47 ,1.60 ,1.65 ,1.68 &
     &        ,1.72 ,1.92 ,1.76, 1.789,1.854,2.01 ,2.19 ,2.39 ,2.60 &
     &        ,0.67 ,0.89, 20*2.0/

!     from common lenjonpar
!     --> Lennard-Jones parameter for all elements
real(kind=8)::ade(94),ade13(94)
!     from common screenovlp
!     --> For EVT torsion profile calculation
real(kind=8)::sscal(6)
!
!     from common atomAT
!     --> For thermochemical properties
!
real(kind=8)::eheat(94)
!
!     correction of QMDFFs inducted from CORR_NONB, in order to
!     archive better results in k(T) calculations
!
integer::corr_nonb
integer::corr_number
integer,allocatable::corr_at1(:),corr_at2(:)
real(kind=8),allocatable::corr_factor(:)

!    for detailed printout of QMDFF and EVB parameters etc.
logical::details
!    If in the case of VASP reference selective dynamics has been used and parts 
!    of the hessian are zero, print a warning
logical::vasp_hessian_sel
!    If qmdffgen shall not generate a QMDFF but only determine coordinates based 
!    on Wiberg-Mayer bond orders and print out Wilson matrix/derivatives
logical::check_coord

logical::ff_mod_noncov   ! If the noncovalent QMDFF parameter were optimized separately
real(kind=8),allocatable::mn_par(:)  ! the parameters for the noncovalent optimization

real(kind=8)::vir_ten(3,3)  ! the virial tensor for barostat 
logical::calc_vir   ! if the virial shall be calculated at all
end module qmdff

