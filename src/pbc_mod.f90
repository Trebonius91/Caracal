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
!     ##  module pbc_mod  --  parameters for periodic boundaries      ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "evb_mod" contains all important global informations for calculations 
!     on potential energy surfaces with periodic boundary conditions (PBC)
!

module pbc_mod
implicit none
save
!
!     Activates the periodic boundary conditions (always 3D!)
!
logical::periodic
!
!     Activates a box with hard walls
!
logical::box_nonp
!
!     Center of the box in cartesian space
!
real(kind=8)::box_center(3)
real(kind=8)::box_size
!
!     Length of the box (currently, only rectangular boxes allowed)
!
real(kind=8)::boxlen_x,boxlen_y,boxlen_z,boxlen_x2,boxlen_y2,boxlen_z2
real(kind=8)::box_vol  ! volume of the box

!
!     If a POSCAR file is given for geometry, assume VASP format for output
!
logical::coord_vasp
logical::vasp_selective   ! if selective dynamics is activated
logical::vasp_direct     ! if direct coordinates are given
character(len=2)::vasp_names(20)  ! the element names header section (VASP)
integer::vasp_numbers(20) ! the element numbers header section (VASP)
integer::nelems_vasp  ! number of elements in the header line
!
!     VASP unit cell vectors
!
real(kind=8)::vasp_scale
real(kind=8)::vasp_a_vec(3),vasp_b_vec(3),vasp_c_vec(3)
! for Hessian matrix and IR spectra calculations
logical::calc_freq_int
integer::int_incr
real(kind=8),allocatable::int_pos_vecs(:,:,:)  ! numerical elongations 
real(kind=8),allocatable::dip_list(:,:)  ! numerical dipole derivatives
!
!     Parameters for periodic potentials (cutoffs, Ewald, etc)
!
real(kind=8)::coul_cut,cut_low,vdw_cut  ! the value for the Coulomb/dispersion cutoff
logical::ewald,zahn ! which type for Ewald summation
logical::ewald_brute
real(kind=8)::zahn_a,zahn_acut,zahn_par  ! parameters for Zahn method
real(kind=8)::a_ewald,r_ew_cut
!real(kind=8)::a_ewald  ! the adjustable alpha parameter for Ewald sum
real(kind=8)::sqrt_p  ! the square root of the accuracy parameter
real(kind=8)::trtf  ! relative timings of real and reciprocal parts
real(kind=8)::ebuffer ! electrostatic buffering constant
real(kind=8)::volbox ! actual volume of the periodic simulation box
!   for the FFT in the Ewald sum
integer::nfft  ! number of grid points per dimension for reciprocal
integer::maxprime,maxtable
real(kind=8),allocatable::ffttable(:,:)
integer,allocatable::iprime(:,:)
integer::bsorder   ! the spline polynomial order 
real(kind=8),allocatable::bsmod1(:),bsmod2(:),bsmod3(:)
integer(kind=8)::planb,planf  ! for FFTW: pointer to backward/forward transform 
real(kind=8),allocatable::qgrid(:,:,:,:)  ! charge grid for reciprocal


end module pbc_mod
