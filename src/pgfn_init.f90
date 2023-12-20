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
!     subroutine pgfn_init: initialize the pGFN-FF force field for a certain
!     geometry (and from this, a specific bonding configuration)
!
!     part of EVB
! 
subroutine pgfn_init(natoms,xyz_init,names_init) 
#ifdef GULP
use datatypes
use m_gulp_interface
#endif
use pbc_mod
implicit none 
integer,intent(in)::natoms
real(kind=8),intent(in)::xyz_init(3,natoms)
character(len=2),intent(in)::names_init(natoms)
real(kind=8)::charges(natoms)
real(kind=8),allocatable::coord_cell(:,:)
character(len=4)::atomicsymbols(natoms)
character(len=80)::keywords,libraryfile
logical::lgulpoutput
integer::ndim
real(kind=8),allocatable::cell(:,:)
real(kind=8)::xyz(3,natoms)
!
!    Parameters to be transferred to the GULP library
!
!lgulpoutput          : logical      : If true then GULP will produce output; 
!                                        if false then the output will be suppressed.
!ndim                 : integer i4   : Indicates the number of periodic directions 
!                                        between 0 and 3
!                     :              : 0 => cluster
!                     :              : 1 => polymer (i.e. periodic in x direction)
!                     :              : 2 => slab    (i.e. periodic in x and y directions)
!                     :              : 3 => solid   (i.e. periodic in all directions)
!natoms               : integer i4   : The number of atoms in the system
!cell(3:ndim)         : real    dp   : The ndim Cartesian lattice vectors (in Angstroms)
!atomicsymbols(natoms): character    : The atom types in GULP format (i.e. element symbol, 
!                                        followed optionally by a number)
!xyz(3:natoms)        : real    dp   : The Cartesian coordinates of atoms (in Angstroms)
!charges(natoms)      : real    dp   : The charges of the atoms (in a.u.) unless set by 
!                                        the force field library file
!keywords             : character    : A string that contains the keywords required by
!                                        GULP (if any) NB: This can just be a blank string
!libraryfile          : character    : A string that contains the name of the library file
!                                         containing the force field (if any) NB: This 
!                                         can just be a blank string
!
lgulpoutput=.false.
!
!     Set the periodicity
!     Currently, only rectangular boxes are supported by Caracal
!     TODO: add support for arbitrary box shapes!
!
allocate(coord_cell(3,3))
coord_cell(1:3,1:3) = 0.0d0
if (periodic) then
   ndim=3
   coord_cell(1,1)=boxlen_x
   coord_cell(2,2)=boxlen_y
   coord_cell(3,3)=boxlen_z
else
   ndim=0
end if

!
!  Set atomic coordinates
!
xyz(:,:)=xyz_init(:,:)

atomicsymbols(:)=names_init(:)

charges(1:9) = 0.0d0

keywords="gfnff gwolf"
libraryfile=" "
#ifdef GULP
call init_gulp(lgulpoutput,ndim,natoms,coord_cell,names_init,xyz,charges,keywords,libraryfile)
#endif


return
end subroutine
