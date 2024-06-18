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
!     subroutine aenet_initialize: initialize an artificial 
!     neural network (ANN) potential. Taken from the aenet program
!
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!+ ---------------------------------------------------------------------
!+ If you make use of AENET for your publication, please cite:
!+ [1] N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
!+ [2] J. Behler and M. Parrinello, Phys. Rev. Lett. 98 (2007) 146401.
!+
!+ If you used the Chebyshev descriptor, please cite:
!+ [3] N. Artrith, A. Urban, and G. Ceder, PRB 96 (2017) 014112.
!-----------------------------------------------------------------------
! 2011-11-17 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
!


  subroutine aenet_initialize(ann_elnum,ann_files,ann_elements)


  use aenet_aeio,      only: aeio_header,                    &
                       aeio_timestamp,                 &
                       aeio_print_copyright

  use aenet_aenet,     only: aenet_init,                     &
                       aenet_final,                    &
                       aenet_atomic_energy,            &
                       aenet_atomic_energy_and_forces, &
                       aenet_convert_atom_types,       &
                       aenet_free_atom_energy,         &
                       aenet_load_potential,           &
                       aenet_print_info,               &
                       aenet_Rc_min, aenet_Rc_max,     &
                       aenet_nnb_max

  use aenet_constants, only: PI

  use aenet_geometry,  only: geo_init,                       &
                       geo_final,                      &
                       pbc,                            &
                       latticeVec,                     &
                       recLattVec,                     &
                       geo_update_bounds,              &
                       origin,                         &
                       nAtoms,                         &
                       nTypes,                         &
                       atomType,                       &
                       atomTypeName,                   &
                       cooLatt

  use aenet_input,     only: InputData,                      &
                       read_InpPredict

  use aenet_io,        only: io_adjustl

  use aenet_lclist,    only: lcl_init,                       &
                       lcl_final,                      &
                       lcl_nmax_nbdist,                &
                       lcl_nbdist_cart

  use aenet_optimize,  only: opt_init,                       &
                       opt_final,                      &
                       opt_optimize_coords

  use aenet_parallel,  only: pp_init,                        &
                       pp_final,                       &
                       pp_bcast,                       &
                       pp_bcast_coo,                   &
                       pp_print_info,                  &
                       pp_bcast_InputData,             &
                       pp_bcast_latt,                  &
                       pp_sum,                         &
                       ppMaster, ppRank, ppSize
    implicit none

    integer::ann_elnum  ! number of elements
    character(len=50)::ann_files(50)  ! ANN file names
    character(len=2)::ann_elements(50) ! associated elements


    logical :: fexists
    integer :: nargs
    integer :: stat
    integer :: itype

    call pp_init()

    ! broadcast the information to all MPI processes
    call pp_bcast(ann_elnum)
    do itype=1,ann_elnum
       call pp_bcast(ann_files(itype))
       call pp_bcast(ann_elements(itype))
    end do


    ! initialize aenet (take names of elements)
    call aenet_init(ann_elements(1:ann_elnum), stat)
    if (stat /= 0) then
       write(*,*) 'Error: aenet initialization failed'
       call fatal
       stop
    end if

    ! load ANN potentials (take filenames of ANNs)
    do itype = 1, ann_elnum
       call aenet_load_potential(itype, ann_files(itype), stat)
       if (stat /= 0) then
       write(*,*) 'Error: could not load ANN potentials'
          call fatal
          stop
       end if
    end do

!    if (ppMaster .and. (inp%verbosity > 0)) then
!       ! write header and copyright info
!       call aeio_header("Atomic Energy Network Interpolation", char='=')
!       call aeio_header(aeio_timestamp(), char=' ')
!       write(*,*)
!       call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')
!    end if

  end subroutine aenet_initialize
