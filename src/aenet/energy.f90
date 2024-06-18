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
!     subroutine aenet_energy: calculate energy and forces for the 
!     current structure with an artificial neural network (ANN) 
!     potential. Taken from the aenet program
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

  subroutine aenet_energy(latticeVec, nAtoms, cooLatt, atomType, &
                        pbc, Ecoh, Etot, forCart, atomicEnergy)

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
                       recLattVec,                     &
                       geo_update_bounds,              &
                       origin,                         &
                       nTypes,                         &
                       atomTypeName                   

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

    double precision, dimension(3,3),                intent(in)  :: latticeVec
    integer,                                         intent(in)  :: nAtoms
    double precision, dimension(3,nAtoms),           intent(in)  :: cooLatt
    integer,          dimension(nAtoms),             intent(in)  :: atomType
    logical,                                         intent(in)  :: pbc
    double precision,                                intent(out) :: Ecoh
    double precision,                                intent(out) :: Etot
    double precision, dimension(3,nAtoms), optional, intent(out) :: forCart
    double precision, dimension(nAtoms),   optional, intent(out) :: atomicEnergy

#ifdef CHECK_FORCES
    double precision, dimension(3,nAtoms) :: forCart_num
    double precision                              :: E_i1, E_i2
    double precision :: d
    double precision, dimension(3,3) :: dd
    integer :: i, j
#endif

    logical                                       :: do_F, do_E_atom

    integer                                       :: nnb
    integer,          dimension(aenet_nnb_max)    :: nblist
    double precision, dimension(3,aenet_nnb_max)  :: nbcoo
    double precision, dimension(aenet_nnb_max)    :: nbdist
    integer,          dimension(aenet_nnb_max)    :: nbtype

    integer                                       :: type_i
    double precision, dimension(3)                :: coo_i
    double precision                              :: E_i

    integer                                       :: iatom, stat

    do_F = present(forCart)
    do_E_atom = present(atomicEnergy)
    if (do_F) forCart(1:3,1:nAtoms) = 0.0d0
    if (do_E_atom) atomicEnergy(1:nAtoms) = 0.0d0

    call lcl_init(aenet_Rc_min, aenet_Rc_max, latticeVec, nAtoms, &
                  atomType, cooLatt, pbc)

#ifdef CHECK_FORCES
    d = 0.01d0
    dd(:,1) = [d, 0.0d0, 0.0d0]
    dd(:,2) = [0.0d0, d, 0.0d0]
    dd(:,3) = [0.0d0, 0.0d0, d]
    forCart_num = 0.0d0
#endif

    Ecoh = 0.0d0
    Etot = 0.0d0
    atoms : do iatom = 1, nAtoms

       ! distribute atoms over processes:
       if (mod(iatom-1,ppSize) /= ppRank) cycle

       type_i = atomType(iatom)
       coo_i(1:3) = matmul(latticeVec, cooLatt(1:3,iatom))

       ! get all atoms of species type_i within the cut-off:
       nnb = aenet_nnb_max
       call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, aenet_Rc_max, &
                            nblist=nblist, nbtype=nbtype)

       if (do_F) then
          call aenet_atomic_energy_and_forces( &
               coo_i, type_i, iatom, nnb, nbcoo, nbtype, nblist, &
               nAtoms, E_i, forCart, stat)
#ifdef CHECK_FORCES
          do i = 1, 3
             coo_i = coo_i - dd(:,i)
             call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i1, stat)
             coo_i = coo_i + 2.0d0*dd(:,i)
             call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i2, stat)
             coo_i = coo_i - dd(:,i)
             forCart_num(i,iatom) = forCart_num(i,iatom) - (E_i2 - E_i1)/(2.0d0*d)
          end do
          do j = 1, nnb
             do i = 1, 3
                nbcoo(:,j) = nbcoo(:,j) - dd(:,i)
                call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                         E_i1, stat)
                nbcoo(:,j) = nbcoo(:,j) + 2.0d0*dd(:,i)
                call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                         E_i2, stat)
                nbcoo(:,j) = nbcoo(:,j) - dd(:,i)
                forCart_num(i,nblist(j)) = forCart_num(i,nblist(j)) - (E_i2 - E_i1)/(2.0d0*d)
             end do
          end do
#endif
       else
          call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i, stat)
       end if

       Etot = Etot + E_i
       Ecoh = Ecoh + E_i - aenet_free_atom_energy(type_i)
       if (do_E_atom) atomicEnergy(iatom) = E_i

    end do atoms

#ifdef CHECK_FORCES
    open(99, file='CHECK_FORCES.dat', status='replace', action='write')
    do iatom = 1, nAtoms
       write(99,'(9(1x,ES15.8))') &
            forCart(1:3,iatom), forCart_num(1:3,iatom), &
            forCart(1:3,iatom) - forCart_num(1:3,iatom)
    end do
    close(99)
#endif

    call lcl_final()

    call pp_sum(Ecoh)
    call pp_sum(Etot)

    ! gather results from all processes
    if (do_F) then
       if (ppSize>1) then
          do iatom = 1, nAtoms
             call pp_sum(forCart(1:3,iatom), 3)
          end do
       end if
    end if
    if (do_E_atom) call pp_sum(atomicEnergy(1:nAtoms), nAtoms)

  end subroutine aenet_energy
