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
!     subroutine egrad_aenet: Calculate energy and gradients with an
!      artificial network in the style of the aenet program
!
!     part of EVB
!

subroutine egrad_aenet(xyz2,pot_grad,e_evb)
use evb_mod 
use general
use pbc_mod
implicit none
real(kind=8)::xyz2(3,natoms)
real(kind=8)::xyz_direct(3,natoms)
real(kind=8)::pot_grad(3,natoms,1)
real(kind=8)::e_evb
real(kind=8)::e_coh
real(kind=8)::coord_mat(3,3)
real(kind=8)::coord_mat_inv(3,3)
real(kind=8)::q_act_frac
real(kind=8)::e_per_atom(natoms)
real(kind=8)::xyz_frac_act(3)
character(len=2)::name_act
integer::i,j

!
!     For periodic systems (currently only cubic): initialize it
!
!     If the input was read from a POSCAR file, apply the VASP symmetry directly!
!
if (coord_vasp) then
   periodic = .true.
   coord_mat(:,1) = vasp_a_vec
   coord_mat(:,2) = vasp_b_vec
   coord_mat(:,3) = vasp_c_vec
else
   if (periodic) then
      coord_mat = 0.d0
      coord_mat(1,1)= boxlen_x
      coord_mat(2,2) = boxlen_y
      coord_mat(3,3) = boxlen_z
   else
      coord_mat = 0.d0
      coord_mat(1,1)= 100.0
      coord_mat(2,2) = 100.0
      coord_mat(3,3) = 100.0
   end if
end if

!
!    Convert the coordinates to direct coordinates!
!
call matinv3(coord_mat/bohr,coord_mat_inv)

do i=1,natoms
   xyz_frac_act=matmul(coord_mat_inv,xyz2(:,i))
   xyz_direct(:,i)=xyz_frac_act
end do
!
!    Determine all element indices according to the ANN ordering  
!

if (.not. allocated(at_type_ann)) then
   allocate(at_type_ann(natoms))
   do i=1,natoms
      do j=1,ann_elnum
         if (name(i) .eq. ann_elements(j)) then
            at_type_ann(i)=j
            exit
         end if 
      end do
   end do
end if

!
!    Call the aenet energy calculation routine
!
call aenet_energy(coord_mat,natoms,xyz_direct,at_type_ann,periodic,&
       & e_coh,e_evb,pot_grad,e_per_atom)

!
!    Convert energy and gradient vector to atomic coordinates
!
e_evb=e_evb/evolt
pot_grad=-pot_grad/evolt*bohr

return

end subroutine egrad_aenet

