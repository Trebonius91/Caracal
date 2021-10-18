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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!     subroutine read2int: Read in energies or energies and gradients or 
!     energies, gradients and hessians 
!     and convert to internal coordinates
!     It is used for the dg_evb initialization procedure
!     All infos are read in from ref.input
!     
!     part of EVB
!
subroutine read2int(dg_evb_mode) 
use evb_mod
implicit none
integer::atom_num
integer::state_open  ! if a file was opened successfully
character(len=20)::buffer
real(kind=8)::geo_int(nat6),grad_int(nat6)
real(kind=8)::geo_xyz(3,natoms),geo_xyz1(3*natoms)
real(kind=8)::hess_int(nat6,nat6)
real(kind=8)::energy
integer::i,j,dg_evb_mode
logical::redun
state_open=0
redun=.false.
!
!    For TREQ, define the number of points
!
if (treq) dg_evb_points=rp_evb_points

ref_input_unit=10
open(unit=ref_input_unit,file=dg_ref_file,status="old",iostat=state_open)
if (state_open .ne. 0) then
   write(*,*) "The file ",trim(dg_ref_file)," can´t be found! Exiting.."
   call fatal
end if
do i=1,dg_evb_points
!
!   first, read in the geometries and convert it to internals
!
   call read_geo(natoms,geo_int,geo_xyz1,geo_xyz,energy)
   all_ens(i)=energy
   all_xyz(:,i)=geo_xyz1
   all_int(:,i)=geo_int
!
!   then, convert the gradient to internals
!
!   write(117,*) "Structure No.",i
   if (dg_evb_mode .ge. 2) then
      call read_grad(natoms,grad_int,geo_int,geo_xyz)
      all_grad(:,i)=grad_int(:)
    
!      write(*,*) "grad_int reeaded!",grad_int
   end if 
!
!   last, read the hessian and convert it to internals
!  
   if (dg_evb_mode .eq. 3) then
      call read_hess(natoms,hess_int,geo_int,geo_xyz,grad_int)
      call project_hess(hess_int)
      all_hess(:,:,i)=hess_int(:,:)
   end if
end do
close(ref_input_unit)

return
end subroutine read2int
