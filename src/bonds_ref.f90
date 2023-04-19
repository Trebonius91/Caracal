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
!     subroutine bonds_ref: calculate reference bondlengths 
!       of the forming and breaking bonds at the TS to define 
!       the s1 dividing surface for the rpmd rate calculation
!       program
!       For unimolecular reactions, also define reference 
!       bond lengths of the s0 educt surface 
!     
!     part of EVB
!
subroutine bonds_ref
use general
use evb_mod
implicit none
integer::i,j  ! loop indices
integer::atom1,atom2  ! the current atom numbers
real(kind=8)::Rx,Ry,Rz

allocate(form_ref(form_num),break_ref(break_num))
!
!     calculate the forming bonds
!
do i=1,form_num
   atom1=bond_form(i,1)
   atom2=bond_form(i,2)
   Rx=ts_ref(1,atom1)-ts_ref(1,atom2)
   Ry=ts_ref(2,atom1)-ts_ref(2,atom2) 
   Rz=ts_ref(3,atom1)-ts_ref(3,atom2)
   form_ref(i)=sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
end do

!
!     calculate the breaking bonds
!
do i=1,break_num
   atom1=bond_break(i,1)
   atom2=bond_break(i,2)
   Rx=ts_ref(1,atom1)-ts_ref(1,atom2)
   Ry=ts_ref(2,atom1)-ts_ref(2,atom2) 
   Rz=ts_ref(3,atom1)-ts_ref(3,atom2)
   break_ref(i)=sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
end do

!
!     For unimolecular reactions:
!     calculate the forming and breaking bonds 
!     at the transition state 
!
if ((umbr_type .eq. "CYCLOREVER") .or. (umbr_type .eq. "REARRANGE") .or.&
           &  (umbr_type .eq. "DECOM_1BOND") .or. (umbr_type .eq. "ELIMINATION")) then
   allocate(form_ed(form_num),break_ed(break_num))
   do i=1,form_num
      atom1=bond_form(i,1)
      atom2=bond_form(i,2)
      Rx=ed_ref(1,atom1)-ed_ref(1,atom2)
      Ry=ed_ref(2,atom1)-ed_ref(2,atom2)
      Rz=ed_ref(3,atom1)-ed_ref(3,atom2)
      form_ed(i)=sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
   end do

!
!     calculate the breaking bonds
!
   do i=1,break_num
      atom1=bond_break(i,1)
      atom2=bond_break(i,2)
      Rx=ed_ref(1,atom1)-ed_ref(1,atom2)
      Ry=ed_ref(2,atom1)-ed_ref(2,atom2)
      Rz=ed_ref(3,atom1)-ed_ref(3,atom2)
      break_ed(i)=sqrt(Rx*Rx+Ry*Ry+Rz*Rz)
   end do
end if
return
end subroutine bonds_ref
