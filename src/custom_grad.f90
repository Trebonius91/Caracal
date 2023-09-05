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
!     subroutine custom_grad: Call a subroutine added by the user and 
!     calculate energy and gradient with it. If needed, numerical gradients 
!     can also be calculated. 
!
!     part of EVB
!
subroutine custom_grad (xyz2,e_evb,g_evb)
use general
use evb_mod
implicit none

integer::i,j,k
real(kind=8)::xyz2(3,natoms),e_evb,g_evb(3,natoms)
real(kind=8)::e_upper,e_lower,e_tmp,step

!
!     Caracal gives the current coordinates in bohr and expects the calculated 
!     energy to be in Hartrees and the gradient to be in Hartrees/bohr
!     If your subroutine uses other units, it might be helpful to apply one of the 
!     available global unit conversion parameters (simply multiply or divide 
!     the variables or the real numbers, directly):
! 
!     "bohr" (0.52917721092d0) : Bohr to Angstroms
!     "hartree" (627.5094743d0) : Hartree to kcal/mol
!     "evolt" (27.21138503d0): Hartree to eV
!     "joule" (4.1840d0) : kcal/mol to kJ/mol
!

if (cust_number .eq. 1) then
!
!     If the external routine is able to calculate the energy and analytical 
!     (or numerical gradients), activate this call here and adapt it to your 
!     routine!
!
!   call your_routine(xyz2,e_evb,g_evb)

!
!     If the external routine has no method of calculating analytical gradients,
!     add the keyword NUM_GRAD to the Caracal.key file and add the subroutine
!     call into the numerical gradient loop
!
   if (num_grad) then
!
!     Energgy of the undistorted structure
!     Activate this call here and adapt it to your routine!
!
!   call your_routine(xyz2,e_evb)

      step=num_grad_step
      do i=1,natoms
         do j=1,3
            do k=1,2
               if (k.eq.1) then
                  xyz2(j,i)=xyz2(j,i)+step
               else if (k.eq.2) then
                  xyz2(j,i)=xyz2(j,i)-2*step
               end if
!
!     Energies of the structure elongations for the gradient components
!     Activate this call here and adapt it to your routine!
!
!               call your_routine(xyz2,e_tmp)


               if (k.eq.1) then
                  e_upper=e_tmp
               else if (k.eq.2) then
                  e_lower=e_tmp
               end if
            end do
            g_evb(j,i)=(e_upper-e_lower)/(2*step)
            g_evb(j,i)=g_evb(j,i)
            xyz2(j,i)=xyz2(j,i)+step
         end do
      end do
   end if
else if (cust_number .eq. 2) then
!     
!     If you want to connect a second custom PES subroutine to Caracal, add it here!
!     It can then be called with PES_NUMBER 2
!     Copy the code from above
!

else if (cust_number .eq. 3) then
!
!     For a third custom PES and so on ...
!

end if

end subroutine custom_grad
