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
!     subroutine inter_uma: Dummy routines for the interface to UMA
!      over a C wrapper for a submission of input and output to Python
!      without file-IO      
!
!     part of EVB
!

module inter_uma

interface
!
!     The initialization routine
!
   subroutine init_uma(mlip_file,coord_file,task_name,mlip_len, &
               & coord_len,task_len) bind(C)
      use iso_c_binding
      character(kind=c_char),dimension(*)::mlip_file
      character(kind=c_char),dimension(*)::coord_file
      character(kind=c_char),dimension(*)::task_name
      integer(c_int),value::mlip_len
      integer(c_int),value::coord_len
      integer(c_int),value::task_len

   end subroutine init_uma

!
!     The energy+gradient calculation routine
!
   subroutine ase_uma(coords,unitcell,energy,gradient,natoms) bind(C)
      use iso_c_binding
      real(c_double),dimension(*)::coords
      real(c_double),dimension(*)::unitcell
      real(c_double),dimension(*)::gradient
      real(c_double)::energy
      integer(c_int)::natoms

   end subroutine ase_uma

end interface      

end module inter_uma
