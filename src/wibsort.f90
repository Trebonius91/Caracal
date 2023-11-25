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
! This file is part of xtb.
!     
! Copyright (C) 2017-2020 Stefan Grimme
!        
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!     
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!        
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.   
!     
!     Subroutine wibsort: Sorts Wiberg bond orders for an atom with ascending
!       order. Taken from the xtb program.

!        
!     part of tblite
!  
SUBROUTINE wibsort(ncent,imo,imem,qmo)
   implicit none
   integer  :: ncent
   integer  :: imo
   real(kind=8) :: qmo(ncent,ncent)
   integer  :: imem(ncent)
   integer  :: ii,i,j,k,ihilf
   real(kind=8) :: pp
   
   do ii = 2,ncent
      i = ii - 1
      k = i 
      pp= qmo(imo,i)
      do j = ii, ncent 
         if (qmo(imo,j) .lt. pp) cycle
         k = j
         pp=qmo(imo,j)
      enddo
      if (k .eq. i) cycle
      qmo(imo,k) = qmo(imo,i)
      qmo(imo,i) = pp
   
      ihilf=imem(i)
      imem(i)=imem(k)
      imem(k)=ihilf
   enddo
   
end SUBROUTINE wibsort

