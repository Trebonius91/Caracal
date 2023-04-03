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
!     ##################################################################
!     ##                                                              ##
!     ##  module h2co_mod  --  parameters for H2CO analytical PES     ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "h2co_mod" contains all important global informations that are
!     used within the H2CO potential energy surface 
!

module h2co_mod
implicit none 
save 

real,parameter::total=8
integer,parameter::coff=1561
integer,parameter ::wp=selected_real_kind(12,300)
!real,parameter::alpha(4)=(/2.0,2.0,2.0,2.0/)
real(kind=8)::cof(coff),power_read(6,coff),power(coff,6)



end module h2co_mod

