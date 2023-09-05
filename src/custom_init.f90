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
!     subroutine custom_init: If a new PES has been added by the user, its 
!     initialization (if needed) is performed via this subroutine
!
!     part of EVB
!
subroutine custom_init
use general
use evb_mod
implicit none

if (cust_number .eq. 1) then
!
!     If an initialization is required, remove the comment for this call and 
!     name the routine as it is called into the provided file
!
!   call init_routine


else if (cust_number .eq. 2) then
!
!     If a second subroutine is added, add the call to its init routine here, if 
!     needed
!

else if (cust_number .eq. 3) then
!
!     For a third one ...
!
end if

end subroutine custom_init
