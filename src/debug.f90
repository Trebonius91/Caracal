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
!     ##################################################################
!     ##                                                              ##
!     ##  module debug  --  arrays und variables for debug printing   ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "debug" contains variables used for error fixing or looking into
!     details of QMDFF and/or EVB energy/gradient calculations
!     It can be invoked by adding the keyword "DEBUG" into an egrad.x 
!     calculation (may be extended further..)
!
module debug
implicit none 
save
!     Shall the debugging mode be started?
logical::do_debug  
!     The actual energy component to be filled in the array next
real(kind=8)::act_part 
!     Array for QMDFF energy components along a reactionpath
real(kind=8),dimension(:,:,:),allocatable::qmdff_parts
!     Captions for the different energy components (labels)
character(len=50),dimension(:,:),allocatable::parts_labels
!     Loop indices for filling of qmdff_parts array (for 2 QMDFFs!)
integer::struc_no,comp_no,comp_no2
!     Current number of QMDFF (for H-bond/X-bond terms)
integer::qmdff_num
!     IO unit for debug trajectory (structures of equilibration)
integer::debug_unit 

end module debug
