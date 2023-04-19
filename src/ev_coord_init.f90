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
!     subroutine ev_coord_init: initialize evaluation of coordinate values 
!       during dynamic.x trajectory calculation
!
!     part of EVB
!
subroutine ev_coord_init
use general
implicit none 
integer::i
integer::readstat
integer::nlines
integer::i1,i2,i3,i4
character(len=80)::line_act

!
!    Read in number of lines/coordinates to evaluate 
!
open(unit=138,file="coord_eval.inp",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The EVAL_COORD keyword was found, but the file 'coord_eval.inp'"
   write(*,*) " is not there!"
   call fatal
end if
nlines=0
do 
   read(138,*,iostat=readstat) line_act
   if (readstat .ne. 0) exit   
   nlines=nlines+1
end do
eval_number = nlines
close(138)

allocate(eval_inds(eval_number,5))
eval_inds=0
!
!     Then read in all actual coordinate definitions 
!
open(unit=139,file="coord_eval.inp",status="old")
do i=1,eval_number
   read(139,'(a)') line_act
   read(line_act,*,iostat=readstat) i1,i2,i3,i4
   if (readstat .ne. 0) then
      read(line_act,*,iostat=readstat) i1,i2,i3
      if (readstat .ne. 0) then
         read(line_act,*,iostat=readstat) i1,i2
         if (readstat .ne. 0) then
            read(line_act,*,iostat=readstat) i1
            eval_inds(i,1)=1
            eval_inds(i,2)=i1
         else 
            eval_inds(i,1)=2
            eval_inds(i,2:3)=(/ i1,i2 /)
         end if
      else 
         eval_inds(i,1)=3
         eval_inds(i,2:4)=(/ i1,i2,i3 /)
      end if
   else 
      eval_inds(i,1)=4
      eval_inds(i,2:5)=(/ i1,i2,i3,i4 /)
   end if
end do
close(139)
!
!    Open the output file for coordinate printouts 
!     and print its header line 
!
open(unit=141,file="coord_eval.dat",status="replace")
write(141,'(a)') "# In this file, time-dependent values of different coordinates are listed."
 write(141,'(a)',advance="no") "# (MD step )    "
do i=1,eval_number
   if (eval_inds(i,1) .eq. 1) then
      write(141,'(a,i5,a)',advance="no") "(atom ",eval_inds(i,2),": x-coord,  y-coord, z-coord )    "
   else if (eval_inds(i,1) .eq. 2) then 
      write(141,'(a,i5,a,i5,a)',advance="no") "(bond ",eval_inds(i,2),"-",eval_inds(i,3)," )   "
   else if (eval_inds(i,1) .eq. 3) then
      write(141,'(a,i5,a,i5,a,i5,a)',advance="no") "(angle ",eval_inds(i,2),&
                & "-",eval_inds(i,3),"-",eval_inds(i,4)," )   "         
   else if (eval_inds(i,1) .eq. 4) then
      write(141,'(a,i5,a,i5,a,i5,a,i5,a)',advance="no") "(dihedral ",eval_inds(i,2),&
                &"-",eval_inds(i,3), "-",eval_inds(i,4),"-",eval_inds(i,5)," )    "
   end if
end do
write(141,*)

end subroutine ev_coord_init

