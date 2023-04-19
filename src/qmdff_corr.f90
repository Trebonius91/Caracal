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
!     Subroutine qmdff_corr: Correction of QMDFF parameters to 
!       ensure a better asymptotic behavior in k(T) calculations
!
!     part of QMDFF
!
subroutine qmdff_corr(qmdff_index,rank)
use qmdff 
implicit none 
integer::qmdff_index
integer::i
integer::rank
integer::readstat
integer::corr_at1_loc(100),corr_at2_loc(100)
real(kind=8)::corr_factor_loc(100)
logical::existing
!
!     activate the global QMDFF correction flag :
!     Store the global index with the number of QMDFF (1 or 2) that shall be 
!     corrected 
!
corr_nonb=qmdff_index

!
!     Read in the file with detailed informations about corrections
!
inquire(file="qmdff_corr.dat",exist=existing)
if (.not. existing) then
   write(*,*) "You have chosen the option CORR_NONB, but the file"
   write(*,*) "qmdff_corr.dat with additional infos could not been found!"
   call fatal
end if
if (rank .eq. 0) then
   write(*,*) "You have called the CORR_NONB option! Therefore, nonbonded terms"
   write(*,*) "of QMDFF No.",qmdff_index," will be corrected regarding the commands"
   write(*,*) "in the file qmdff_corr.dat."
end if

i=0
open(unit=97,file="qmdff_corr.dat",status="old")
do
   i=i+1
   read(97,*,iostat=readstat) corr_at1_loc(i),corr_at2_loc(i),corr_factor_loc(i)
   if (readstat .ne. 0) then
      i=i-1
      exit
   end if   

end do
!
!      If no useful line is present in the file; throw an error
!
if (i.eq.0) then
   write(*,*) "The file qmdff_corr.dat contains no useful lines!"
   call fatal
end if

!
!     store total number of correction terms
!
corr_number=i
!
!     store correction informations in global arrays
!
allocate(corr_at1(corr_number),corr_at2(corr_number))
allocate(corr_factor(corr_number))
corr_at1=corr_at1_loc(1:corr_number)
corr_at2=corr_at2_loc(1:corr_number)
corr_factor=corr_factor_loc(1:corr_number)

close(97)

return
end subroutine qmdff_corr
