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
!     subroutine prog_initial: sets up original values for some parameters and
!     variables that might not otherwise get initialized at beginning of programs
!
!     part of EVB
!
subroutine prog_initial(rank)
use general
use qmdff
use evb_mod
use fftw_mod
use omp_lib

implicit none
real*8 precise
integer::rank  ! the current MPI process
logical first
integer::id,threads

first=.true.
!
!     default unit numbers for input and output
!
input = 5
iout = 6
!
!     display program banner and copyright notice
!
if (rank .eq. 0) then
   call promo
end if
!
!     number of lines in the keyfile
!
nkey_lines = 0
!
!     number of atoms in the system
!
n = 0
!
!     flags for information levels within the program
!
debug = .false.
abort = .false.
!
!     flags for temperature and pressure baths
!
isothermal = .false.
isobaric = .false.
!
!     Set initial parameter for Fast Fourier Transform
!
Np=0
!
!     Give information about used number of OMP threads
!
!$OMP Parallel private(id) shared(threads)
threads = omp_get_num_threads()
!$OMP end Parallel
!if (rank .eq. 0) then
!   write(*,'(a,i5,a)') " This calculation will use ",threads," OpenMP threads on common memory."
!end if

return
end subroutine prog_initial
