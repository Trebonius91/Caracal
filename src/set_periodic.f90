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
!     subroutine set_periodic: set module parameters for periodic
!          calculations in the qmdff module
!

subroutine set_periodic(switch_per,length,switch_ew,brute_force, &
   & cut_coul,cut_vdw)
use qmdff
implicit none
integer::i,k
logical::switch_per,switch_ew,brute_force
integer::ifft
real(kind=8)::length
real(kind=8)::ewald_accuracy
real(kind=8),allocatable::array(:),bsarray(:)
real(kind=8)::ratio,xhi,xlo,y_par,x_par,eps
real(kind=8)::dens,delta
real(kind=8)::cut_coul,cut_vdw
integer::minfft,maxfft
integer::maxpower
parameter (maxpower=63)
integer::multi(maxpower)
integer::ifront,iback,error,iguess  !for FFTW routines
integer::nthread
!
!    The possible PME grid sizes
!
data multi  /   2,   4,   6,   8,  10,  12,  16,  18,  20, &
         &    24,  30,  32,  36,  40,  48,  50,  54,  60,  & 
         &    64,  72,  80,  90,  96, 100, 108, 120, 128,  &
         &    144, 150, 160, 162, 180, 192, 200, 216, 240, &
         &    250, 256, 270, 288, 300, 320, 324, 360, 384, & 
         &    400, 432, 450, 480, 486, 500, 512, 540, 576, &
         &    600, 640, 648, 720, 750, 768, 800, 810, 864 /


periodic=switch_per
ewald=switch_ew
coul_cut=cut_coul
vdw_cut=cut_vdw/0.52917721092d0
ewald_brute=brute_force
box_len=length
box_len2=length*0.5d0!
!
!     set the Coulomb cutoff to half the box length if its value 
!     is too small (set it directly to bohr)
!
if (coul_cut .lt. 5d0) then
   coul_cut=box_len2-0.1d0
end if
!
!     If the Coulomb cutoff is larger than half the box length,
!     set it to half the box length
!
if (coul_cut .gt. box_len2) then
   coul_cut=box_len2-0.1d0
   write(*,*) "Warning: The Coulomb cutoff was chosen too large and will"
   write(*,*) " be set to half the periodic box length."
end if
!
!     Set the value where the smooth switch off of the Coulomb energy 
!     near the Coulomb cutoff begins
!
cut_low=coul_cut-4.d0
!
!     set the global Ewald parameters
!  
if (ewald) then

   r_ew_cut=7/0.52917721092d0
   sqrt_p=sqrt(-log(ewald_accuracy))

!
!     default grid size for reciprocal part from system dimensions
!
   minfft=16
   maxfft=864
   dens=1.2d0
   delta=1D-8
   ifft=int(box_len*0.52917721092d0*dens-delta)+1
   write(*,*) "Ifft",ifft,ewald

   nfft = maxfft
   do i = maxpower, 1, -1
      k = multi(i)
      if (k .le. maxfft) then
         if (k .ge. ifft)  nfft = k
      end if
   end do
   minfft = 16
   if (nfft .lt. minfft)  nfft = minfft

   write(*,*) "PME grid points per dimension: ",nfft

!
!     Set the fft-parameters for the Smooth Particle Mesh Ewald Method
!
   maxtable=4*nfft
   maxprime=15
   allocate(ffttable(maxtable,3))
   allocate(iprime(maxprime,3))

!
!     Initialize the FFT routines
!
   ifront = -1
   iback = 1
   error = 0
   iguess = 0
   nthread=1
   allocate(qgrid(2,nfft,nfft,nfft))
!call dfftw_init_threads (error)
!call dfftw_plan_with_nthreads (nthread)
   call dfftw_plan_dft_3d (planf,nfft,nfft,nfft,qgrid, &
      &                              qgrid,ifront,iguess)
   call dfftw_plan_dft_3d (planb,nfft,nfft,nfft,qgrid, &
      &                              qgrid,iback,iguess)
!
!     DETERMINATION OF ALPHA EWALD PARAMETER
!
!     set the tolerance value; use of 1.0d-8 results in large
!     Ewald coefficients that ensure continuity in the gradient
!
   eps = 1.0d-8
!
!     get approximate value from cutoff and tolerance
!
   ratio = eps + 1.0d0
   x_par = 0.5d0
   i = 0
   do while (ratio .ge. eps)
      i = i + 1
      x_par = 2.0d0 * x_par
      y_par = x_par * r_ew_cut
      ratio = erfc(y_par) / r_ew_cut
   end do
!
!     use a binary search to refine the coefficient
!
   k = i + 60
   xlo = 0.0d0
   xhi = x_par
   do i = 1, k
      x_par = (xlo+xhi) / 2.0d0
      y_par = x_par * r_ew_cut
      ratio = erfc(y_par) / r_ew_cut
      if (ratio .ge. eps) then
         xlo = x_par
      else
         xhi = x_par
      end if
   end do
   a_ewald = x_par

   write(*,*) "alpha",a_ewald
!
!     Set further parameters for the FFT in SPME
!
   bsorder=5
   allocate (array(bsorder))
   allocate (bsarray(nfft))
!
!     compute and load the moduli values
!
   call bspline (0.0d0,bsorder,array)
   do i = 1, nfft
      bsarray(i) = 0.0d0
   end do
   do i = 1, bsorder
      bsarray(i+1) = array(i)
   end do
   allocate(bsmod1(nfft),bsmod2(nfft),bsmod3(nfft))
   call dftmod (bsmod1,bsarray,nfft,bsorder)
   call dftmod (bsmod2,bsarray,nfft,bsorder)
   call dftmod (bsmod3,bsarray,nfft,bsorder)
!
!     perform deallocation of some local arrays
!
   deallocate (array)
   deallocate (bsarray)

end if


return

end subroutine set_periodic
