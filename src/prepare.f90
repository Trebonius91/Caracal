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
!     subroutine prepare: reads in the input of already stated QMDFF´s and coupling
!     informations and stores them in the program for being used later
!
!     part of EVB
!     
subroutine prepare (fffilen1,fffilen2,fffilen3,qmdffnumber)
use general 
use qmdff
use evb_mod
use pbc_mod

implicit none

character(len=*)::fffilen1,fffilen2,fffilen3
integer::iz1_two,iz2_two,ndummy
integer::iz1_three,iz2_three
integer::parameters,qmdffnumber
real(kind=8)::dens_two,c6_two,dens_three,c6_three
real(kind=8)::vz(94)
real(kind=8)::scalehb(94),scalexb(94)
real(kind=8)::dens
real(kind=8)::r42,c6
integer::i,j,i1,i2,iz1,iz2,lin,n2
character(len=80)::fname
parameters=3
!
!     Set the periodic switch to false (unless activated by dynamic.x)
!
!periodic=.false.
!
!     read number of atoms, allocate and check that all qmdffs have the 
!     same number of atoms
!
call rdsolvff0(n_one,fffilen1)
if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
   call rdsolvff0(n_two,fffilen2)      
   if (n_one.ne.n_two) then
      write(*,*)'solvent1 has ',n_one,' atoms'
      write(*,*)'solvent2 has ',n_two,' atoms'
      write(*,*)'   This will not work!'
      stop 'atom number mismatch between FF1 and FF2'
   end if
end if
if (qmdffnumber.eq.3) then
   call rdsolvff0(n_three,fffilen3)
   if (n_one.ne.n_three) then
      write(*,*)'solvent1 has ',n_one,' atoms'
      write(*,*)'solvent3 has ',n_three,' atoms'
      write(*,*)'   This will not work!'
      stop 'atom number mismatch between FF1 and FF3'
   end if
end if
natoms=n_one
!
!     now both atomnumbers are assumed to be identical
!

allocate(at(n_one),xyz(3,n_one),at2(n_one), &
   &     g_one(3,n_one),q(n_one),c6xy(n_one,n_one),cn(n_one), &
   &     imass(n_one))
if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
   allocate(at_two(n_one),xyz_two(3,n_one),g_two(3,n_one), &
      &     q_two(n_one),c6xy_two(n_one,n_one),cn_two(n_one), &
      &     imass_two(n_one))
end if
if (qmdffnumber.eq.3) then
   allocate(at_three(n_one),xyz_three(3,n_one),g_three(3,n_one), &
     &     q_three(n_one),c6xy_three(n_one,n_one),cn_three(n_one), &
     &     imass_three(n_one))
end if
call copyc6 ! not molecule/FF specific
call setnonb(scalehb,scalexb,vz,sr42,zab,r0ab)
r094_mod=r094
!
!     read FF file 1
!
call rdsolvff(n_one,xyz,at,q,imass,dens,scalehb,scalexb,fffilen1) 
!
!     read FF file 2
!
if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
   call rdsolvff_two(n_one,xyz_two,at_two,q_two,imass_two,dens_two, &
           &     scalehb,scalexb,fffilen2)
end if
if (qmdffnumber.eq.3) then
   call rdsolvff_three(n_one,xyz_three,at_three,q_three,imass_three, &
           &     dens_three,scalehb,scalexb,fffilen3)
end if
!
!     compute CN and pair C6
!
call ncoord_qmdff(n_one,rcov,at,xyz,cn,5000.0d0)
do i1=1,n_one
   iz1=at(i1)
   do i2=1,i1
      iz2=at(i2)
      call getc6(5,94,c6ab,maxci,iz1,iz2,cn(i1),cn(i2),c6)
      c6xy(i2,i1)=c6 
      c6xy(i1,i2)=c6 
   end do
end do 
if (qmdffnumber.eq.2 .or. qmdffnumber.eq.3) then
   call ncoord_qmdff(n_one,rcov,at_two,xyz_two,cn_two,5000.0d0)
   do i1=1,n_one
      iz1_two=at_two(i1)
      do i2=1,i1
         iz2_two=at_two(i2)
         call getc6(5,94,c6ab,maxci,iz1_two,iz2_two, &
                &     cn_two(i1),cn_two(i2),c6_two)
         c6xy_two(i1,i2)=c6_two
         c6xy_two(i2,i1)=c6_two
      end do
   end do
end if 
if (qmdffnumber.eq.3) then
   call ncoord_qmdff(n_one,rcov,at_three,xyz_three,cn_three,5000.0d0)
   do i1=1,n_one
      iz1_three=at_three(i1)
      do i2=1,i1
         iz2_three=at_three(i2)
         call getc6(5,94,c6ab,maxci,iz1_three,iz2_three, &
                &     cn_three(i1),cn_three(i2),c6_three)
         c6xy_three(i1,i2)=c6_three
         c6xy_three(i2,i1)=c6_three
      end do
   end do
end if 
!
!     Default: no correction will be applied to QMDFFs
!        "QMDFF zero will be corrected"
! 
corr_nonb=0
!
!     If the noncovalent QMDFF parameters were optimized separately
!
ff_mod_noncov=.false.
if (fffilen2 .eq. "dummy" .and. fffilen3 .eq. "dummy") then
   ff_mod_noncov=.true.
end if

!
!     COVALENT RADII of all elements
!     based on "Atomic Radii of the Elements," M. Mantina, R. Valero, C. J. Cramer, and D. G. Truhlar,
!     in CRC Handbook of Chemistry and Physics, 91st Edition (2010-2011),
!     edited by W. M. Haynes (CRC Press, Boca Raton, FL, 2010), pages 9-49-9-50;
!     corrected Nov. 17, 2010 for the 92nd edition.
!
rad = (/ &
 & 0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0, &
 & 0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0, &
 & 2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0, &
 & 1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0, &
 & 2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0, &
 & 1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0, &
 & 2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0, &
 & 1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0, &
 & 1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0, &
 & 1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0, &
 & 1.84D0,1.83D0,1.80D0,1.80D0 /)
!
!     Convert radii to bohr
!
rad=rad/bohr


return
end subroutine prepare
