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
!     subroutine setnonb: Set various arrays for QMDFF non-bonding potential, not 
!                molecule specific
!
!     part of QMDFF
!
subroutine setnonb(scalehb,scalexb,vz,sr42,zab,r0ab)
use qmdff
implicit none
real(kind=8)::vz(94),zab(94,94),r0ab(94,94)
real(kind=8)::sr42(94,94),scalehb(94),scalexb(94)

integer::i,j
!
!     D3 (Becke-Johnson damping) parameters
!     eq. (11) in QMDFF paper
!
a1 = 0.45
s8 = 2.7 
a2 = 4.0
!
!     unknown parameters for damping function (?)
!

r2r4 = (/ &
 & 2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
 & 3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
 & 6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
 & 4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
 & 6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
 & 5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
 & 6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
 & 4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
 & 6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
 & 5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
 & 5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
 &10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
 & 8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
 & 8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
 & 7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
 & 6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
 & 7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
 & 6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
 & 8.77140725,  8.65402716,  8.53923501,  8.85024712 /)
!
!     R0AB parameters in damping function 
!
rcov = (/ &
 & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, & 
 & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
 & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
 & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
 & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
 & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
 & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
 & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, & 
 & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
 & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
 & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
 & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
 & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
 & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
 & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
 & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
 & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
 & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
 & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

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
!     as in qmsolv (D3 is in block.f):
!     scaling parameters for hydrogen bonds (HB)
!     only for H,O,N,... unequal zero
!
scalehb=0
scalehb(7  )=0.8
scalehb(8  )=0.3
scalehb(9  )=0.1
scalehb(15 )=2.0
scalehb(16 )=2.0
scalehb(17 )=2.0
scalehb(34 )=2.0
scalehb(35 )=2.0
!
!     scaling parameters for halogen bonds (XB)
!
scalexb=0
scalexb(17 )=0.30
scalexb(35 )=0.60
scalexb(53 )=0.80
scalexb(85 )=1.00
!
!     valence charges (REP), 4 fit parameters included there
!
call valel(vz)
!
!     determine D3 parameters 
!
do i=1,94
   do j=1,94
      sr42(j,i)=3.0d0*s8*r2r4(i)*r2r4(j)
      r094(j,i)=a1*sqrt(3.0d0*r2r4(i)*r2r4(j))+a2
      zab (j,i)=vz(i)*vz(j)
   end do
end do
!
!     get D3 radii (used in REP), alpha included for calc. of REP energy 
!     Pauli repulsion formula (eq. (12))
!
call setr0(94,r0ab)
r0ab = 16.5 / r0ab**1.5
!
!     intramolecular damping parameters (Becke-Johnson)
!     fitted to aconf,pconf,cyconf,sconf,mconf      
!
eps1(1)=0
eps1(2)=0
eps1(3)=0.85
eps1(4)=1
eps1(5)=1

!
!     rot34 geometries very sensitive to scaling     
! 
eps2(1)=0
eps2(2)=0
eps2(3)=0.5
eps2(4)=0.5
eps2(5)=1

!
!     repulsion-dispersion for 1,3-metal. 
!
eps2(6)=1.00
eps1(6)=0
!
!     overlap scaling    
!  
eps3(1:3)=1
eps3(4)=0.5
eps3(5)=0
eps3(6)=0

end subroutine setnonb
