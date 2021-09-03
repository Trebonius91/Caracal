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
!     subroutine heat: calculate thermochemical properties for 
!           dissociation etc.
!
!     following code produces correct zero HOF for the
!     gas diatomics if the correct D_e is supplied as input
!     Mon May 26 12:18:39 CEST 2014
!
!     part of QMDFF
!
subroutine heat(n,at,de_mol,zpve_h298m,freq)
use qmdff
implicit none

real(kind=8)::freq(3*n) 
real(kind=8)::de_mol
integer::n,at(n)
integer::nelem(94)
real(kind=8)::so(94),h298a(94)
real(kind=8)::au2cal, so_sum, zpve_h298m, corr
real(kind=8)::dum1, hof, pv               
integer::i
logical::pr
!
!     EXPERIMENTAL ATOMIC DATA (in kcal/mol) from
!     Assessment of Gaussian-2 and density functional theories
!     for the computation of enthalpies of formation
!     Larry A. Curtiss
!     Krishnan Raghavachari
!     Paul C. Redfern
!     John A. Pople
!     J. Chem. Phys. 106 (3), 1997, 1063
!
pr=.false.
so   =0
h298a=0

au2cal=627.509541
pv    =0.001*298.15*8.31451/4.184
!
!     H298-H0 
!
h298a(1)=1.01
h298a(3)=1.10
h298a(4)=0.46
h298a(5)=0.29
h298a(6)=0.25
h298a(7)=1.04
h298a(8)=1.04
h298a(9)=1.05
h298a(11)=1.54
h298a(12)=1.19
h298a(13)=1.08
h298a(14)=0.76
h298a(15)=1.28
h298a(16)=1.05
h298a(17)=1.10
!
!      SO (not not used)
!
so(5) =-0.05
so(6) =-0.14
so(8) =-0.36
so(9) =-0.61
so(13)=-0.34
so(14)=-0.68
so(16)=-0.89
so(17)=-1.34
write(10,*)
write(10,*) "Thermochemical properties (HOF=heat of formation):"
write(10,*)
write(10,'('' ZPVE+H0-H298 (kcal) = '',f8.3)')zpve_h298m     

nelem=0
do i=1,n
   nelem(at(i))=nelem(at(i))+1
end do

so_sum=0
do i=1,94 
   so_sum=so_sum+nelem(i)*so(i)
end do

de_mol=au2cal*de_mol
if (pr) write(10,'(''De(incl.SO )'',f10.3)')de_mol+so_sum

dum1=0
!
!     atomic 298 K enthalpies 
!
do i=1,94 
   dum1=dum1+nelem(i)*(eheat(i)-pv-h298a(i))
end do
hof=dum1-de_mol+zpve_h298m
corr=dum1+zpve_h298m
!
!     this value should be compared to exp:
!
write(10,'('' HOF->De (kcal) =  '',f10.3)')corr
write(10,'('' HOF@298 (kcal) =  '',f10.3)')hof

!open(unit=137,file='.EAT')
!write(137,*) de_mol
!close(137)

return
end subroutine heat 

