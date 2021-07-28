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
!     suboutine g98fake: write a Gaussian98 fake output such that 
!             GaussView can display the optimized QMDFF vibrations
!
!     part of QMDFF
!
subroutine g98fake(fname,n,at,xyz,freq,u2,u)
implicit none
integer::n,at(n)
real(kind=8)::u(3*n,3*n),freq(3*n),xyz(3,n),u2(3*n,3*n)
character(len=*) fname
integer::gu,i,j,ka,kb,kc,la,lb,k
character(len=2)::irrep
real(kind=8)::red_mass(3*n)
real(kind=8)::force(3*n)
real(kind=8)::ir_int(3*n)
real(kind=8)::f2(3*n)
real(kind=8)::zero

irrep='a'
red_mass=99.0
force   =99.0
ir_int  =99.0
zero    =0.0

k=0
do i=1,3*n
   if (abs(freq(i)).gt.1.d-2) then
      k=k+1
      u(1:3*n,k)=u2(1:3*n,i)
      f2(k)=freq(i)
   end if
end do

write (10,'(" writing <",a,"> molden file with frequencies/normal modes")') trim(fname)
write(10,*) "Open it with molden to see spectrum and normal modes of the QMDFF!"
gu=55

open(unit=gu,file=fname)
write (gu,'('' Entering Gaussian System'')')
write (gu,'('' *********************************************'')')
write (gu,'('' Gaussian 98:'')')
write (gu,'('' frequency fake output'')')
write (gu,'('' *********************************************'')')

write (gu,*) '                        Standard orientation:'
write (gu,*) '---------------------------------------------', &
     &           '-----------------------'
write (gu,*) ' Center     Atomic     Atomic', &
     &         '              Coordinates (Angstroms)'
write (gu,*) ' Number     Number      Type ', &
     &         '             X           Y           Z'
write (gu,*) '-----------------------', &
     &         '---------------------------------------------'
j=0
do i=1,n
   write(gu,111) i,at(i),j,xyz(1:3,i)*0.52917726
end do
write (gu,*) '----------------------', &
     &           '----------------------------------------------'
write (gu,*) '    1 basis functions        1 primitive gaussians'
write (gu,*) '    1 alpha electrons        1 beta electrons'
write (gu,*)

111   format(i5,i11,i14,4x,3f12.6)

write (gu,*) 'Harmonic frequencies (cm**-1), IR intensities', &
     &         ' (KM/Mole),'
write (gu,*) 'Raman scattering activities (A**4/amu),', &
     &         ' Raman depolarization ratios,'
write (gu,*) 'reduced masses (AMU), force constants (mDyne/A)', &
     &         ' and normal coordinates:'

ka=1
kc=3
60 continue
kb=min0(kc,k)  
write (gu,100) (j,j=ka,kb)
write (gu,105) (irrep,j=ka,kb)
write (gu,110) ' Frequencies --',(f2(j),j=ka,kb)
write (gu,110) ' Red. masses --',(red_mass(j),j=ka,kb)
write (gu,110) ' Frc consts  --',(force(j),j=ka,kb)
write (gu,110) ' IR Inten    --',(ir_int(j),j=ka,kb)
write (gu,110) ' Raman Activ --',(zero,j=ka,kb)
write (gu,110) ' Depolar     --',(zero,j=ka,kb)
write (gu,*)'Atom AN      X      Y      Z        X      Y', &
     &      '      Z        X      Y      Z'
la=1
70 continue
lb=n
do i=la,lb
   write (gu,130) i,at(i), &
      &  (u(i*3-2,j), &
      &   u(i*3-1,j), &
      &   u(i*3  ,j),j=ka,kb)
end do
if (lb.eq.n) goto 90
go to 70
90 continue

if (kb.eq.k) then
   return
end if

ka=kc+1
kc=kc+3
goto 60

100  format (3(20x,i3))
105  format (3x,3(18x,a5))
110  format (a15,f11.4,12x,f11.4,12x,f11.4)
130  format (2i4,3(2x,3f7.2))

close(gu)

return
end subroutine 
