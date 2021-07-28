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
!     subroutine mdstat: called at each molecular dynamics time step to
!     form statistics on various average values and fluctuations,
!     and to periodically save the state of the trajectory
!
!     part of EVB
!
subroutine mdstat (istep,dt,etot,epot,ekin,temp,pres)
use general
use evb_mod
implicit none
integer::istep,modstep
real(kind=8)::dt,temp,pres
real(kind=8)::etot,epot,ekin
real(kind=8)::pico,dens
real(kind=8)::fluctuate,fluctuate2
real(kind=8)::intfluct,intfluct2
real(kind=8)::potfluct,potfluct2
real(kind=8)::kinfluct,kinfluct2
real(kind=8)::tfluct,pfluct,dfluct
real(kind=8)::tfluct2,pfluct2,dfluct2
real(kind=8)::etot_sum,etot2_sum
real(kind=8)::eint_sum,eint2_sum
real(kind=8)::etot_ave,etot2_ave
real(kind=8)::eint_ave,eint2_ave
real(kind=8)::epot_sum,epot2_sum
real(kind=8)::ekin_sum,ekin2_sum
real(kind=8)::epot_ave,epot2_ave
real(kind=8)::ekin_ave,ekin2_ave
real(kind=8)::temp_sum,temp2_sum
real(kind=8)::temp_ave,temp2_ave
real(kind=8)::pres_sum,pres2_sum
real(kind=8)::pres_ave,pres2_ave
real(kind=8)::dens_sum,dens2_sum
real(kind=8)::dens_ave,dens2_ave
save etot_sum,etot2_sum
save eint_sum,eint2_sum
save epot_sum,epot2_sum
save ekin_sum,ekin2_sum
save temp_sum,temp2_sum
save pres_sum,pres2_sum
save dens_sum,dens2_sum

!
!     set number of steps for block averages of properties
!
modstep = mod(istep,iprint)
!
!     zero out summation variables for new averaging period
!
if (modstep.eq.1 .or. iprint.eq.1) then
   etot_sum = 0.0d0
   etot2_sum = 0.0d0
   epot_sum = 0.0d0
   epot2_sum = 0.0d0
   ekin_sum = 0.0d0
   ekin2_sum = 0.0d0
   eint_sum = 0.0d0
   eint2_sum = 0.0d0
   temp_sum = 0.0d0
   temp2_sum = 0.0d0
   pres_sum = 0.0d0
   pres2_sum = 0.0d0
   dens_sum = 0.0d0
   dens2_sum = 0.0d0
end if

!
!     print header for the averages over a group of recent steps
!
if (modstep .eq. 0) then
   pico = dble(istep) * dt
   if (writestat) then
      write (iout,'(/," Average Values for the Last",i6," Out of", &
          &      i9," Dynamics Steps")')  iprint,istep
      write (iout,'(/," Simulation Time",5x,f15.4," Picosecond")')  pico
   end if
end if
!
!     compute total energy and fluctuation for recent steps
!
etot_sum = etot_sum + etot
etot2_sum = etot2_sum + etot**2
if (modstep .eq. 0) then
   etot_ave = etot_sum / dble(iprint)
   etot2_ave = etot2_sum / dble(iprint)
   fluctuate2 = etot2_ave - etot_ave**2
   if (fluctuate2 .gt. 0.0d0) then
      fluctuate = sqrt(fluctuate2)
   else
      fluctuate = 0.0d0
   end if
   if (writestat) then
      write (iout,'(" Total Energy",8x,f15.4," Kcal/mole",3x, &
          &       "(+/-",f9.4,")")')  etot_ave,fluctuate
   end if
end if
!
!     compute average potential energy and its fluctuation
!
epot_sum = epot_sum + epot
epot2_sum = epot2_sum + epot**2
if (modstep .eq. 0) then
   epot_ave = epot_sum / dble(iprint)
   epot2_ave = epot2_sum / dble(iprint)
   potfluct2 = epot2_ave - epot_ave**2
   if (potfluct2 .gt. 0.0d0) then
      potfluct = sqrt(potfluct2)
   else
      potfluct = 0.0d0
   end if
   if (writestat) then
      write (iout,'(" Potential Energy",4x,f15.4," Kcal/mole",3x, &
         &    "(+/-",f9.4,")")')  epot_ave,potfluct
   end if
end if
!
!     compute average kinetic energy and its fluctuation
!
ekin_sum = ekin_sum + ekin
ekin2_sum = ekin2_sum + ekin**2
if (modstep .eq. 0) then
   ekin_ave = ekin_sum / dble(iprint)
   ekin2_ave = ekin2_sum / dble(iprint)
   kinfluct2 = ekin2_ave - ekin_ave**2
   if (kinfluct2 .gt. 0.0d0) then
      kinfluct = sqrt(kinfluct2)
   else
      kinfluct = 0.0d0
   end if
   if (writestat) then
      write (iout,'(" Kinetic Energy",6x,f15.4," Kcal/mole",3x, &
         &       "(+/-",f9.2,")")')  ekin_ave,kinfluct
   end if
end if
!
!     compute the average temperature and its fluctuation
!
temp_sum = temp_sum + temp
temp2_sum = temp2_sum + temp**2
if (modstep .eq. 0) then
   temp_ave = temp_sum / dble(iprint)
   temp2_ave = temp2_sum / dble(iprint)
   tfluct2 = temp2_ave - temp_ave**2
   if (tfluct2 .gt. 0.0d0) then
      tfluct = sqrt(tfluct2)
   else
      tfluct = 0.0d0
   end if
   if (writestat) then
      write (iout,'(" Temperature",9x,f15.2," Kelvin",6x, &
         &       "(+/-",f9.2,")")')  temp_ave,tfluct
   end if
end if

! TEST: print average temperatures!
!if (modstep .eq. 0) then
!write(99,*) temp_ave
!if (istep .gt. 10000) then
!   t_avg=t_avg+temp_ave
!end if
!end if
return
end subroutine mdstat 
