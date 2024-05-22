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
!     ###############################################################
!     ##                                                           ##
!     ##  subroutine atommass_loc  --  initialize atomic masses    ##
!     ##                                                           ##
!     ###############################################################
!
!     Normally the atomic masses needed for dynamics are read in 
!     from the ffield file.
!     But in the case of EVB-QMDFF no masses are specified in the ffield
!     so they need to defined here in this file.
!     In this routine, the masses are given back as single reals!
!     The masses from H up to Rn are availiable 
!     They are taken from "Das große Tafelwerk, 1. Auflage 2003"
!
!     part of EVB
!
subroutine atommass_loc(adum,mass_act)
use general
use evb_mod
implicit none
integer::i
character(len=2)::adum
real(kind=8)::mass_act
!
!     initialize the element-masses
!     Convert all letters to upcase letters in order to avoid 
!     errors with atom symbols!
!
call upcase(adum)
if (adum .eq. "H") then
   mass_act=1.00782503207d0
else if (adum .eq. "D") then
   mass_act=2.0141017778d0
else if (adum .eq. "HE") then
   mass_act=4.00260d0
else if (adum .eq. "LI") then
   mass_act=6.94000d0
else if (adum .eq. "BE") then
   mass_act=9.01218d0
else if (adum .eq. "B") then
   mass_act=10.81000d0
else if (adum .eq. "C") then
   mass_act=12.00000d0
else if (adum .eq. "N") then
   mass_act=14.0030740048d0
else if (adum .eq. "O") then
   mass_act=15.99491461956d0
else if (adum .eq. "F") then
   mass_act=18.99840d0
else if (adum .eq. "NE") then
   mass_act=20.18d0
else if (adum .eq. "NA") then
   mass_act=23.0d0
else if (adum .eq. "MG") then
   mass_act=24.31d0
else if (adum .eq. "AL") then
   mass_act=26.98d0
else if (adum .eq. "SI") then
   mass_act=28.09d0
else if (adum .eq. "P") then
   mass_act=30.97d0
else if (adum .eq. "S") then
   mass_act=32.06000d0
else if (adum .eq. "CL") then
   mass_act=34.96885268d0
else if (adum .eq. "AR") then
   mass_act=39.95d0
else if (adum .eq. "K") then
   mass_act=39.1d0
else if (adum .eq. "CA") then
   mass_act=40.08d0
else if (adum .eq. "SC") then
   mass_act=44.96d0
else if (adum .eq. "TI") then
   mass_act=47.88d0
else if (adum .eq. "V") then
   mass_act=50.94d0
else if (adum .eq. "CR") then
   mass_act=52.0d0
else if (adum .eq. "MN") then
   mass_act=54.94d0
else if (adum .eq. "FE") then
   mass_act=55.85d0
else if (adum .eq. "CO") then
   mass_act=58.93d0
else if (adum .eq. "NI") then
   mass_act=58.69d0
else if (adum .eq. "CU") then
   mass_act=63.55d0
else if (adum .eq. "ZN") then
   mass_act=65.39d0
else if (adum .eq. "GA") then
   mass_act=69.72d0
else if (adum .eq. "GE") then
   mass_act=72.61d0
else if (adum .eq. "AS") then
   mass_act=74.92d0
else if (adum .eq. "SE") then
   mass_act=78.96d0
else if (adum .eq. "BR") then
   mass_act=78.9183371d0
else if (adum .eq. "KR") then
   mass_act=83.80d0
else if (adum .eq. "RB") then
   mass_act=85.47d0
else if (adum .eq. "SR") then
   mass_act=87.62d0
else if (adum .eq. "Y") then
   mass_act=88.91d0
else if (adum .eq. "ZR") then
   mass_act=91.22d0
else if (adum .eq. "NB") then
   mass_act=92.91d0
else if (adum .eq. "MO") then
   mass_act=95.94d0
else if (adum .eq. "TC") then
   mass_act=98d0
else if (adum .eq. "RU") then
   mass_act=101.07d0
else if (adum .eq. "RH") then
   mass_act=102.91d0
else if (adum .eq. "PD") then
   mass_act=106.42d0
else if (adum .eq. "AG") then
   mass_act=107.87d0
else if (adum .eq. "CD") then
   mass_act=112.41d0
else if (adum .eq. "IN") then
   mass_act=114.82d0
else if (adum .eq. "SN") then
   mass_act=118.71d0
else if (adum .eq. "SB") then
   mass_act=121.76d0
else if (adum .eq. "TE") then
   mass_act=127.60d0
else if (adum .eq. "I") then
   mass_act=125.90d0
else if (adum .eq. "XE") then
   mass_act=131.29d0
else if (adum .eq. "CS") then
   mass_act=132.91d0
else if (adum .eq. "BA") then
   mass_act=137.33d0
else if (adum .eq. "LA") then
   mass_act=138.91d0
else if (adum .eq. "HF") then
   mass_act=178.49d0
else if (adum .eq. "TA") then
   mass_act=180.95d0
else if (adum .eq. "W") then
   mass_act=183.95d0
else if (adum .eq. "RE") then
   mass_act=186.21d0
else if (adum .eq. "OS") then
   mass_act=190.23d0
else if (adum .eq. "IR") then
   mass_act=192.22d0
else if (adum .eq. "PT") then
   mass_act=195.08d0
else if (adum .eq. "AU") then
   mass_act=196.97d0
else if (adum .eq. "HG") then
   mass_act=200.59d0
else if (adum .eq. "TL") then
   mass_act=204.38d0
else if (adum .eq. "PB") then
   mass_act=207.2d0
else if (adum .eq. "BI") then
   mass_act=208.98d0
else if (adum .eq. "PO") then
   mass_act=209d0
else if (adum .eq. "AT") then
   mass_act=210d0
else if (adum .eq. "RN") then
   mass_act=222d0
end if
!
!     If masses shall be changed manually
!
if (change_mass) then
   if (adum .eq. elem_mass) then
      mass_act=newmass
   end if
end if


if (use_rpmd) then
   mass_act=mass_act/emass
!  test for RPMDrate
!   mass_act=mass_act*0.001d0/6.02214179d+23/9.1093826e-31
end if
 
return
end subroutine atommass_loc
