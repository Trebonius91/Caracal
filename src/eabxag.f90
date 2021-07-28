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
!     subroutine eabxag: calculate energy/gradient of single halogen bond
!          (numerical gradient)
!
!     part of QMDFF
!
subroutine eabxag(n,A,B,H,xyz,ca,energy,gdr)
use debug
implicit none
integer::A,B,H,n
real(kind=8)::xyz(3,n),ca,energy,gdr(3,n)
real(kind=8)::cosabh,aterm,rdampl
real(kind=8)::rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
real(kind=8)::drah(3),drbh(3),drab(3)
real(kind=8)::ga(3),gb(3),gh(3),dg(3),dg2(3)
real(kind=8)::gi,denom,ratio,apref,aprod
real(kind=8)::eabh
REAL(kind=8)::longcut=120.d0
REAL(kind=8)::alp= 6.d0
REAL(kind=8)::alp3=6.d0
real(kind=8)::step,er,el,dum
integer::i,j,i1,i2

i1=A
i2=B
i =H

call eabx(n,i1,i2,i,xyz,ca,er)
energy=energy+er
!
!     If debugging is activated: store single energy parts
!
!if (do_debug) then
!   if (qmdff_num.eq.1) then
!      comp_no=comp_no+1
!      act_part=er
!      qmdff_parts(1,struc_no,comp_no)=act_part
!      write(parts_labels(1,comp_no),'(a,i4,i4,i4,i4)') "x-bond",A,B,H
!  else if (qmdff_num.eq.2) then
!      comp_no2=comp_no2+1
!      act_part=er
!      qmdff_parts(2,struc_no,comp_no2)=act_part
!      write(parts_labels(2,comp_no2),'(a,i4,i4,i4,i4)') "x-bond",A,B,H
!   end if
!end if

!
!     calculate numerical gradient: elogante coordinates!
!
step=1.d-6
dum=1./(2.0d0*step)
do j=1,3
   xyz(j,i1)=xyz(j,i1)+step
   call eabx(n,i1,i2,i,xyz,ca,er)
   xyz(j,i1)=xyz(j,i1)-step*2.0d0
   call eabx(n,i1,i2,i,xyz,ca,el)
   xyz(j,i1)=xyz(j,i1)+step
   gdr(j,i1)=gdr(j,i1)+(er-el)*dum
end do
do j=1,3
   xyz(j,i2)=xyz(j,i2)+step
   call eabx(n,i1,i2,i,xyz,ca,er)
   xyz(j,i2)=xyz(j,i2)-step*2.0d0
   call eabx(n,i1,i2,i,xyz,ca,el)
   xyz(j,i2)=xyz(j,i2)+step
   gdr(j,i2)=gdr(j,i2)+(er-el)*dum 
end do
do j=1,3
   xyz(j,i )=xyz(j,i )+step
   call eabx(n,i1,i2,i,xyz,ca,er)
   xyz(j,i )=xyz(j,i )-step*2.0d0
   call eabx(n,i1,i2,i,xyz,ca,el)
   xyz(j,i )=xyz(j,i )+step
   gdr(j,i )=gdr(j,i )+(er-el)*dum
end do

return
end subroutine eabxag
