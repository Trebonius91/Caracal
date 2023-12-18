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
!     subroutine eabhag: calculate energy/gradiet of single H-bond
!     --> analytical gradients!
!
!     part of QMDFF
!
subroutine eabhag(n,A,B,H,xyz,ca,cb,energy,gdr)
use debug
use qmdff
use pbc_mod
implicit none
integer::A,B,H,n,i
real(kind=8)::xyz(3,n),ca,cb,energy,gdr(3,n)
real(kind=8)::cosabh,aterm,rdampl,da
real(kind=8)::rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
real(kind=8)::drah(3),drbh(3),drab(3)
real(kind=8)::ga(3),gb(3),gh(3),dg(3),dg2(3)
real(kind=8)::gi,denom,ratio,apref,aprod
real(kind=8)::eabh
REAL(kind=8)::longcut=8.d0
REAL(kind=8)::alp7= 12.d0
REAL(kind=8)::alp3= 6.d0
!
!     compute distances
!
drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
drab(1:3)=xyz(1:3,A)-xyz(1:3,B)
!
!     apply periodic boundaries, if needed
!
if (periodic) then
   call box_image(drah)
end if

!
!     A-B distance
!
rab2=sum(drab**2)
rab =sqrt(rab2)
!
!     A-H distance
!
rah2=sum(drah**2)
rah =sqrt(rah2)
!
!     B-H distance
!
rbh2=sum(drbh**2)
rbh =sqrt(rbh2)
!
!     long  damping
!
ratio = (rab/longcut)**alp7
rdampl=1.d0/(1.d0+ratio)
rdampl=rdampl/rab2/rab
!
!     cos angle A-B and H (or B-H H) is term
!
if (rah2.gt.rbh2) then
  aprod = 1.d0/rbh/rab
  cosabh = -sum(drbh*drab)*aprod
else
  aprod = 1.d0/rah/rab
  cosabh = sum(drah*drab)*aprod
end if
!
!     angle damping term
!
aterm = 0.5d0*(cosabh+1.d0)
apref = aterm**(alp3-1)
aterm = aterm*apref
apref = alp3*0.5d0*apref
!
!     donor-acceptor term
!
rah4 = rah2*rah2
rbh4 = rbh2*rbh2
denom = 1.d0/(rah4+rbh4)
da = (ca*rah4 + cb*rbh4)*denom
!
!     r^3 only slightly better than r^4 (SG, 8/12)
!
eabh = -da*rdampl*aterm

if (eabh.gt.-1.d-8) return

energy = energy + eabh
!
!     For debugging: print out the single components
!
!if (do_debug) then
!   if (qmdff_num.eq.1) then
!      comp_no=comp_no+1
!      act_part=eabh
!      qmdff_parts(1,struc_no,comp_no)=act_part
!      write(parts_labels(1,comp_no),'(a,i4,i4,i4,i4)') "h-bond",A,B,H
!   else if (qmdff_num.eq.2) then
!      comp_no2=comp_no2+1
!      act_part=eabh
!      qmdff_parts(2,struc_no,comp_no2)=act_part
!      write(parts_labels(2,comp_no2),'(a,i4,i4,i4,i4)') "h-bond",A,B,H
!   end if
!end if
!
!     gradient  
!     donor-acceptor part  
!  
gi = 4.d0*(ca-cb)*rah2*rbh4*denom*denom
gi = -gi*rdampl*aterm
ga(1:3) = gi*drah(1:3)
gi = 4.d0*(cb-ca)*rbh2*rah4*denom*denom
gi = -gi*rdampl*aterm
gb(1:3) = gi*drbh(1:3)
gh(1:3) = -ga(1:3)-gb(1:3)
!
!     long-range damping part
!
gi = rdampl*rdampl*rab*(3.d0+(3.d0+alp7)*ratio)
gi = gi*da*aterm
dg(1:3) = gi*drab(1:3)
ga(1:3) = ga(1:3) + dg(1:3)
gb(1:3) = gb(1:3) - dg(1:3)
!
!     angle part
!
if (rah2.gt.rbh2) then
  dg(1:3) = -aprod*drbh(1:3) - cosabh*drab(1:3)/rab2
  dg2(1:3) = aprod*drab(1:3) + cosabh*drbh(1:3)/rbh2
  gi = -da*rdampl*apref
  dg(1:3) = gi*dg(1:3)
  dg2(1:3) = gi*dg2(1:3)

  ga(1:3) = ga(1:3) + dg(1:3)
  gh(1:3) = gh(1:3) + dg2(1:3)
  gb(1:3) = gb(1:3) - dg(1:3) - dg2(1:3)

else
  dg(1:3) = -aprod*drah(1:3) + cosabh*drab(1:3)/rab2
  dg2(1:3) = -aprod*drab(1:3) + cosabh*drah(1:3)/rah2
  gi = -da*rdampl*apref
  dg(1:3) = gi*dg(1:3)
  dg2(1:3) = gi*dg2(1:3)

  gb(1:3) = gb(1:3) + dg(1:3)
  gh(1:3) = gh(1:3) + dg2(1:3)
  ga(1:3) = ga(1:3) - dg(1:3) - dg2(1:3)

end if
!
!     move gradients into place
!
gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
gdr(1:3,H) = gdr(1:3,H) + gh(1:3)

!
!     Calculate the virial tensor components, if needed!   
!

if (calc_vir) then
   vir_ten(1,1)=vir_ten(1,1)+drah(1)*ga(1)+drbh(1)*gb(1)
   vir_ten(2,1)=vir_ten(2,1)+drah(2)*ga(1)+drbh(2)*gb(1)
   vir_ten(3,1)=vir_ten(3,1)+drah(3)*ga(1)+drbh(3)*gb(1)
   vir_ten(1,2)=vir_ten(1,2)+drah(1)*ga(2)+drbh(1)*gb(2)
   vir_ten(2,2)=vir_ten(2,2)+drah(2)*ga(2)+drbh(2)*gb(2)
   vir_ten(3,2)=vir_ten(3,2)+drah(3)*ga(2)+drbh(3)*gb(2)
   vir_ten(1,3)=vir_ten(1,3)+drah(1)*ga(3)+drbh(1)*gb(3)
   vir_ten(2,3)=vir_ten(2,3)+drah(2)*ga(3)+drbh(2)*gb(3)
   vir_ten(3,3)=vir_ten(3,3)+drah(3)*ga(3)+drbh(3)*gb(3)
end if

return
end subroutine eabhag
