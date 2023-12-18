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
!     function dihed: calculate dihedral angle between at1-at2-at3-at4
!
!     part of EVB
!

function dihed(at1,att2,at3,at4,xyz5)
use evb_mod
implicit none
integer::at1,att2,at3,at4
integer::i,j
integer::natms
real(kind=8)::dihed
real(kind=8)::ax,ay,az,bx,by,bz,cx,cy,cz
real(kind=8)::xyz5(3,natoms)
real(kind=8)::cos_val,nan,nbn,valijk,vecnorm_loc
real(kind=8)::ra(3),rb(3),rc(3),na(3),nb(3)
real(kind=8)::axb(3),bxc(3),abxbc(3)
real(kind=8)::deter,eps,snanb,thab,thbc
real(kind=8)::sin_val
integer::ic

real(kind=8)::norm
real(kind=8)::q1(3),q2(3),q3(3)
real(kind=8)::q12(3),q23(3),n1(3),n2(3),u1(3),u2(3),u3(3)

!     Difference vectors (involved bonds)
real(kind=8)::u_vec(3),v_vec(3),w_vec(3),u_len,v_len,w_len
!     Normalized difference vectors 
real(kind=8)::u_norm(3),v_norm(3),w_norm(3)
!     Cross products of involved bonds
real(kind=8)::uxw(3),vxw(3)
real(kind=8)::sin_phi_u,sin_phi_v
!     Actual atom indices
integer::atm1,atm2,atm3,atm4


!
!     The following strategy to calculate a dihedral is taken from:
!     http://azevedolab.net/resources/dihedral_angle.pdf
!     (also stored into literatur/utilities)
!
!     1. calculate differece vectors of the points
!        P1,P2,P3 and P4
!
!q1(1)=xyz5(1,att2)-xyz5(1,at1)
!q1(2)=xyz5(2,att2)-xyz5(2,at1)
!q1(3)=xyz5(3,att2)-xyz5(3,at1)
!q2(1)=xyz5(1,at3)-xyz5(1,att2)
!q2(2)=xyz5(2,at3)-xyz5(2,att2)
!q2(3)=xyz5(3,at3)-xyz5(3,att2)
!q3(1)=xyz5(1,at4)-xyz5(1,at3)
!q3(2)=xyz5(2,at4)-xyz5(2,at3)
!q3(3)=xyz5(3,at4)-xyz5(3,at3)
!
!     2. calulate the cross vectors
!
!call crossprod(q1,q2,q12)
!call crossprod(q2,q3,q23)
!
!     3. calculate vectors n1 and n2 normal to planes defined by 
!         points P1,P2,P3,P4
!
!norm=vecnorm_loc(q12,3,1)
!n1=q12
!norm=vecnorm_loc(q23,3,1)
!n2=q23

!
!     4. Calculate orthogonal unit vectors 
!
!u1=n2
!norm=vecnorm_loc(q2,3,1)
!u3=q2
!call crossprod(u3,u1,u2)

!
!     5. calculate dihedral angle
!
!cos_val=dot_product(n1,u1)
!sin_val=dot_product(n1,u2)

!dihed=-atan2(sin_val,cos_val)
!if (dihed .ne. dihed) then
!   write(*,*) "Error in dihedral angle calculation!"
!   write(*,*) "Restart the current trajectory"
!   traj_error=1
!end if
!write(*,*) dihed
!
!      New: 
!
dihed=0.d0
atm1=at1
atm2=att2
atm3=at3
atm4=at4
u_vec(1)=xyz5(1,atm1)-xyz5(1,atm2)
u_vec(2)=xyz5(2,atm1)-xyz5(2,atm2)
u_vec(3)=xyz5(3,atm1)-xyz5(3,atm2)
v_vec(1)=xyz5(1,atm4)-xyz5(1,atm3)
v_vec(2)=xyz5(2,atm4)-xyz5(2,atm3)
v_vec(3)=xyz5(3,atm4)-xyz5(3,atm3)
w_vec(1)=xyz5(1,atm3)-xyz5(1,atm2)
w_vec(2)=xyz5(2,atm3)-xyz5(2,atm2)
w_vec(3)=xyz5(3,atm3)-xyz5(3,atm2)

u_len=sqrt(dot_product(u_vec,u_vec))
v_len=sqrt(dot_product(v_vec,v_vec))
w_len=sqrt(dot_product(w_vec,w_vec))

u_vec=u_vec/u_len
v_vec=v_vec/v_len
w_vec=w_vec/w_len


call crossprod(u_vec,w_vec,uxw)
call crossprod(v_vec,w_vec,vxw)
sin_phi_u=sqrt(1.d0-dot_product(u_vec,w_vec)*dot_product(u_vec,w_vec))
sin_phi_v=sqrt(1.d0-dot_product(v_vec,w_vec)*dot_product(v_vec,w_vec))


!dihed=sin_phi_u*sin_phi_v
!dihed=dot_product(uxw,vxw)
cos_val=dot_product(uxw,vxw)/(sin_phi_u*sin_phi_v)
!!write(*,*) acos(cos_val)
if (cos_val .ge. 1.d0) then
   cos_val=1.d0
else if (cos_val .le. -1.d0) then!   cos_val=-1.d0
   cos_val=-1.d0
end if
dihed=acos(cos_val)
!dihed=u_vec(1)*v_vec(1)
!write(*,*) "product_num",dot_product(u_vec,w_vec)**2
!write(*,*) "square",dot_product(u_norm,w_norm)**2
!dihed=sin_phi_u*sin_phi_v+dot_product(uxw,vxw)

!dihed=sqrt(1.d0-dot_product(u_norm,w_norm)**2)
!stop "Jjgf"
return
end function dihed
