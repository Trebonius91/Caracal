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
!     subroutine calc_dwilson: Calculate the derivative of the Wilson B matrix 
!     with respect to the cartesian coordinates needed for coordinate transfor-
!     maton of the hessian matrix (internal <--> cartesian).
!     Optionally, the matrix can be calculated analytially or numerically.
! 
!     part of EVB
!

subroutine calc_dwilson(xyz5,internal,dB_mat)
use evb_mod
implicit none
!     The actual structure in cartesian coordinates
real(kind=8), intent(in) :: xyz5(3,natoms)
!     The actual structure in internal coordinates
real(kind=8), intent(in) :: internal(nat6)
!     The returned derivative of the Wilson matrix (result) 
real(kind=8), intent(out) :: dB_mat(nat6,3*natoms,3*natoms)

!     The numerical derivative of the wilson
real(kind=8) :: dB_num(nat6,3*natoms,3*natoms)
!     Loop indices
integer::i,j,k,l,m,o,p,r,s
!     Number of atoms in structure 
integer::nat,nat3 
!     Temporarily stored cartesian coordinates
real(kind=8)::xyz_tmp(3,natoms)
!     component differences (cartesian) 
real(kind=8)::ax,ay,az,bx,by,bz
!     Shift and results of coordinate routines
real(kind=8)::shift,dist,ang,dihed,oop
!     For numerical calculation
real(kind=8)::hi,lo,coord
!     Formula parts for analytical calculation
real(kind=8)::denom,deriv1,deriv2
real(kind=8)::arg,root2,root1,roots
!     Actual atom indices
integer::atm1,atm2,atm3,atm4
integer::pos1,pos2
real(kind=8)::actlen,actlen3,num
!     for numerical: storage of 2D gradients
real(kind=8)::difftmp(3*natoms,3*natoms,2)
!     For new analytical angle method
integer::at_act,sign_func
!     Kronecker deltas for second derivatives
real(kind=8)::krondel
!     Difference vectors (involved bonds)
real(kind=8)::u_vec(3),v_vec(3),w_vec(3),u_len,v_len,w_len
!     The vector length calculation function
real(kind=8)::vlen
!     Angle coordinate values 
real(kind=8)::angval,angvec(3)
!     Prefactors for angle derivatives 
real(kind=8)::prefac1,prefac2
!     Auxiliary for angles: bond first derivative matrices
real(kind=8)::du_mat(3,3),dv_mat(3,3),dw_mat(3,3)
!     Auxiliary: bond second derivative matrices 
real(kind=8)::d2u_mat(3,3,3),d2v_mat(3,3,3),d2w_mat(3,3,3)
!     Cross products of involved bonds
real(kind=8)::uxw(3),wxv(3),vxw(3),uxv(3),wxu(3)
!     Dot products of involved bonds 
real(kind=8)::udotw,vdotw
!     Cross products of bond derivatives 
real(kind=8)::duxw(3,3),uxdw(3,3),dvxw(3,3),vxdw(3,3)
!     Cross products of second bond derivatives 
real(kind=8)::d2uxw(3,3,3),d2vxw(3,3,3),uxd2w(3,3,3),vxd2w(3,3,3)
real(kind=8)::duxdw(3,3,3),dvxdw(3,3,3)
!     Bend angle sines and cosines
real(kind=8)::sin_phi,cos_phi
!     Bend angle values for dihedral calculation
real(kind=8)::sin_phi_u,sin_phi_v
!     factors for dihedral calculation
real(kind=8)::divide
real(kind=8)::dnumdai,dnumdbj,ddenomdai,ddenomdbj,d2num,d2denom
!     matrix with first derivatives for bend angle
real(kind=8)::dbends(3,3)
!     Angular values for out of plane coordinate
real(kind=8)::cosphi1,cosphi2,cosphi3
!     Derivatives of single angles in out of plane 
real(kind=8)::dcosphi1(3),dcosphi2(3),dcosphi3(3)
!     Plane normation factor h for out of plane 
real(kind=8)::normvec(3),oopnorm
!     Cross products for oop coordinate 
real(kind=8)::uv_vec(3),vw_vec(3),wu_vec(3)
!     Out of plane value for derivative calculation
real(kind=8)::psi
!     First derivative of out of plane normation factor h to bondlengths 
real(kind=8)::dh_mat(3,3),auxcos(3)
!     Parts of the total second out of plane derivative 
real(kind=8)::d2norm,dinorm,djnorm
real(kind=8)::prod,diprod,djprod,d2prod
!     The Pi
real(kind=8),parameter :: pi=3.141592653589793238462643383d0

!
!     Restore the number of atoms
!
nat=natoms
nat3=3*nat
!
!     Fill matrix with zeros
!
dB_mat=0.d0
dB_num=0.d0
!
!     OPTION 1: Calculate the numerical Wilson matrix
!
if (num_wilson .eqv. .true.) then
!
!     Save the input structure 
!
   xyz_tmp=xyz5
!
!     Fill the matrix with zeros
!
   shift=0.0001d0

   do i=1,nat6
      l=0
      do j=1,nat3
         m=j-((j-1)/3*3)
         if (mod(j+2,3) .eq. 0) then
            l=l+1
         end if
         do o=1,2
            if (o.eq.1) then
               xyz_tmp(m,l)=xyz_tmp(m,l)+shift
            else if (o.eq.2) then
               xyz_tmp(m,l)=xyz_tmp(m,l)-2.d0*shift
            end if
            r=0
            do k=1,nat3
               s=k-((k-1)/3*3)
               if (mod(k+2,3) .eq. 0) then
                  r=r+1
               end if
               do p=1,2
                  if (p .eq. 1) then
                     xyz_tmp(s,r)=xyz_tmp(s,r)+shift
                  else
                     xyz_tmp(s,r)=xyz_tmp(s,r)-2.d0*shift
                  end if
                  if (coord_def(i,1) .eq. 1) then
                     coord=dist(coord_def(i,2),coord_def(i,3),xyz_tmp)
                  else if (coord_def(i,1) .eq. 2) then
                     coord=ang(coord_def(i,2),coord_def(i,3),coord_def(i,4),xyz_tmp)
                  else if (coord_def(i,1) .eq. 3) then
                     coord=dihed(coord_def(i,2),coord_def(i,3),coord_def(i,4), & 
                            & coord_def(i,5),xyz_tmp)
                  else if (coord_def(i,1) .eq. 4) then
                     coord=oop(coord_def(i,2),coord_def(i,3),coord_def(i,4), &
                            & coord_def(i,5),xyz_tmp)
                  end if
                  if (p.eq.1) then
                     hi=coord
                  else if (p.eq.2) then
                     lo=coord
                  end if
               end do
               xyz_tmp(s,r)=xyz_tmp(s,r)+shift
               difftmp(j,k,o)=(hi-lo)/(2.d0*shift)
            end do
         end do
         xyz_tmp(m,l)=xyz_tmp(m,l)+shift
         xyz_tmp=xyz5
      end do
      do j=1,nat3
         do k=1,nat3
!            write(*,*) i,j,k,(difftmp(j,k,1)-difftmp(j,k,2))/(2.d0*shift)
            dB_mat(i,j,k)=(difftmp(j,k,1)-difftmp(j,k,2))/(2.d0*shift)
         end do
      end do
   end do
   db_num=dB_mat
!
!     TEST 26.11.2018: calculate Wilson matrix derivative numerically
!       in order to avoid hessian explosions..
!
!   return

else 
!
!    OPTION 2: calculate it analytically!
!
!    outer loop: internal coordinates
!    inner loop(s): cartesian coordinates
!
   db_mat=0.d0
   do i=1,nat6
!
!     1: BOND LENGTHS
!
      if (coord_def(i,1) .eq. 1) then
!
!     Calculate entries of involved bond vector, take its 
!     absolute value from the coordinate array
!
         u_len=internal(i)
         atm1=coord_def(i,2)
         atm2=coord_def(i,3)
         u_vec(1)=xyz5(1,atm1)-xyz5(1,atm2)
         u_vec(2)=xyz5(2,atm1)-xyz5(2,atm2)
         u_vec(3)=xyz5(3,atm1)-xyz5(3,atm2)
         u_vec=u_vec/u_len
!
!     Now calculate the values of all second derivatives for involved atoms
!
         do j=1,2
            pos1=coord_def(i,j+1)
            do k=1,2
               pos2=coord_def(i,k+1)
               do l=1,3
                  do m=1,3

                     dB_mat(i,(pos1-1)*3+l,(pos2-1)*3+m)=(-1)**krondel(pos1,pos2)*&
                               & (u_vec(l)*u_vec(m)-krondel(l,m))/u_len
                  end do
               end do
            end do
         end do
!
!     2: BOND ANGLES
!
      else if (coord_def(i,1) .eq. 2) then
!
!     Calculate entries of involved bond vectors
!
         atm1=coord_def(i,2)
         atm2=coord_def(i,3)
         atm3=coord_def(i,4)
         u_vec(1)=xyz5(1,atm1)-xyz5(1,atm2)
         u_vec(2)=xyz5(2,atm1)-xyz5(2,atm2)
         u_vec(3)=xyz5(3,atm1)-xyz5(3,atm2)
         v_vec(1)=xyz5(1,atm3)-xyz5(1,atm2)
         v_vec(2)=xyz5(2,atm3)-xyz5(2,atm2)
         v_vec(3)=xyz5(3,atm3)-xyz5(3,atm2)

         u_len=vlen(u_vec)
         v_len=vlen(v_vec)

!
!     Auxiliary matrix for first derivatives of normalized bonds 
!
!
         du_mat=0.d0
         dv_mat=0.d0
         do j=1,3
            do k=1,3
               if (j .eq. k) then
                  do l=1,3
                     if (l .ne.j) then
                        du_mat(j,k)=du_mat(j,k)+u_vec(l)*u_vec(l)
                        dv_mat(j,k)=dv_mat(j,k)+v_vec(l)*v_vec(l)
                     end if
                  end do
               else
                  du_mat(j,k)=-u_vec(j)*u_vec(k)
                  dv_mat(j,k)=-v_vec(j)*v_vec(k)
               end if
            end do
         end do
         du_mat=du_mat/(u_len*u_len*u_len)
         dv_mat=dv_mat/(v_len*v_len*v_len)
!
!     Auxiliary matrix for second derivatives of normalized bonds 
!

         d2u_mat=0.d0
         call calc_d2vec(d2u_mat,u_vec,u_len)
         d2v_mat=0.d0
         call calc_d2vec(d2v_mat,v_vec,v_len)
!
!     Normalize vectors 
!
         u_vec=u_vec/u_len
         v_vec=v_vec/v_len
!
!     value of the angle coordinate 
!
         angvec=0.5d0*(u_vec(:)-v_vec(:))
         angval=vlen(angvec)
!
!     prefactors for all derivatives 
!
         prefac1=-1d0/(4d0*angval)
         prefac2=1d0/(4d0*angval*angval)
!
!     Now calculate the actual derivatives 
!
!        atom1,atom1
!            
         pos1=coord_def(i,2)
         pos2=pos1
         do j=1,3
            do k=1,3
               dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*&
                       & (prefac2*dot_product(du_mat(j,:),v_vec)*&
                       & dot_product(du_mat(k,:),v_vec)+dot_product(d2u_mat(j,k,:),v_vec))
            end do
         end do
!
!        atom1,atom2
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,3)
         do j=1,3
            do k=1,3
                dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(-prefac2*&
                       & dot_product(du_mat(j,:),v_vec)*(dot_product(du_mat(k,:),v_vec)+&
                       & dot_product(u_vec,dv_mat(k,:))) - dot_product(d2u_mat(j,k,:),v_vec)-&
                       & dot_product(du_mat(j,:),dv_mat(k,:)))
            end do
         end do
!
!        atom1,atom3
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,4)
         do j=1,3
            do k=1,3
               dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(prefac2*&
                       & dot_product(du_mat(j,:),v_vec)*dot_product(dv_mat(k,:),u_vec)+&
                       & dot_product(du_mat(j,:),dv_mat(k,:)))
            end do
         end do
!
!        atom2,atom1
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,2)
         do j=1,3
            do k=1,3
                dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(-prefac2*&
                       & dot_product(du_mat(k,:),v_vec)*(dot_product(du_mat(j,:),v_vec)+&
                       & dot_product(u_vec,dv_mat(j,:))) - dot_product(d2u_mat(j,k,:),v_vec)-&
                       & dot_product(du_mat(k,:),dv_mat(j,:)))
            end do
         end do

!
!        atom2,atom2
!
         pos1=coord_def(i,3)
         pos2=pos1

         do j=1,3
            do k=1,3
                dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(prefac2*&
                       & (dot_product(du_mat(j,:),v_vec)+dot_product(u_vec,dv_mat(j,:)))*&
                       & (dot_product(du_mat(k,:),v_vec)+dot_product(u_vec,dv_mat(k,:)))+&
                       & (dot_product(d2u_mat(j,k,:),v_vec)+dot_product(du_mat(j,:),dv_mat(k,:))+&
                       & dot_product(du_mat(k,:),dv_mat(j,:))+dot_product(d2v_mat(j,k,:),u_vec)))
            end do
         end do
!
!        atom2,atom3
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,4)
         do j=1,3
            do k=1,3
                dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(-prefac2*&
                       & dot_product(dv_mat(k,:),u_vec)*(dot_product(du_mat(j,:),v_vec)+&
                       & dot_product(u_vec,dv_mat(j,:))) - dot_product(d2v_mat(j,k,:),u_vec)-&
                       & dot_product(du_mat(j,:),dv_mat(k,:)))
            end do
         end do
!
!        atom3,atom1
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,2)
         do j=1,3
            do k=1,3
               dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(prefac2*&
                       & dot_product(du_mat(k,:),v_vec)*dot_product(dv_mat(j,:),u_vec)+&
                       & dot_product(du_mat(k,:),dv_mat(j,:)))
            end do
         end do
!
!        atom3,atom2
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,3)
         do j=1,3
            do k=1,3
                dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*(-prefac2*&
                       & dot_product(dv_mat(j,:),u_vec)*(dot_product(du_mat(k,:),v_vec)+&
                       & dot_product(u_vec,dv_mat(k,:))) - dot_product(d2v_mat(j,k,:),u_vec)-&
                       & dot_product(du_mat(k,:),dv_mat(j,:)))
            end do
         end do
!
!        atom3,atom3
! 
         pos1=coord_def(i,4)
         pos2=pos1
         do j=1,3
            do k=1,3
               dB_mat(i,(pos1-1)*3+j,(pos2-1)*3+k)=prefac1*&
                       & (prefac2*dot_product(dv_mat(j,:),u_vec)*&
                       & dot_product(dv_mat(k,:),u_vec)+dot_product(d2v_mat(j,k,:),u_vec))
            end do
         end do

!
!     3: TORSIONAL ANGLES
!
      else if (coord_def(i,1) .eq. 3) then

!
!     Calculate entries of involved bond vectors
!
         atm1=coord_def(i,2)
         atm2=coord_def(i,3)
         atm3=coord_def(i,4)
         atm4=coord_def(i,5)
         u_vec(1)=xyz5(1,atm1)-xyz5(1,atm2)
         u_vec(2)=xyz5(2,atm1)-xyz5(2,atm2)
         u_vec(3)=xyz5(3,atm1)-xyz5(3,atm2)
         v_vec(1)=xyz5(1,atm4)-xyz5(1,atm3)
         v_vec(2)=xyz5(2,atm4)-xyz5(2,atm3)
         v_vec(3)=xyz5(3,atm4)-xyz5(3,atm3)
         w_vec(1)=xyz5(1,atm3)-xyz5(1,atm2)
         w_vec(2)=xyz5(2,atm3)-xyz5(2,atm2)
         w_vec(3)=xyz5(3,atm3)-xyz5(3,atm2)

!
!      Calculate the bondlengths of relevant vectors 
!
         u_len=sqrt(dot_product(u_vec,u_vec))
         v_len=sqrt(dot_product(v_vec,v_vec))
         w_len=sqrt(dot_product(w_vec,w_vec))

!
!      Calculate auxiliary matrices for first bond derivatives 
!

         du_mat=0.d0
         dv_mat=0.d0
         dw_mat=0.d0
         do j=1,3
            do k=1,3
               if (j .eq. k) then
                  do l=1,3
                     if (l .ne.j) then
                        du_mat(j,k)=du_mat(j,k)+u_vec(l)*u_vec(l)
                        dv_mat(j,k)=dv_mat(j,k)+v_vec(l)*v_vec(l)
                        dw_mat(j,k)=dw_mat(j,k)+w_vec(l)*w_vec(l)
                     end if
                  end do
               else
                  du_mat(j,k)=-u_vec(j)*u_vec(k)
                  dv_mat(j,k)=-v_vec(j)*v_vec(k)
                  dw_mat(j,k)=-w_vec(j)*w_vec(k)
               end if
            end do
         end do
         du_mat=du_mat/(u_len*u_len*u_len)
         dv_mat=dv_mat/(v_len*v_len*v_len)
         dw_mat=dw_mat/(w_len*w_len*w_len)
!
!     Auxiliary matrix for second derivatives of normalized bonds 
!

         d2u_mat=0.d0
         call calc_d2vec(d2u_mat,u_vec,u_len)
         d2v_mat=0.d0
         call calc_d2vec(d2v_mat,v_vec,v_len)
         d2w_mat=0.d0
         call calc_d2vec(d2w_mat,w_vec,w_len)
!
!      Normalize bond vectors 
!
         u_vec=u_vec/u_len
         v_vec=v_vec/v_len
         w_vec=w_vec/w_len
!
!      calculate composite vectors needed for cross product derivatives 
!
         call crossprod(u_vec,w_vec,uxw)
         call crossprod(v_vec,w_vec,vxw)
         do k=1,3
            call crossprod(du_mat(k,:),w_vec,duxw(k,:))
            call crossprod(dv_mat(k,:),w_vec,dvxw(k,:))
            call crossprod(u_vec,dw_mat(k,:),uxdw(k,:))
            call crossprod(v_vec,dw_mat(k,:),vxdw(k,:))
         end do
!
!      composite vectors for second cross product derivatives
!
         do k=1,3
            do l=1,3
               call crossprod(d2u_mat(k,l,:),w_vec,d2uxw(k,l,:))
               call crossprod(u_vec,d2w_mat(k,l,:),uxd2w(k,l,:))
               call crossprod(v_vec,d2w_mat(k,l,:),vxd2w(k,l,:))
               call crossprod(d2v_mat(k,l,:),w_vec,d2vxw(k,l,:))

               call crossprod(du_mat(k,:),dw_mat(l,:),duxdw(k,l,:))
               call crossprod(dv_mat(k,:),dw_mat(l,:),dvxdw(k,l,:))
            end do
         end do
!
!      calculate dot products needed for dot denominator derivatives 
!
         udotw=dot_product(u_vec,w_vec)
         vdotw=dot_product(v_vec,w_vec)
!
!      calculate pseudo sinus expressions for both combinations of bond vectors 
!
         sin_phi_u=sqrt(1.d0-dot_product(u_vec,w_vec)*dot_product(u_vec,w_vec))
         sin_phi_v=sqrt(1.d0-dot_product(v_vec,w_vec)*dot_product(v_vec,w_vec))

!
!      calculate expressions for numerator and denominator (global)
!
         num=dot_product(uxw,vxw)
         denom=sin_phi_u*sin_phi_v
         divide=denom**3*(denom**2-num**2)*sqrt(1d0-(num/denom)**2)
!
!     Now calculate the actual derivatives 
!
!        atom1,atom1
!            
         pos1=coord_def(i,2)
         pos2=pos1
         do k=1,3
            dnumdai=dot_product(duxw(k,:),vxw(:))
            ddenomdai=-dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*sin_phi_v
            do l=1,3
               dnumdbj=dot_product(duxw(l,:),vxw(:))
               ddenomdbj=-dot_product(du_mat(l,:),w_vec)*udotw/sin_phi_u*sin_phi_v
               d2num=dot_product(d2uxw(k,l,:),vxw(:))
               d2denom=-(dot_product(d2u_mat(k,l,:),w_vec(:))*udotw+&
                     & dot_product(du_mat(k,:),w_vec(:))*dot_product(du_mat(l,:),w_vec(:)))*&
                     & sin_phi_v/sin_phi_u-(dot_product(du_mat(k,:),w_vec(:))*&
                     & dot_product(du_mat(l,:),w_vec(:))*udotw**2)*&
                     & sin_phi_v/sin_phi_u**3
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom1,atom2
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,3)
         do k=1,3
            dnumdai=dot_product(duxw(k,:),vxw(:))
            ddenomdai=-dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*sin_phi_v
            do l=1,3
               dnumdbj=-dot_product(duxw(l,:),vxw(:))-dot_product(uxdw(l,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(l,:))
               ddenomdbj=(dot_product(du_mat(l,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v)
               d2num=-dot_product(d2uxw(k,l,:),vxw(:))-&
                     & dot_product(duxdw(k,l,:),vxw(:))-dot_product(duxw(k,:),vxdw(l,:))
               d2denom= + (udotw*(dot_product(d2u_mat(k,l,:),w_vec)+&
                     & dot_product(du_mat(k,:),dw_mat(l,:)))+dot_product(du_mat(k,:),w_vec)* &
                     & (dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:))))/sin_phi_u*sin_phi_v+&
                     & udotw*dot_product(du_mat(k,:),w_vec)* &
                     &((dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u**3*sin_phi_v-&
                     & (dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v/sin_phi_u)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide            
             end do
         end do
!
!        atom1,atom3
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,4)
         do k=1,3
            dnumdai=dot_product(duxw(k,:),vxw(:))
            ddenomdai=-dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*sin_phi_v
            do l=1,3
               dnumdbj=dot_product(uxdw(l,:),vxw(:))+dot_product(uxw(:),vxdw(l,:))-&
                     & dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-(dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(l,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(l,:)))*vdotw)/sin_phi_v)
               d2num=+dot_product(duxdw(k,l,:),vxw(:))+&
                     & dot_product(duxw(k,:),vxdw(l,:))-dot_product(duxw(k,:),dvxw(l,:))
               d2denom=-((dot_product(du_mat(k,:),dw_mat(l,:))*udotw)+&
                     & dot_product(du_mat(k,:),w_vec)*dot_product(u_vec,dw_mat(l,:)))*sin_phi_v/sin_phi_u-&
                     & dot_product(du_mat(k,:),w_vec)*udotw*(dot_product(u_vec,dw_mat(l,:))*&
                     & udotw/sin_phi_u**3*sin_phi_v+&
                     & (dot_product(dv_mat(l,:),w_vec)-dot_product(v_vec,dw_mat(l,:)))*&
                     & vdotw/sin_phi_v/sin_phi_u)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom1,atom4
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,5)
         do k=1,3
            dnumdai=dot_product(duxw(k,:),vxw(:))
            ddenomdai=-dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*sin_phi_v
            do l=1,3
               dnumdbj=dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-dot_product(dv_mat(l,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
               d2num=dot_product(duxw(k,:),dvxw(l,:))
               d2denom=dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*&
                     & dot_product(dv_mat(l,:),w_vec)*vdotw/sin_phi_v
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom2,atom1
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,2)
         do k=1,3
            dnumdai=-dot_product(duxw(k,:),vxw(:))-dot_product(uxdw(k,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(k,:))
            ddenomdai=(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(duxw(l,:),vxw(:))
               ddenomdbj=-dot_product(du_mat(l,:),w_vec)*udotw/sin_phi_u*sin_phi_v
               d2num=-dot_product(d2uxw(k,l,:),vxw(:))-&
                     & dot_product(duxdw(l,k,:),vxw(:))-dot_product(duxw(l,:),vxdw(k,:))
               d2denom=(udotw*(dot_product(d2u_mat(l,k,:),w_vec)+&
                     & dot_product(du_mat(l,:),dw_mat(k,:)))+dot_product(du_mat(l,:),w_vec)* &
                     & (dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:))))/sin_phi_u*sin_phi_v+&
                     & udotw*dot_product(du_mat(l,:),w_vec)* &
                     &((dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u**3*sin_phi_v-&
                     & (dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v/sin_phi_u)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom2,atom2
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,3)
         do k=1,3
            dnumdai=-dot_product(duxw(k,:),vxw(:))-dot_product(uxdw(k,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(k,:))
            ddenomdai=(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=-dot_product(duxw(l,:),vxw(:))-dot_product(uxdw(l,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(l,:))
               ddenomdbj=(dot_product(du_mat(l,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v)
               d2num=dot_product(d2uxw(k,l,:),vxw(:))+&
                     & dot_product(duxdw(k,l,:),vxw(:))+dot_product(duxw(k,:),vxdw(l,:))+&
                     & dot_product(duxdw(l,k,:),vxw(:))+dot_product(uxd2w(k,l,:),vxw(:))+&
                     & dot_product(uxdw(k,:),vxdw(l,:))+dot_product(duxw(l,:),vxdw(k,:))+&
                     & dot_product(uxdw(l,:),vxdw(k,:))+dot_product(uxw(:),vxd2w(k,l,:))
!        derivative denominator:
!        first term: d2sin_phi_u/dn*sin_phi_v
               d2denom=(-(udotw*(dot_product(d2u_mat(k,l,:),w_vec)+&
                     & dot_product(du_mat(k,:),dw_mat(l,:))+dot_product(du_mat(l,:),dw_mat(k,:))+&
                     & dot_product(u_vec,d2w_mat(k,l,:)))+(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*((dot_product(u_vec,dw_mat(l,:))+ &
                     & dot_product(du_mat(l,:),w_vec))))/sin_phi_u-&
                     & (dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:)))*udotw*&
                     & (dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:)))*udotw/&
                     & (sin_phi_u**3))*sin_phi_v+&
!        second term: dsin_phi_v/dni*dsin_phi_u/dnj
                     & (dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*&
                     & ((dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v)+&
!        third term: d2sin_phi_v/dn*sin_phi_u
                     &(-(dot_product(v_vec,d2w_mat(k,l,:))*vdotw+&
                     & dot_product(v_vec,dw_mat(k,:))*&
                     & dot_product(v_vec,dw_mat(l,:)))/sin_phi_v- &
                     & vdotw**2*(dot_product(v_vec,dw_mat(k,:))*dot_product(v_vec,dw_mat(l,:)))/&
                     & sin_phi_v**3)*sin_phi_u+&
!        fourth term: dsin_phi_v/dnj*dsin_phi_u/dni
                     & (dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*&
                     & ((dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide            
            end do
         end do
!
!        atom2,atom3
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,4)
         do k=1,3
            dnumdai=-dot_product(duxw(k,:),vxw(:))-dot_product(uxdw(k,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(k,:))
            ddenomdai=(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(uxdw(l,:),vxw(:))+dot_product(uxw(:),vxdw(l,:))-&
                     & dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-(dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(l,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(l,:)))*vdotw)/sin_phi_v)
               d2num=-dot_product(duxdw(k,l,:),vxw(:))-&
                     & dot_product(duxw(k,:),vxdw(l,:))+dot_product(duxw(k,:),dvxw(l,:))-&
                     & dot_product(uxd2w(k,l,:),vxw(:))-dot_product(uxdw(k,:),vxdw(l,:))+&
                     & dot_product(uxdw(k,:),dvxw(l,:))-dot_product(uxdw(l,:),vxdw(k,:))-&
                     & dot_product(uxw(:),vxd2w(k,l,:))+dot_product(uxw(:),dvxdw(l,k,:))
!        derivative denominator:
!        first term: d2sin_phi_u/dodn*sin_phi_v
               d2denom=((udotw*(dot_product(du_mat(k,:),dw_mat(l,:))+&
                     & dot_product(u_vec,d2w_mat(l,k,:)))+dot_product(u_vec,dw_mat(l,:))*&
                     &(dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:))))/sin_phi_u+&
                     & dot_product(u_vec,dw_mat(l,:))*udotw**2*((dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:))))/sin_phi_u**3)*sin_phi_v-&
!        second term: dsin_phi_v/doi*dsin_phi_u/dnj
                     & dot_product(u_vec,dw_mat(l,:))*udotw*dot_product(v_vec,dw_mat(k,:))*&
                     & vdotw/sin_phi_u/sin_phi_v+ &
!        third term: d2sin_phi_v/dodn*sin_phi_u
                     & ((vdotw*(dot_product(-dv_mat(l,:),dw_mat(k,:))+&
                     & dot_product(v_vec,d2w_mat(l,k,:)))+&
                     & (dot_product(v_vec,dw_mat(k,:))*(-dot_product(dv_mat(l,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(l,:)))))/sin_phi_v+ &
                     & vdotw**2*(dot_product(-dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:)))*&
                     & (dot_product(v_vec,dw_mat(k,:)))/sin_phi_v**3)*sin_phi_u-&
!        fourth term: dsin_phi_u/doi*dsin_phi_v/dnj
                     & (vdotw*(dot_product(-dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:))))*&
                     & (udotw*(dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:))))/&
                     & sin_phi_u/sin_phi_v
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom2,atom4
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,5)
         do k=1,3
            dnumdai=-dot_product(duxw(k,:),vxw(:))-dot_product(uxdw(k,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(k,:))
            ddenomdai=(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-dot_product(dv_mat(l,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
               d2num=-dot_product(duxw(k,:),dvxw(l,:))-&
                     & dot_product(uxdw(k,:),dvxw(l,:))-dot_product(uxw(:),dvxdw(l,k,:))
               d2denom=((dot_product(dv_mat(l,:),dw_mat(k,:))*vdotw)+&
                     & dot_product(dv_mat(l,:),w_vec)*dot_product(v_vec,dw_mat(k,:)))*sin_phi_u/sin_phi_v+&
                     & dot_product(dv_mat(l,:),w_vec)*vdotw*(dot_product(v_vec,dw_mat(k,:))*&
                     & vdotw/sin_phi_v**3*sin_phi_u-&
                     & (dot_product(du_mat(k,:),w_vec)+dot_product(u_vec,dw_mat(k,:)))*&
                     & udotw/sin_phi_u/sin_phi_v)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide            
            end do
         end do
!
!        atom3,atom1
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,2)
         do k=1,3
            dnumdai=dot_product(uxdw(k,:),vxw(:))+dot_product(uxw(:),vxdw(k,:))-&
                     & dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-(dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(duxw(l,:),vxw(:))
               ddenomdbj=-dot_product(du_mat(l,:),w_vec)*udotw/sin_phi_u*sin_phi_v
               d2num=+dot_product(duxdw(l,k,:),vxw(:))+&
                     & dot_product(duxw(l,:),vxdw(k,:))-dot_product(duxw(l,:),dvxw(k,:))
               d2denom=-((dot_product(du_mat(l,:),dw_mat(k,:))*udotw)+&
                     & dot_product(du_mat(l,:),w_vec)*dot_product(u_vec,dw_mat(k,:)))*sin_phi_v/sin_phi_u-&
                     & dot_product(du_mat(l,:),w_vec)*udotw*(dot_product(u_vec,dw_mat(k,:))*&
                     & udotw/sin_phi_u**3*sin_phi_v+&
                     & (dot_product(dv_mat(k,:),w_vec)-dot_product(v_vec,dw_mat(k,:)))*&
                     & vdotw/sin_phi_v/sin_phi_u)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom3,atom2
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,3)
         do k=1,3
            dnumdai=dot_product(uxdw(k,:),vxw(:))+dot_product(uxw(:),vxdw(k,:))-&
                     & dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-(dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=-dot_product(duxw(l,:),vxw(:))-dot_product(uxdw(l,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(l,:))
               ddenomdbj=(dot_product(du_mat(l,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v)
               d2num=-dot_product(duxdw(l,k,:),vxw(:))-&
                     & dot_product(duxw(l,:),vxdw(k,:))+dot_product(duxw(l,:),dvxw(k,:))-&
                     & dot_product(uxd2w(l,k,:),vxw(:))-dot_product(uxdw(l,:),vxdw(k,:))+&
                     & dot_product(uxdw(l,:),dvxw(k,:))-dot_product(uxdw(k,:),vxdw(l,:))-&
                     & dot_product(uxw(:),vxd2w(l,k,:))+dot_product(uxw(:),dvxdw(k,l,:))
!        derivative denominator:
!        first term: d2sin_phi_u/dodn*sin_phi_v
               d2denom=((udotw*(dot_product(du_mat(l,:),dw_mat(k,:))+&
                     & dot_product(u_vec,d2w_mat(k,l,:)))+dot_product(u_vec,dw_mat(k,:))*&
                     &(dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:))))/sin_phi_u+&
                     & dot_product(u_vec,dw_mat(k,:))*udotw**2*((dot_product(du_mat(l,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(l,:))))/sin_phi_u**3)*sin_phi_v-&
!        second term: dsin_phi_v/doi*dsin_phi_u/dnj
                     & dot_product(u_vec,dw_mat(k,:))*udotw*dot_product(v_vec,dw_mat(l,:))*&
                     & vdotw/sin_phi_u/sin_phi_v+ &
!        third term: d2sin_phi_v/dodn*sin_phi_u
                     & ((vdotw*(dot_product(-dv_mat(k,:),dw_mat(l,:))+&
                     & dot_product(v_vec,d2w_mat(k,l,:)))+&
                     & (dot_product(v_vec,dw_mat(l,:))*(-dot_product(dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))))/sin_phi_v+ &
                     & vdotw**2*(dot_product(-dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:)))*&
                     & (dot_product(v_vec,dw_mat(l,:)))/sin_phi_v**3)*sin_phi_u-&
!        fourth term: dsin_phi_u/doi*dsin_phi_v/dnj
                     & (vdotw*(dot_product(-dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:))))*&
                     & (udotw*(dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:))))/&
                     & sin_phi_u/sin_phi_v
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom3,atom3
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,4)
         do k=1,3
            dnumdai=dot_product(uxdw(k,:),vxw(:))+dot_product(uxw(:),vxdw(k,:))-&
                     & dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-(dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(uxdw(l,:),vxw(:))+dot_product(uxw(:),vxdw(l,:))-&
                     & dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-(dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(l,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(l,:)))*vdotw)/sin_phi_v)
               d2num=dot_product(uxd2w(k,l,:),vxw(:))+&
                     & dot_product(uxdw(k,:),vxdw(l,:))-dot_product(uxdw(k,:),dvxw(l,:))+&
                     & dot_product(uxdw(l,:),vxdw(k,:))+dot_product(uxw(:),vxd2w(k,l,:))-&
                     & dot_product(uxw(:),dvxdw(l,k,:))-dot_product(uxdw(l,:),dvxw(k,:))-&
                     & dot_product(uxw(:),dvxdw(k,l,:))+dot_product(uxw(:),d2vxw(k,l,:))
!        derivative denominator:
!        first term: d2sin_phi_v/dn*sin_phi_u
               d2denom=(-(vdotw*(dot_product(d2v_mat(k,l,:),w_vec)+&
                     & dot_product(-dv_mat(k,:),dw_mat(l,:))+dot_product(-dv_mat(l,:),dw_mat(k,:))+&
                     & dot_product(v_vec,d2w_mat(k,l,:)))+(dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*((dot_product(v_vec,dw_mat(l,:))+ &
                     & dot_product(-dv_mat(l,:),w_vec))))/sin_phi_v-&
                     & (dot_product(-dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:)))*vdotw*&
                     & (dot_product(-dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:)))*vdotw/&
                     & (sin_phi_v**3))*sin_phi_u+&
!        second term: dsin_phi_v/dni*dsin_phi_u/dnj
                     & (dot_product(-dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:)))*vdotw/sin_phi_v*&
                     & ((dot_product(u_vec,dw_mat(l,:))*udotw)/sin_phi_u)+&
!        third term: d2sin_phi_u/dn*sin_phi_v
                     & (-(dot_product(u_vec,d2w_mat(k,l,:))*udotw+&
                     & dot_product(u_vec,dw_mat(k,:))*&
                     & dot_product(u_vec,dw_mat(l,:)))/sin_phi_u- &
                     & udotw**2*(dot_product(u_vec,dw_mat(k,:))*dot_product(u_vec,dw_mat(l,:)))/&
                     & sin_phi_u**3)*sin_phi_v+&
!        fourth term: dsin_phi_v/dnj*dsin_phi_u/dni
                     & (dot_product(-dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:)))*vdotw/sin_phi_v*&
                     & ((dot_product(u_vec,dw_mat(k,:))*udotw)/sin_phi_u)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide            
            end do
         end do
!
!        atom3,atom4
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,5)
         do k=1,3
            dnumdai=dot_product(uxdw(k,:),vxw(:))+dot_product(uxw(:),vxdw(k,:))-&
                     & dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-(dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*vdotw)/sin_phi_v)
            do l=1,3
               dnumdbj=dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-dot_product(dv_mat(l,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
               d2num=dot_product(uxdw(k,:),dvxw(l,:))+&
                     & dot_product(uxw(:),dvxdw(l,k,:))-dot_product(uxw(:),d2vxw(k,l,:))
               d2denom= - (vdotw*(-dot_product(d2v_mat(l,k,:),w_vec)+&
                     & dot_product(dv_mat(l,:),dw_mat(k,:)))+dot_product(dv_mat(l,:),w_vec)* &
                     & (-dot_product(dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:))))/sin_phi_v*sin_phi_u-&
                     & vdotw*dot_product(dv_mat(l,:),w_vec)* &
                     &((-dot_product(dv_mat(k,:),w_vec)+dot_product(v_vec,dw_mat(k,:)))*vdotw/sin_phi_v**3*sin_phi_u-&
                     & (dot_product(u_vec,dw_mat(k,:))*udotw)/sin_phi_u/sin_phi_v)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom4,atom1
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,2)
         do k=1,3
            dnumdai=dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
            do l=1,3
               dnumdbj=dot_product(duxw(l,:),vxw(:))
               ddenomdbj=-dot_product(du_mat(l,:),w_vec)*udotw/sin_phi_u*sin_phi_v
               d2num=dot_product(duxw(l,:),dvxw(k,:))
               d2denom=dot_product(du_mat(l,:),w_vec)*udotw/sin_phi_u*&
                     & dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom4,atom2
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,3)
         do k=1,3
            dnumdai=dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
            do l=1,3
               dnumdbj=-dot_product(duxw(l,:),vxw(:))-dot_product(uxdw(l,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(l,:))
               ddenomdbj=(dot_product(du_mat(l,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(l,:))*vdotw)/sin_phi_v)
               d2num=-dot_product(duxw(l,:),dvxw(k,:))-&
                     & dot_product(uxdw(l,:),dvxw(k,:))-dot_product(uxw(:),dvxdw(k,l,:))
               d2denom=((dot_product(dv_mat(k,:),dw_mat(l,:))*vdotw)+&
                     & dot_product(dv_mat(k,:),w_vec)*dot_product(v_vec,dw_mat(l,:)))*sin_phi_u/sin_phi_v+&
                     & dot_product(dv_mat(k,:),w_vec)*vdotw*(dot_product(v_vec,dw_mat(l,:))*&
                     & vdotw/sin_phi_v**3*sin_phi_u-&
                     & (dot_product(du_mat(l,:),w_vec)+dot_product(u_vec,dw_mat(l,:)))*&
                     & udotw/sin_phi_u/sin_phi_v)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom4,atom3
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,4)
         do k=1,3
            dnumdai=dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
            do l=1,3
               dnumdbj=dot_product(uxdw(l,:),vxw(:))+dot_product(uxw(:),vxdw(l,:))-&
                     & dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-(dot_product(u_vec,dw_mat(l,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(l,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(l,:)))*vdotw)/sin_phi_v)
               d2num=dot_product(uxdw(l,:),dvxw(k,:))+&
                     & dot_product(uxw(:),dvxdw(k,l,:))-dot_product(uxw(:),d2vxw(l,k,:))
               d2denom= - (vdotw*(-dot_product(d2v_mat(k,l,:),w_vec)+&
                     & dot_product(dv_mat(k,:),dw_mat(l,:)))+dot_product(dv_mat(k,:),w_vec)* &
                     & (-dot_product(dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:))))/sin_phi_v*sin_phi_u-&
                     & vdotw*dot_product(dv_mat(k,:),w_vec)* &
                     &((-dot_product(dv_mat(l,:),w_vec)+dot_product(v_vec,dw_mat(l,:)))*vdotw/sin_phi_v**3*sin_phi_u-&
                     & (dot_product(u_vec,dw_mat(l,:))*udotw)/sin_phi_u/sin_phi_v)
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do
!
!        atom4,atom4
!            
         pos1=coord_def(i,5)
         pos2=pos1
         do k=1,3
            dnumdai=dot_product(uxw(:),dvxw(k,:))
            ddenomdai=-dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
            do l=1,3
               dnumdbj=dot_product(uxw(:),dvxw(l,:))
               ddenomdbj=-dot_product(dv_mat(l,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
               d2num=dot_product(uxw(:),d2vxw(k,l,:))
               d2denom=-(dot_product(d2v_mat(k,l,:),w_vec(:))*vdotw+&
                     & dot_product(dv_mat(k,:),w_vec(:))*dot_product(dv_mat(l,:),w_vec(:)))*&
                     & sin_phi_u/sin_phi_v-(dot_product(dv_mat(k,:),w_vec(:))*&
                     & dot_product(dv_mat(l,:),w_vec(:))*vdotw**2)*&
                     & sin_phi_u/sin_phi_v**3
               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(num**3*ddenomdbj*ddenomdai-denom*num**3*d2denom+&
                     & denom**3*(dnumdbj*ddenomdai+ddenomdbj*dnumdai+num*d2denom)-denom**4*&
                     & d2num+denom**2*num*(-2d0*ddenomdbj*ddenomdai-dnumdbj*dnumdai+num*d2num))/divide
            end do
         end do


!
!     4: OUT OF PLANE ANGLES
!
      else if (coord_def(i,1) .eq. 4) then
!
!     Calculate entries of involved bond vectors
!
         atm1=coord_def(i,2)
         atm2=coord_def(i,3)
         atm3=coord_def(i,4)
         atm4=coord_def(i,5)
         u_vec(1)=xyz5(1,atm1)-xyz5(1,atm4)
         u_vec(2)=xyz5(2,atm1)-xyz5(2,atm4)
         u_vec(3)=xyz5(3,atm1)-xyz5(3,atm4)
         v_vec(1)=xyz5(1,atm2)-xyz5(1,atm4)
         v_vec(2)=xyz5(2,atm2)-xyz5(2,atm4)
         v_vec(3)=xyz5(3,atm2)-xyz5(3,atm4)
         w_vec(1)=xyz5(1,atm3)-xyz5(1,atm4)
         w_vec(2)=xyz5(2,atm3)-xyz5(2,atm4)
         w_vec(3)=xyz5(3,atm3)-xyz5(3,atm4)
!
!      Normalize vectors and calculate the bondlengths
!
         u_len=sqrt(dot_product(u_vec,u_vec))
         v_len=sqrt(dot_product(v_vec,v_vec))
         w_len=sqrt(dot_product(w_vec,w_vec))
!
!      Calculate auxiliary matrices for first bond derivatives 
!

         du_mat=0.d0
         dv_mat=0.d0
         dw_mat=0.d0
         do j=1,3
            do k=1,3
               if (j .eq. k) then
                  do l=1,3
                     if (l .ne.j) then
                        du_mat(j,k)=du_mat(j,k)+u_vec(l)*u_vec(l)
                        dv_mat(j,k)=dv_mat(j,k)+v_vec(l)*v_vec(l)
                        dw_mat(j,k)=dw_mat(j,k)+w_vec(l)*w_vec(l)
                     end if
                  end do
               else
                  du_mat(j,k)=-u_vec(j)*u_vec(k)
                  dv_mat(j,k)=-v_vec(j)*v_vec(k)
                  dw_mat(j,k)=-w_vec(j)*w_vec(k)
               end if
            end do
         end do
         du_mat=du_mat/(u_len*u_len*u_len)
         dv_mat=dv_mat/(v_len*v_len*v_len)
         dw_mat=dw_mat/(w_len*w_len*w_len)
!
!     Auxiliary matrix for second derivatives of normalized bonds 
!

         d2u_mat=0.d0
         call calc_d2vec(d2u_mat,u_vec,u_len)
         d2v_mat=0.d0
         call calc_d2vec(d2v_mat,v_vec,v_len)
         d2w_mat=0.d0
         call calc_d2vec(d2w_mat,w_vec,w_len)
!
!        normalize vectors for separate usage 
!
         u_vec=u_vec/u_len
         v_vec=v_vec/v_len
         w_vec=w_vec/w_len
!
!     calculate norm vector for coordinate (cross products)
!
         call crossprod(u_vec,v_vec,uxv)
         call crossprod(v_vec,w_vec,vxw)
         call crossprod(w_vec,u_vec,wxu)
!
!     calculate needed cross product derivatives 
!
         do k=1,3
            call crossprod(dv_mat(k,:),w_vec,dvxw(k,:))
            call crossprod(v_vec,dw_mat(k,:),vxdw(k,:))
         end do
!
!      composite vectors for second cross product derivatives
!
         do k=1,3
            do l=1,3
               call crossprod(d2u_mat(k,l,:),w_vec,d2uxw(k,l,:))
               call crossprod(u_vec,d2w_mat(k,l,:),uxd2w(k,l,:))
               call crossprod(v_vec,d2w_mat(k,l,:),vxd2w(k,l,:))
               call crossprod(d2v_mat(k,l,:),w_vec,d2vxw(k,l,:))

               call crossprod(du_mat(k,:),dw_mat(l,:),duxdw(k,l,:))
               call crossprod(dv_mat(k,:),dw_mat(l,:),dvxdw(k,l,:))
            end do
         end do
!
!     calculate the artifical vector for the normalization
!  
         normvec=uxv+vxw+wxu
!
!     calculate the length of this vector
!
         oopnorm=vlen(normvec)
!
!     calculate the oop coordinate (normalized triple product)
!
         psi=dot_product(u_vec,vw_vec)/oopnorm
!
!     calculate the three cosinus values needed for derivatives of oopnorm
!
         cosphi1=dot_product(w_vec,v_vec)
         cosphi2=dot_product(w_vec,u_vec)
         cosphi3=dot_product(u_vec,v_vec)
!
!     calculate derivatives of normation factors to the three involved bonds 
!
!     first: build auxiliary vector for derivatives 
!
         auxcos(1)=-cosphi1-1d0+cosphi2+cosphi3
         auxcos(2)=-cosphi2-1d0+cosphi1+cosphi3
         auxcos(3)=-cosphi3-1d0+cosphi2+cosphi1

         do k=1,3
            dcosphi2(k)=dot_product(w_vec,du_mat(k,:))
            dcosphi3(k)=dot_product(v_vec,du_mat(k,:))
            dh_mat(1,k)=(dcosphi2(k)*auxcos(2)+dcosphi3(k)*auxcos(3))
            dcosphi1(k)=dot_product(w_vec,dv_mat(k,:))
            dcosphi3(k)=dot_product(u_vec,dv_mat(k,:))
            dh_mat(2,k)=(dcosphi1(k)*auxcos(1)+dcosphi3(k)*auxcos(3))
            dcosphi1(k)=dot_product(v_vec,dw_mat(k,:))
            dcosphi2(k)=dot_product(u_vec,dw_mat(k,:))
            dh_mat(3,k)=(dcosphi1(k)*auxcos(1)+dcosphi2(k)*auxcos(2))
         end do
!
!     The triple product of the original oop formulation
!
         prod=dot_product(u_vec,vxw)
!
!     Now calculate the actual derivatives 
!
!        atom1,atom1
!           

         pos1=coord_def(i,2)
         pos2=pos1
         do k=1,3
            dinorm=dh_mat(1,k)/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)
            do l=1,3
               djnorm=dh_mat(1,l)/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)
               d2norm=3d0*dh_mat(1,k)*dh_mat(1,l)/oopnorm**5-&
                  & ((dot_product(du_mat(k,:),w_vec)-dot_product(du_mat(k,:),v_vec))*&
                  & (-dot_product(du_mat(l,:),w_vec)+dot_product(du_mat(l,:),v_vec))+&
                  & dot_product(d2u_mat(k,l,:),w_vec)*auxcos(2)+dot_product(d2u_mat(k,l,:),v_vec)*auxcos(3))/&
                  & oopnorm**3
               d2prod=dot_product(d2u_mat(k,l,:),vxw)               

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom1,atom2
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,3)
         do k=1,3
            dinorm=dh_mat(1,k)/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)
            do l=1,3
               djnorm=dh_mat(2,l)/oopnorm**3
               djprod=dot_product(u_vec,dvxw(l,:))
               d2norm=3d0*dh_mat(1,k)*dh_mat(2,l)/oopnorm**5-&
                  & (dot_product(du_mat(k,:),dv_mat(l,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(k,:),w_vec)*(dot_product(dv_mat(l,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(l,:))+dot_product(dv_mat(l,:),w_vec)))/oopnorm**3
               d2prod=dot_product(du_mat(k,:),dvxw(l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom1,atom3
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,4)
         do k=1,3
            dinorm=dh_mat(1,k)/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)
            do l=1,3
               djnorm=dh_mat(3,l)/oopnorm**3
               djprod=dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*dh_mat(1,k)*dh_mat(3,l)/oopnorm**5-&
                  & (dot_product(du_mat(k,:),dw_mat(l,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(k,:),w_vec)*(-dot_product(u_vec,dw_mat(l,:))+&
                  & dot_product(v_vec,dw_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(l,:))+dot_product(v_vec,dw_mat(l,:))))/oopnorm**3
               d2prod=dot_product(du_mat(k,:),vxdw(l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom1,atom4
!
         pos1=coord_def(i,2)
         pos2=coord_def(i,5)
         do k=1,3
            dinorm=dh_mat(1,k)/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)
            do l=1,3
               djnorm=(dh_mat(1,l)+dh_mat(2,l)+dh_mat(3,l))/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)+dot_product(u_vec,dvxw(l,:))+&
                  & dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*(dh_mat(1,k)*dh_mat(1,l)+dh_mat(1,k)*dh_mat(2,l)+dh_mat(1,k)*dh_mat(3,l))/&
                  & oopnorm**5-&
                  & (((dot_product(du_mat(k,:),w_vec)-dot_product(du_mat(k,:),v_vec))*&
                  & (-dot_product(du_mat(l,:),w_vec)+dot_product(du_mat(l,:),v_vec))+&
                  & dot_product(d2u_mat(k,l,:),w_vec)*auxcos(2)+dot_product(d2u_mat(k,l,:),v_vec)*auxcos(3))+&
                  & (dot_product(du_mat(k,:),dv_mat(l,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(k,:),w_vec)*(dot_product(dv_mat(l,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(l,:))+dot_product(dv_mat(l,:),w_vec)))+ &
                  & (dot_product(du_mat(k,:),dw_mat(l,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(k,:),w_vec)*(-dot_product(u_vec,dw_mat(l,:))+&
                  & dot_product(v_vec,dw_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(l,:))+dot_product(v_vec,dw_mat(l,:)))))/oopnorm**3
               d2prod=dot_product(d2u_mat(k,l,:),vxw)+dot_product(du_mat(k,:),dvxw(l,:))+&
                  & dot_product(du_mat(k,:),vxdw(l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)

            end do
         end do
!
!        atom2,atom1
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,2)
         do k=1,3
            dinorm=dh_mat(2,k)/oopnorm**3
            diprod=dot_product(u_vec,dvxw(k,:))
            do l=1,3
               djnorm=dh_mat(1,l)/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)
               d2norm=3d0*dh_mat(1,l)*dh_mat(2,k)/oopnorm**5-&
                  & (dot_product(du_mat(l,:),dv_mat(k,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(l,:),w_vec)*(dot_product(dv_mat(k,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(k,:))+dot_product(dv_mat(k,:),w_vec)))/oopnorm**3
               d2prod=dot_product(du_mat(l,:),dvxw(k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm
            end do
         end do
!
!        atom2,atom2
!           

         pos1=coord_def(i,3)
         pos2=pos1
         do k=1,3
            dinorm=dh_mat(2,k)/oopnorm**3
            diprod=dot_product(u_vec,dvxw(k,:))
            do l=1,3
               djnorm=dh_mat(2,l)/oopnorm**3
               djprod=dot_product(u_vec,dvxw(l,:))
               d2norm=3d0*dh_mat(2,k)*dh_mat(2,l)/oopnorm**5-&
                  & ((dot_product(dv_mat(k,:),w_vec)-dot_product(dv_mat(k,:),u_vec))*&
                  & (-dot_product(dv_mat(l,:),w_vec)+dot_product(dv_mat(l,:),u_vec))+&
                  & dot_product(d2v_mat(k,l,:),w_vec)*auxcos(1)+dot_product(d2v_mat(k,l,:),u_vec)*auxcos(3))/&
                  & oopnorm**3
               d2prod=dot_product(u_vec,d2vxw(k,l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom2,atom3
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,4)
         do k=1,3
            dinorm=dh_mat(2,k)/oopnorm**3
            diprod=dot_product(u_vec,dvxw(k,:))
            do l=1,3
               djnorm=dh_mat(3,l)/oopnorm**3
               djprod=dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*dh_mat(2,k)*dh_mat(3,l)/oopnorm**5-&
                  & (dot_product(dv_mat(k,:),dw_mat(l,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(k,:),w_vec)*(-dot_product(dw_mat(l,:),v_vec)+&
                  & dot_product(dw_mat(l,:),u_vec))+dot_product(dv_mat(k,:),u_vec)*(&
                  & dot_product(dw_mat(l,:),v_vec)+dot_product(u_vec,dw_mat(l,:))))/oopnorm**3
               d2prod=dot_product(u_vec,dvxdw(k,l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm
            end do
         end do
!
!        atom2,atom4
!
         pos1=coord_def(i,3)
         pos2=coord_def(i,5)
         do k=1,3
            dinorm=dh_mat(2,k)/oopnorm**3
            diprod=dot_product(u_vec,dvxw(k,:))
            do l=1,3
               djnorm=(dh_mat(1,l)+dh_mat(2,l)+dh_mat(3,l))/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)+dot_product(u_vec,dvxw(l,:))+dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*(dh_mat(1,l)*dh_mat(2,k)+dh_mat(2,k)*dh_mat(2,l)+dh_mat(2,k)*&
                  & dh_mat(3,l))/oopnorm**5-&
                  & ((dot_product(du_mat(l,:),dv_mat(k,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(l,:),w_vec)*(dot_product(dv_mat(k,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(k,:))+dot_product(dv_mat(k,:),w_vec)))+&
                  & ((dot_product(dv_mat(k,:),w_vec)-dot_product(dv_mat(k,:),u_vec))*&
                  & (-dot_product(dv_mat(l,:),w_vec)+dot_product(dv_mat(l,:),u_vec))+&
                  & dot_product(d2v_mat(k,l,:),w_vec)*auxcos(1)+dot_product(d2v_mat(k,l,:),u_vec)*auxcos(3))+&
                  & (dot_product(dv_mat(k,:),dw_mat(l,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(k,:),w_vec)*(-dot_product(dw_mat(l,:),v_vec)+&
                  & dot_product(dw_mat(l,:),u_vec))+dot_product(dv_mat(k,:),u_vec)*(&
                  & dot_product(dw_mat(l,:),v_vec)+dot_product(u_vec,dw_mat(l,:)))))/oopnorm**3
               d2prod=dot_product(du_mat(l,:),dvxw(k,:))+dot_product(u_vec,d2vxw(k,l,:))+&
                  & dot_product(u_vec,dvxdw(k,l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)

            end do
         end do

!
!        atom3,atom1
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,2)
         do k=1,3
            dinorm=dh_mat(3,k)/oopnorm**3
            diprod=dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(1,l)/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)
               d2norm=3d0*dh_mat(1,l)*dh_mat(3,k)/oopnorm**5-&
                  & (dot_product(du_mat(l,:),dw_mat(k,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(l,:),w_vec)*(-dot_product(u_vec,dw_mat(k,:))+&
                  & dot_product(v_vec,dw_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(k,:))+dot_product(v_vec,dw_mat(k,:))))/oopnorm**3
               d2prod=dot_product(du_mat(l,:),vxdw(k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom3,atom2
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,3)
         do k=1,3
            dinorm=dh_mat(3,k)/oopnorm**3
            diprod=dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(2,l)/oopnorm**3
               djprod=dot_product(u_vec,dvxw(l,:))

               d2norm=3d0*dh_mat(2,l)*dh_mat(3,k)/oopnorm**5-&
                  & (dot_product(dv_mat(l,:),dw_mat(k,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(l,:),w_vec)*(-dot_product(dw_mat(k,:),v_vec)+&
                  & dot_product(dw_mat(k,:),u_vec))+dot_product(dv_mat(l,:),u_vec)*(&
                  & dot_product(dw_mat(k,:),v_vec)+dot_product(u_vec,dw_mat(k,:))))/oopnorm**3
               d2prod=dot_product(u_vec,dvxdw(l,k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm
            end do
         end do
!
!        atom3,atom3
!           

         pos1=coord_def(i,4)
         pos2=pos1
         do k=1,3
            dinorm=dh_mat(3,k)/oopnorm**3
            diprod=dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(3,l)/oopnorm**3
               djprod=dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*dh_mat(3,k)*dh_mat(3,l)/oopnorm**5-&
                  & ((dot_product(dw_mat(k,:),v_vec)-dot_product(dw_mat(k,:),u_vec))*&
                  & (-dot_product(dw_mat(l,:),v_vec)+dot_product(dw_mat(l,:),u_vec))+&
                  & dot_product(d2w_mat(k,l,:),v_vec)*auxcos(1)+dot_product(d2w_mat(k,l,:),u_vec)*auxcos(2))/&
                  & oopnorm**3
               d2prod=dot_product(u_vec,vxd2w(k,l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm

            end do
         end do
!
!        atom3,atom4
!
         pos1=coord_def(i,4)
         pos2=coord_def(i,5)
         do k=1,3
            dinorm=dh_mat(3,k)/oopnorm**3
            diprod=dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=(dh_mat(1,l)+dh_mat(2,l)+dh_mat(3,l))/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)+dot_product(u_vec,dvxw(l,:))+dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*(dh_mat(1,l)*dh_mat(3,k)+dh_mat(2,l)*dh_mat(3,k)+dh_mat(3,k)*dh_mat(3,l))/&
                  & oopnorm**5-&
                  & ((dot_product(du_mat(l,:),dw_mat(k,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(l,:),w_vec)*(-dot_product(u_vec,dw_mat(k,:))+&
                  & dot_product(v_vec,dw_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(k,:))+dot_product(v_vec,dw_mat(k,:))))+&
                  & (dot_product(dv_mat(l,:),dw_mat(k,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(l,:),w_vec)*(-dot_product(dw_mat(k,:),v_vec)+&
                  & dot_product(dw_mat(k,:),u_vec))+dot_product(dv_mat(l,:),u_vec)*(&
                  & dot_product(dw_mat(k,:),v_vec)+dot_product(u_vec,dw_mat(k,:))))+&
                  & ((dot_product(dw_mat(k,:),v_vec)-dot_product(dw_mat(k,:),u_vec))*&
                  & (-dot_product(dw_mat(l,:),v_vec)+dot_product(dw_mat(l,:),u_vec))+&
                  & dot_product(d2w_mat(k,l,:),v_vec)*auxcos(1)+&
                  & dot_product(d2w_mat(k,l,:),u_vec)*auxcos(2)))/oopnorm**3
               d2prod=dot_product(du_mat(l,:),vxdw(k,:))+dot_product(u_vec,dvxdw(l,k,:))+&
                  & dot_product(u_vec,vxd2w(k,l,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)

            end do
         end do
!
!        atom4,atom1
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,2)
         do k=1,3
            dinorm=(dh_mat(1,k)+dh_mat(2,k)+dh_mat(3,k))/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)+dot_product(u_vec,dvxw(k,:))+&
                  & dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(1,l)/oopnorm**3
               djprod=dot_product(du_mat(l,:),vxw)

               d2norm=3d0*(dh_mat(1,l)*dh_mat(1,k)+dh_mat(1,l)*dh_mat(2,k)+dh_mat(1,l)*dh_mat(3,k))/&
                  & oopnorm**5-&
                  & (((dot_product(du_mat(l,:),w_vec)-dot_product(du_mat(l,:),v_vec))*&
                  & (-dot_product(du_mat(k,:),w_vec)+dot_product(du_mat(k,:),v_vec))+&
                  & dot_product(d2u_mat(l,k,:),w_vec)*auxcos(2)+dot_product(d2u_mat(l,k,:),v_vec)*auxcos(3))+&
                  & (dot_product(du_mat(l,:),dv_mat(k,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(l,:),w_vec)*(dot_product(dv_mat(k,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(k,:))+dot_product(dv_mat(k,:),w_vec)))+ &
                  & (dot_product(du_mat(l,:),dw_mat(k,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(l,:),w_vec)*(-dot_product(u_vec,dw_mat(k,:))+&
                  & dot_product(v_vec,dw_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(k,:))+dot_product(v_vec,dw_mat(k,:)))))/oopnorm**3
               d2prod=dot_product(d2u_mat(l,k,:),vxw)+dot_product(du_mat(l,:),dvxw(k,:))+&
                  & dot_product(du_mat(l,:),vxdw(k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)
            end do
         end do
!
!        atom4,atom2
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,3)
         do k=1,3
            dinorm=(dh_mat(1,k)+dh_mat(2,k)+dh_mat(3,k))/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)+dot_product(u_vec,dvxw(k,:))+dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(2,l)/oopnorm**3
               djprod=dot_product(u_vec,dvxw(l,:))
               d2norm=3d0*(dh_mat(1,k)*dh_mat(2,l)+dh_mat(2,l)*dh_mat(2,k)+dh_mat(2,l)*&
                  & dh_mat(3,k))/oopnorm**5-&
                  & ((dot_product(du_mat(k,:),dv_mat(l,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(k,:),w_vec)*(dot_product(dv_mat(l,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(l,:))+dot_product(dv_mat(l,:),w_vec)))+&
                  & ((dot_product(dv_mat(l,:),w_vec)-dot_product(dv_mat(l,:),u_vec))*&
                  & (-dot_product(dv_mat(k,:),w_vec)+dot_product(dv_mat(k,:),u_vec))+&
                  & dot_product(d2v_mat(l,k,:),w_vec)*auxcos(1)+dot_product(d2v_mat(l,k,:),u_vec)*auxcos(3))+&
                  & (dot_product(dv_mat(l,:),dw_mat(k,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(l,:),w_vec)*(-dot_product(dw_mat(k,:),v_vec)+&
                  & dot_product(dw_mat(k,:),u_vec))+dot_product(dv_mat(l,:),u_vec)*(&
                  & dot_product(dw_mat(k,:),v_vec)+dot_product(u_vec,dw_mat(k,:)))))/oopnorm**3
               d2prod=dot_product(du_mat(k,:),dvxw(l,:))+dot_product(u_vec,d2vxw(l,k,:))+&
                  & dot_product(u_vec,dvxdw(l,k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)
            end do
         end do
!
!        atom4,atom3
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,4)
         do k=1,3
            dinorm=(dh_mat(1,k)+dh_mat(2,k)+dh_mat(3,k))/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)+dot_product(u_vec,dvxw(k,:))+dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=dh_mat(3,l)/oopnorm**3
               djprod=dot_product(u_vec,vxdw(l,:))
               d2norm=3d0*(dh_mat(1,k)*dh_mat(3,l)+dh_mat(2,k)*dh_mat(3,l)+dh_mat(3,l)*dh_mat(3,k))/&
                  & oopnorm**5-&
                  & ((dot_product(du_mat(k,:),dw_mat(l,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(k,:),w_vec)*(-dot_product(u_vec,dw_mat(l,:))+&
                  & dot_product(v_vec,dw_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(l,:))+dot_product(v_vec,dw_mat(l,:))))+&
                  & (dot_product(dv_mat(k,:),dw_mat(l,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(k,:),w_vec)*(-dot_product(dw_mat(l,:),v_vec)+&
                  & dot_product(dw_mat(l,:),u_vec))+dot_product(dv_mat(k,:),u_vec)*(&
                  & dot_product(dw_mat(l,:),v_vec)+dot_product(u_vec,dw_mat(l,:))))+&
                  & ((dot_product(dw_mat(l,:),v_vec)-dot_product(dw_mat(l,:),u_vec))*&
                  & (-dot_product(dw_mat(k,:),v_vec)+dot_product(dw_mat(k,:),u_vec))+&
                  & dot_product(d2w_mat(l,k,:),v_vec)*auxcos(1)+&
                  & dot_product(d2w_mat(l,k,:),u_vec)*auxcos(2)))/oopnorm**3
               d2prod=dot_product(du_mat(k,:),vxdw(l,:))+dot_product(u_vec,dvxdw(k,l,:))+&
                  & dot_product(u_vec,vxd2w(l,k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=-(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)

            end do
         end do
!
!        atom4,atom4
!
         pos1=coord_def(i,5)
         pos2=coord_def(i,5)
         do k=1,3
            dinorm=(dh_mat(1,k)+dh_mat(2,k)+dh_mat(3,k))/oopnorm**3
            diprod=dot_product(du_mat(k,:),vxw)+dot_product(u_vec,dvxw(k,:))+dot_product(u_vec,vxdw(k,:))
            do l=1,3
               djnorm=(dh_mat(1,l)+dh_mat(2,l)+dh_mat(3,l))/oopnorm**3
               djprod=dot_product(u_vec,vxdw(l,:))+dot_product(du_mat(l,:),vxw)+dot_product(u_vec,dvxw(l,:))
               d2norm=3d0*(dh_mat(1,k)*dh_mat(3,l)+dh_mat(2,k)*dh_mat(3,l)+dh_mat(3,l)*&
                  & dh_mat(3,k)+dh_mat(1,l)*dh_mat(1,k)+dh_mat(1,l)*dh_mat(2,k)+&
                  & dh_mat(1,l)*dh_mat(3,k)+dh_mat(1,k)*dh_mat(2,l)+dh_mat(2,l)*dh_mat(2,k)+&
                  & dh_mat(2,l)*dh_mat(3,k))/&
                  & oopnorm**5-&
                  & (((dot_product(du_mat(k,:),dw_mat(l,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(k,:),w_vec)*(-dot_product(u_vec,dw_mat(l,:))+&
                  & dot_product(v_vec,dw_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(l,:))+dot_product(v_vec,dw_mat(l,:))))+&
                  & (dot_product(dv_mat(k,:),dw_mat(l,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(k,:),w_vec)*(-dot_product(dw_mat(l,:),v_vec)+&
                  & dot_product(dw_mat(l,:),u_vec))+dot_product(dv_mat(k,:),u_vec)*(&
                  & dot_product(dw_mat(l,:),v_vec)+dot_product(u_vec,dw_mat(l,:))))+&
                  & ((dot_product(dw_mat(l,:),v_vec)-dot_product(dw_mat(l,:),u_vec))*&
                  & (-dot_product(dw_mat(k,:),v_vec)+dot_product(dw_mat(k,:),u_vec))+&
                  & dot_product(d2w_mat(l,k,:),v_vec)*auxcos(1)+&
                  & dot_product(d2w_mat(l,k,:),u_vec)*auxcos(2)))+&
                  & (((dot_product(du_mat(l,:),w_vec)-dot_product(du_mat(l,:),v_vec))*&
                  & (-dot_product(du_mat(k,:),w_vec)+dot_product(du_mat(k,:),v_vec))+&
                  & dot_product(d2u_mat(l,k,:),w_vec)*auxcos(2)+dot_product(d2u_mat(l,k,:),v_vec)*auxcos(3))+&
                  & (dot_product(du_mat(l,:),dv_mat(k,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(l,:),w_vec)*(dot_product(dv_mat(k,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(k,:))+dot_product(dv_mat(k,:),w_vec)))+ &
                  & (dot_product(du_mat(l,:),dw_mat(k,:))*(-cosphi2-1d0+cosphi1+cosphi3)+&
                  & dot_product(du_mat(l,:),w_vec)*(-dot_product(u_vec,dw_mat(k,:))+&
                  & dot_product(v_vec,dw_mat(k,:)))+dot_product(du_mat(l,:),v_vec)*(&
                  & dot_product(u_vec,dw_mat(k,:))+dot_product(v_vec,dw_mat(k,:)))))+&
                  & ((dot_product(du_mat(k,:),dv_mat(l,:))*(-cosphi3-1d0+cosphi1+cosphi2)+&
                  & dot_product(du_mat(k,:),w_vec)*(dot_product(dv_mat(l,:),w_vec)+&
                  & dot_product(u_vec,dv_mat(l,:)))+dot_product(du_mat(k,:),v_vec)*(-&
                  & dot_product(u_vec,dv_mat(l,:))+dot_product(dv_mat(l,:),w_vec)))+&
                  & ((dot_product(dv_mat(l,:),w_vec)-dot_product(dv_mat(l,:),u_vec))*&
                  & (-dot_product(dv_mat(k,:),w_vec)+dot_product(dv_mat(k,:),u_vec))+&
                  & dot_product(d2v_mat(l,k,:),w_vec)*auxcos(1)+dot_product(d2v_mat(l,k,:),u_vec)*auxcos(3))+&
                  & (dot_product(dv_mat(l,:),dw_mat(k,:))*(-cosphi1-1d0+cosphi2+cosphi3)+&
                  & dot_product(dv_mat(l,:),w_vec)*(-dot_product(dw_mat(k,:),v_vec)+&
                  & dot_product(dw_mat(k,:),u_vec))+dot_product(dv_mat(l,:),u_vec)*(&
                  & dot_product(dw_mat(k,:),v_vec)+dot_product(u_vec,dw_mat(k,:))))))/oopnorm**3
               d2prod=dot_product(du_mat(k,:),vxdw(l,:))+dot_product(u_vec,dvxdw(k,l,:))+&
                  & dot_product(u_vec,vxd2w(l,k,:))+dot_product(d2u_mat(l,k,:),vxw)+&
                  & dot_product(du_mat(l,:),dvxw(k,:))+dot_product(du_mat(l,:),vxdw(k,:))+&
                  & dot_product(du_mat(k,:),dvxw(l,:))+dot_product(u_vec,d2vxw(l,k,:))+&
                  & dot_product(u_vec,dvxdw(l,k,:))

               dB_mat(i,(pos1-1)*3+k,(pos2-1)*3+l)=(-d2norm*prod+dinorm*djprod+djnorm*diprod-&
                         & d2prod/oopnorm)

            end do
         end do


      end if 
   end do
!write(*,*) "analytical:"
!do i=1,3*natoms
!   write(*,*) B_mat(:,i)
!end do
!write(*,*) "numerical, analytical"
!do k=1,nat6
!   do i=1,3*natoms
!      do j=1,3*natoms
!         write(*,*) k,i,j,dB_num(k,i,j),dB_mat(k,i,j),dB_num(k,i,j)/dB_mat(k,i,j)
!      end do
!   end do
!end do 
end if
!db_mat=db_num
!stop "In dwilseon!"
return
end subroutine calc_dwilson
