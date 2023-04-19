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
!     subroutine calc_wilson: Calculate the Wilson B matrix needed 
!     for coordinate transformaton (internal <--> cartesian).
!     Optionally, the matrix can be calculated analytially or numerically.
! 
!     part of EVB
!

subroutine calc_wilson(xyz5,internal,B_mat)
use evb_mod
implicit none
!     The actual structure in cartesian coordinates
real(kind=8), intent(in) :: xyz5(3,natoms)
!     The actual structure in internal coordinates
real(kind=8), intent(in) :: internal(nat6)
!     The returned Wilson matrix (result) 
real(kind=8), intent(out) :: B_mat(nat6,3*natoms)
!     Loop indices
integer::i,j,k,l,m
!     Number of atoms in structure 
integer::nat 
!     Active atoms for elongations in numerical calculations
integer::active(4)
!     Number of ative atoms and the currently active atom 
integer::act_num,act_atom
!     Temporarily stored cartesian coordinates
real(kind=8)::xyz_tmp(3,natoms)
!     Shift and results of coordinate routines
real(kind=8)::shift,dist,ang,dihed,oop
!     For numerical calculation
real(kind=8)::hi,lo,coord
!     Actual atom indices
integer::atm1,atm2,atm3,atm4
!     For new analytical angle method
integer::at_act,sign_func
!     Difference vectors (involved bonds)
real(kind=8)::u_vec(3),v_vec(3),w_vec(3),u_len,v_len,w_len
!     Function for vector lengths 
real(kind=8)::vlen
!     For correct sign of dihedral angle derivative 
real(kind=8)::cos_val,dihed_test,sign
!     Cross products of involved bonds
real(kind=8)::uxw(3),wxv(3),vxw(3)
!     Dot products of involved bonds 
real(kind=8)::udotw,vdotw
!     Cross products of bond derivatives 
real(kind=8)::duxw(3,3),uxdw(3,3),dvxw(3,3),vxdw(3,3)
!     Cross products for oop coordinate 
real(kind=8)::uv_vec(3),vw_vec(3),wu_vec(3)
!     Normalized bond vectors 
real(kind=8)::u_norm(3),v_norm(3),w_norm(3)
!     Numerator and denominator for dihedrals and of derivative
real(kind=8)::num,denom,divide
!     Derivatives of numerator and denominator for dihedrals 
real(kind=8)::diff_num,diff_denom
!     Prefactor for dihedral calculation
real(kind=8)::fact
!     Angular values for out of plane coordinate
real(kind=8)::cosphi1,cosphi2,cosphi3
!     Derivatives of single angles in out of plane 
real(kind=8)::dcosphi1(3),dcosphi2(3),dcosphi3(3)
!     For angle calculations 
real(kind=8)::angvec(3),angval,prefac
!     Auxiliary for angles: bond derivative matrices
real(kind=8)::du_mat(3,3),dv_mat(3,3),dw_mat(3,3)
!     Derivative of out of plane normation factor h to bondlengths 
real(kind=8)::dh_mat(3,3),auxcos(3)
!     Bend angle values for dihedral calculation
real(kind=8)::sin_phi_u,sin_phi_v,cos_phi_u,cos_phi_v
!     Additional values and results for OOP calculations
real(kind=8)::normvec(3),oopnorm,psi
!     The Pi
real(kind=8),parameter :: pi=3.141592653589793238462643383d0
!   TEST
real(kind=8):: B_num(nat6,3*natoms)
real(kind=8)::test
logical::bad
real(kind=8)::y_vec(3),y_len
!
!     Restore the number of atoms
!
nat=natoms
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
   B_mat=0.d0
   shift=0.001d0
   do i=1,nat6
      l=0
!
!     Select atoms for which elongations shall be done (active ones)
!
!     Bond lengths
      if (coord_def(i,1) .eq. 1) then
         act_num=2
         active(1:2)=coord_def(i,2:3)
!     Bend angles 
      else if (coord_def(i,1) .eq. 2) then 
         act_num=3
         active(1:3)=coord_def(i,2:4)
!     Dihedrals or out of planes
      else 
         act_num=4 
         active(1:4)=coord_def(i,2:5)
      end if
      act_atom=active(1)
      do j=1,3*act_num
         m=j-((j-1)/3*3)
         if (mod(j+2,3) .eq. 0) then
            l=l+1
            act_atom=active(l) 
         end if
         do k=1,2
            if (k .eq. 1) then
               xyz_tmp(m,act_atom)=xyz_tmp(m,act_atom)-shift
            else 
               xyz_tmp(m,act_atom)=xyz_tmp(m,act_atom)+2*shift
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
            if (k.eq.1) then
               lo=coord
            else 
               hi=coord
            end if
         end do
         xyz_tmp(m,act_atom)=xyz_tmp(m,act_atom)-shift
         B_mat(i,(act_atom-1)*3+m)=(hi-lo)/(2*shift) 
         xyz_tmp=xyz5
      end do
   end do
   B_num=B_mat
!   write(*,*) "Num-wilson:"
!   do i=1,3*natoms
!      write(*,*) b_mat(:,i)
!   end do
!   write(*,*) sum(b_mat)
!   stop "GUOgou"
else 
!
!    OPTION 2: calculate it analytically!
!
!    outer loop: internal coordinates
!    inner loop(s): cartesian coordinates
!
   B_mat=0.d0
   do i=1,nat6
!
!    Bondlengths
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
!     Finally calculate the values of the derivative for all involved atoms
!
         do j=1,2
            at_act=coord_def(i,j+1)
            do k=1,3
               B_mat(i,(at_act-1)*3+k)=sign_func(at_act,atm1,atm2)*u_vec(k)
            end do
         end do
!
!    Bend angles
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

         u_vec=u_vec!/u_len
         v_vec=v_vec!/v_len
!
!     value of the angle coordinate 
!
         angvec=0.5d0*(u_vec(:)/u_len-v_vec(:)/v_len)
         angval=vlen(angvec) 
!
!     auxiliary matrices for derivative 
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
         v_vec=v_vec/v_len
         u_vec=u_vec/u_len
!
!     fill Wilson matrix for each atom, respectively
! 
         prefac=1.d0/(4.d0*angval)

!     atom no. 1
!
         at_act=coord_def(i,2)
         do k=1,3
          !  B_mat(i,(at_act-1)*3+k)=-dot_product(du_mat(:,k),v_vec)
            B_mat(i,(at_act-1)*3+k)=-prefac*dot_product(du_mat(:,k),v_vec)
         end do
!
!     atom no. 2
!
         at_act=coord_def(i,3)
         do k=1,3
          !   B_mat(i,(at_act-1)*3+k)=(dot_product(du_mat(:,k),v_vec)+dot_product(dv_mat(:,k),u_vec))
            B_mat(i,(at_act-1)*3+k)=prefac*(dot_product(du_mat(:,k),v_vec)+dot_product(dv_mat(:,k),u_vec))
         end do
!
!     atom no. 2
!
         at_act=coord_def(i,4)
         do k=1,3
        !     B_mat(i,(at_act-1)*3+k)=-dot_product(dv_mat(:,k),u_vec)
            B_mat(i,(at_act-1)*3+k)=-prefac*dot_product(dv_mat(:,k),u_vec)
         end do

         
!----------------------------------------------------
!      PART 3 - TORSIONAL ANGLES
!----------------------------------------------------
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
!      Calculate auxiliary matrices for bond derivatives 
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
!      Calculate expression for numerator and denominator of original dihedral expression
!
         num=dot_product(uxw,vxw)
         denom=sin_phi_u*sin_phi_v
         divide=sqrt(1d0-(num/denom)**2)*denom*denom
!
!      Set derivative to zero if bad behavior might be expected 
!
         cos_val=num/denom
         if ((cos_val .ge. 1.d0) .or. (cos_val .le. -1.d0) .or. (denom .lt. 1E-10)) then
            fact=0.d0
         else 
            fact=1.d0
         end if
!
!        atom 1 
!
         at_act=coord_def(i,2)
         do k=1,3
            diff_num=dot_product(duxw(k,:),vxw(:))
            diff_denom=-dot_product(du_mat(k,:),w_vec)*udotw/sin_phi_u*sin_phi_v
            B_mat(i,(at_act-1)*3+k)=-(diff_num*denom-num*diff_denom)/(divide)
         end do
!
!        atom 2
!
         at_act=coord_def(i,3)
         do k=1,3
            diff_num=-dot_product(duxw(k,:),vxw(:))-dot_product(uxdw(k,:),vxw(:))-&
                     & dot_product(uxw(:),vxdw(k,:))
            diff_denom=(dot_product(du_mat(k,:),w_vec)+&
                     & dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v+&
                     & (sin_phi_u*(dot_product(v_vec,dw_mat(k,:))*vdotw)/sin_phi_v)
            B_mat(i,(at_act-1)*3+k)=-(diff_num*denom-num*diff_denom)/(divide)
         end do
!
!        atom 3 
!
         at_act=coord_def(i,4)
         do k=1,3
            diff_num=dot_product(uxdw(k,:),vxw(:))+dot_product(uxw(:),vxdw(k,:))-&
                     & dot_product(uxw(:),dvxw(k,:))
            diff_denom=-(dot_product(u_vec,dw_mat(k,:)))*udotw/sin_phi_u*sin_phi_v-&
                     & (sin_phi_u*((dot_product(-dv_mat(k,:),w_vec)+&
                     & dot_product(v_vec,dw_mat(k,:)))*vdotw)/sin_phi_v)
            B_mat(i,(at_act-1)*3+k)=-(diff_num*denom-num*diff_denom)/(divide)
         end do
!
!        atom 4
!
         at_act=coord_def(i,5)
         do k=1,3
            diff_num=dot_product(uxw(:),dvxw(k,:))
            diff_denom=-dot_product(dv_mat(k,:),w_vec)*vdotw/sin_phi_v*sin_phi_u
            B_mat(i,(at_act-1)*3+k)=-(diff_num*denom-num*diff_denom)/(divide)
         end do


!----------------------------------------------------
!      PART 4 - OUT OF PLANE ANGLES
!----------------------------------------------------
          

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
!      Calculate auxiliary matrices for bond derivatives 
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
!        normalize vectors for separate usage 
!
         u_vec=u_vec/u_len
         v_vec=v_vec/v_len
         w_vec=w_vec/w_len
!
!     calculate norm vector for coordinate (cross products)
!
         call crossprod(u_vec,v_vec,uv_vec)
         call crossprod(v_vec,w_vec,vw_vec)
         call crossprod(w_vec,u_vec,wu_vec)
!
!     calculate the artifical vector for the normalization
!  
         normvec=uv_vec+vw_vec+wu_vec
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
            dh_mat(1,k)=(dcosphi2(k)*auxcos(2)+&
                         &dcosphi3(k)*auxcos(3))/oopnorm
            dcosphi1(k)=dot_product(w_vec,dv_mat(k,:))
            dcosphi3(k)=dot_product(u_vec,dv_mat(k,:))
            dh_mat(2,k)=(dcosphi1(k)*auxcos(1)+&
                         &dcosphi3(k)*auxcos(3))/oopnorm
            dcosphi1(k)=dot_product(v_vec,dw_mat(k,:))
            dcosphi2(k)=dot_product(u_vec,dw_mat(k,:))
            dh_mat(3,k)=(dcosphi1(k)*auxcos(1)+&
                         &dcosphi2(k)*auxcos(2))/oopnorm
         end do         
!
!        atom 1 
!
         at_act=coord_def(i,2)
         do k=1,3
!        test: only rho
!         B_mat(i,(at_act-1)*3+k)=-dh_mat(1,k)/oopnorm**2
!         step 1:
   !       B_mat(i,(at_act-1)*3+k)=(dot_product(du_mat(k,:),w_vec)*(-cosphi2-1d0+cosphi1+cosphi3)+&
    !             & dot_product(du_mat(k,:),v_vec)*(-cosphi3-1d0+cosphi2+cosphi1))/oopnorm

            B_mat(i,(at_act-1)*3+k)=-1.d0/oopnorm*(-psi*dh_mat(1,k)+&
                   & dot_product(du_mat(k,:),vw_vec))
         end do
!
!       atom 2
!         
         at_act=coord_def(i,3)
         do k=1,3
            B_mat(i,(at_act-1)*3+k)=-1.d0/oopnorm*(-psi*dh_mat(2,k)+&
                   & dot_product(dv_mat(k,:),wu_vec))
         end do
!
!       atom 3
!         
         at_act=coord_def(i,4)
         do k=1,3
            B_mat(i,(at_act-1)*3+k)=-1.d0/oopnorm*(-psi*dh_mat(3,k)+&
                   & dot_product(dw_mat(k,:),uv_vec))
         end do
!
!       atom 4
!    
         at_act=coord_def(i,5)
         do k=1,3
            B_mat(i,(at_act-1)*3+k)=1.d0/oopnorm*(-psi*(dh_mat(1,k)+dh_mat(2,k)+&
                   & dh_mat(3,k))+dot_product(du_mat(k,:),vw_vec)+dot_product(dv_mat(k,:)&
                   & ,wu_vec)+dot_product(dw_mat(k,:),uv_vec))
         end do

      end if

   end do
!write(*,*) "analytical:"
!do i=1,3*natoms
!   write(*,*) B_mat(:,i)
!end do
!bad=.false.
!write(*,*) "numerical, analytical:"
!do i=1,nat6
!   do j=1,natoms 
!      do k=1,3
!         write(*,*)i,j, B_num(i,(j-1)*3+k),B_mat(i,(j-1)*3+k),B_num(i,(j-1)*3+k)/B_mat(i,(j-1)*3+k)
!      end do
!   end do
!end do

!do i=nat6-2,nat6
!   do j=1,4
!      do k=1,3
!         write(198,*) i,coord_def(i,j+1)*3-3+k,B_num(i,coord_def(i,j+1)*3-3+k),B_mat(i,coord_def(i,j+1)*3-3+k),&
!                     & B_num(i,coord_def(i,j+1)*3-3+k)-B_mat(i,coord_def(i,j+1)*3-3+k)
!         if (B_num(i,coord_def(i,j+1)*3-3+k)-B_mat(i,coord_def(i,j+1)*3-3+k) .gt. 1E-5) then
!            bad=.true.
!         end if
!      end do
!   end do
!end do
!if (bad) then
!   write(*,*) "false analytical wilson!"
!   call fatal
!end if
!write(*,*) internal(7)

end if

!stop "In wilseon!"
return
end subroutine calc_wilson
