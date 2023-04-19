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
!     subroutine sum_dv12: calculate analytical gradient of DG-EVB 
!          coupling strengh for a given structure 
!
!     part of EVB
!

subroutine sum_dv12(dg_evb_mode,mat_size,alph,act_struc,g_V12)
use evb_mod
implicit none
integer::i,j,k,l,m,inc,block,first
integer::dg_evb_mode,mat_size,nat,nat3
real(kind=8)::g_V12(nat6)
real(kind=8)::expo,d_p,V12
real(kind=8)::alph(dg_evb_points)
real(kind=8)::q_qts(nat6)
real(kind=8)::act_struc(nat6) !The structure for 
                                    !distance measurement
nat=natoms
nat3=3*natoms
!
!     If only energies shall be used
!

if (dg_evb_mode .eq. 1) then
   g_V12=0.d0
   do j=1,dg_evb_points
      call deltaq(act_struc,point_int(:,j),q_qts)
      d_p=dot_product(q_qts,q_qts)
      expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
      if (expo .lt. g_thres) cycle

      do k=1,nat6
         g_V12(k)=g_V12(k)-0.5d0*alph(j)*alph(j)*b_vec(j)*d_p*q_qts(k)*expo
      end do
   end do
!
!     If energies and gradients shall be used
!
else if (dg_evb_mode .eq. 2) then
   g_V12=0.d0
!
!     Calculate "global" parameters for each reference point
!
   do j=1,dg_evb_points
      call deltaq(act_struc,point_int(:,j),q_qts)
      d_p=dot_product(q_qts,q_qts)
      expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
      if (expo .lt. g_thres) cycle

      do l=1,nat6
!
!     At first, add the g(0,0)-derivatives
!

         g_V12(l)=g_V12(l)-0.5d0*alph(j)*alph(j)*b_vec((j-1)*nat6+j)*d_p*&
                 &q_qts(l)*expo
!
!     Then, add the g(i,0)-derivatives
!     distinguish two cases: q_i in prefactor is variable q_i after
!     that is differentiated or q_i is one of the other entries 
!

         do k=1,nat6
             if (k .eq. l) then
                g_V12(l)=g_V12(l)-b_vec(j+(j-1)*nat6+k)*(alph(j)*q_qts(k)*&
                     & q_qts(l)-1.d0)*expo
             else 
                g_V12(l)=g_V12(l)-b_vec(j+(j-1)*nat6+k)*alph(j)*q_qts(l)*q_qts(k)*expo
             end if
         end do
      end do
   end do

!
!     If energies, gradients and hessians shall be used
!

else if (dg_evb_mode .eq. 3) then
   g_V12=0.d0
   block=1+nat6+nat6*(nat6+1)/2
   first=(1+nat6)
!
!     Calculate "global" parameters for each reference point
!
   do j=1,dg_evb_points
      call deltaq(act_struc,point_int(:,j),q_qts)
      d_p=dot_product(q_qts,q_qts)
      expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
      if (expo .lt. g_thres) cycle

      do l=1,nat6
!
!     At first, add the g(0,0)-derivatives
!
         g_V12(l)=g_V12(l)-0.5d0*alph(j)*alph(j)*b_vec((j-1)*(block-1)+j)*d_p*&
                 & q_qts(l)*expo
!
!     Then, add the g(i,0)-derivatives
!     distinguish two cases: q_i in prefactor is variable q_i after
!     that is differentiated or q_i is one of the other entries 
!

         do k=1,nat6
             if (k .eq. l) then
                g_V12(l)=g_V12(l)-b_vec(k+(j-1)*block+1)*(alph(j)*q_qts(k)*&
                        & q_qts(l)-1.d0)*expo
             else
                g_V12(l)=g_V12(l)-b_vec(k+(j-1)*block+1)*alph(j)*q_qts(l)*q_qts(k)*expo
             end if         
         end do

!
!     At least, add the g(i,j)-energy-terms
!     distinguish five cases: both q_k indices equal (q_l equal to them or not)
!     both q_k indices different (q_l equal to first, equal to second, equal to none)
!
         do k=1,nat6
            do m=k,nat6
               inc=inc+1
               if (m .eq. k) then   
                  if (l .eq. k) then
                     g_V12(l)=g_V12(l)-0.5d0*b_vec(inc+(j-1)*block+first)*q_qts(k)*&
                             & (alph(j)*q_qts(k)*q_qts(k)-2.d0)*expo     
                  else                      
                     g_V12(l)=g_V12(l)-0.5d0*b_vec(inc+(j-1)*block+first)*alph(j)*&
                             & q_qts(k)*q_qts(k)*q_qts(l)*expo
                  end if
               else 
                  if (l .eq. k) then
                     g_V12(l)=g_V12(l)-b_vec(inc+(j-1)*block+first)*q_qts(m)*(alph(j)*&
                             & q_qts(k)*q_qts(k)-1.d0)*expo
                  else if (l .eq. m) then
                     g_V12(l)=g_V12(l)-b_vec(inc+(j-1)*block+first)*q_qts(k)*(alph(j)*&
                             & q_qts(m)*q_qts(m)-1.d0)*expo
                  else 
                     g_V12(l)=g_V12(l)-b_vec(inc+(j-1)*block+first)*alph(j)*q_qts(m)*&
                             & q_qts(k)*q_qts(l)*expo
                  end if
               end if
            end do
         end do
         inc=0
      end do
   end do
end if

return
end subroutine sum_dv12
