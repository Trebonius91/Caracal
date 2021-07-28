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
!     "build_dmat" calculates all elements of the coefficient matrix 
!     for the linar system of equations whose solutions are the needed 
!     prefactors of the DG-EVB coupling terms
!
!     part of EVB
!

subroutine build_dmat(dg_evb_mode,mat_size,alph)
use evb_mod
implicit none
integer::i,j,l,k,mp,m,o,r,p,zzz
integer::indii,indj,inc,lin
integer::col,row,block
integer::dg_evb_mode,mat_size
integer::aa ! short script for add_alph
real(kind=8)::expo_aa ! separate exponent for add_alph
real(kind=8)::starttime,endtime ! test: time_measurement
real(kind=8)::d_p,expo,det
real(kind=8)::alph(dg_evb_points)
real(kind=8)::q_qts(nat6)
real(kind=8)::step,hi,lo,val,step2  ! numerical calculation
real(kind=8)::difftmp(nat6,nat6,2)
real(kind=8)::numdiff(nat6,nat6)
real(kind=8)::all_int_store(nat6,dg_evb_points)
real(kind=8)::vals(4)
real(kind=8)::d_mat_num(mat_size,mat_size)
mp=dg_evb_points
block=1+nat6+nat6*(nat6+1)/2  ! one horizonal dg_evb-block for mode 3
inc=dg_evb_points*(1+nat6)   ! size of g and dg blocks for mode 3

!
!    VARIANT A: NUMERICAL CALCULATION  
!
if (dg_mat_num) then
!   call cpu_time(starttime)
!   do zzz=1,1000

   d_mat=0.d0
   step=0.000001d0   ! size of numerical differentiation step
   step2=0.000001d0  ! size of inner numericla differentiation step
!
!    store the all_int array to avoid slight decreasing of values due to
!    numerous applications of numerical derivatives
! 
   all_int_store=all_int
   
   if (dg_evb_mode .eq. 1) then
      do i=1,dg_evb_points
         do j=1,dg_evb_points
         call deltaq(all_int(:,i),point_int(:,j),q_qts)
         d_p=dot_product(q_qts,q_qts)
         expo=exp(-0.5d0*alph(j)*d_p)
         d_mat(i,j)=(1+0.5d0*alph(j)*d_p)*expo
      end do
   end do
   else if (dg_evb_mode .eq. 2) then
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            call deltaq(all_int(:,i),point_int(:,j),q_qts)
            d_p=dot_product(q_qts,q_qts)
         expo=exp(-0.5d0*alph(j)*d_p)
   
            do l=0,nat6
               if (l .eq. 0) then
!
!    The normal g(0,0) terms     
! 
                  d_mat(i,(j-1)*(nat6+1)+1)=(1+0.5d0*alph(j)*d_p)*expo
               else
!
!    The normal g(i,0) terms
!
               d_mat(i,(j-1)*(nat6+1)+1+l)=q_qts(l)*expo
            end if
         end do
      end do
   end do
!
!   Then, fill the lower block with the derivatives
!
      do i=1,dg_evb_points
         do j=1,dg_evb_points
         do l=0,nat6
               do k=1,nat6
                  if (l .eq. 0) then
!
!    The d(g(0,0))/dq terms
!                
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step 
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2d0*step
                        end if
                        call deltaq(all_int(:,i),point_int(:,j),q_qts)
                        d_p=dot_product(q_qts,q_qts)
                        expo=exp(-0.5d0*alph(j)*d_p)
                        val=(1+0.5d0*alph(j)*d_p)*expo
                        if (m .eq. 1) then
                           hi=val
                        else if (m .eq. 2) then
                           lo=val
                        end if
                     end do   
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*(nat6+1)+1)=(hi-lo)/(2*step)
                     all_int(k,i)=all_int(k,i)+step          
                  else
!
!    The d(g(i,0))/dq terms
!
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2d0*step
                        end if
                        call deltaq(all_int(:,i),point_int(:,j),q_qts)
                        d_p=dot_product(q_qts,q_qts)
                        expo=exp(-0.5d0*alph(j)*d_p)
                        val=q_qts(l)*expo
                        if (m .eq. 1) then
                        hi=val
                        else if (m .eq. 2) then
                           lo=val
                        end if
                     end do
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*(nat6+1)+l+1)=(hi-lo)/(2*step)
                     all_int(k,i)=all_int(k,i)+step
                  end if
               end do
            end do
         end do
      end do
   else if (dg_evb_mode .eq. 3) then
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            do l=0,nat6+1
               call deltaq(all_int(:,i),point_int(:,j),q_qts)
               d_p=dot_product(q_qts,q_qts)
               if (l .eq. 0) then
!
!    The normal g(0,0) terms     
!
                  expo=exp(-0.5d0*alph(j)*d_p)
                  d_mat(i,(j-1)*block+1)=(1+0.5d0*alph(j)*d_p)*expo
               else if (l .le. nat6) then
!
!    The normal g(i,0) terms
!
                  expo=exp(-0.5d0*alph(j)*d_p)
                  d_mat(i,(j-1)*block+1+l)=q_qts(l)*expo
               else 
!
!    The normal g(i,j) terms
!               
                  expo=exp(-0.5d0*alph(j+add_alph)*d_p)
                  col=0   ! The actual column index to fill
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
                        col=col+1
                           d_mat(i,(j-1)*block+l+col)=0.5d0*q_qts(indii)*q_qts(indj)*expo
                        else 
                           col=col+1
                           d_mat(i,(j-1)*block+l+col)=q_qts(indii)*q_qts(indj)*expo
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
!
!   Then, fill the lower block with the first derivatives
!
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            do l=0,nat6+1
               if (l .eq. 0) then
!
!    The d(g(0,0))/dq terms
!                 
                  do k=1,nat6 
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2d0*step
                        end if
                        call deltaq(all_int(:,i),point_int(:,j),q_qts)
                        d_p=dot_product(q_qts,q_qts)
                        expo=exp(-0.5d0*alph(j)*d_p)
                        val=(1+0.5d0*alph(j)*d_p)*expo
                        if (m .eq. 1) then
                        hi=val
                        else if (m .eq. 2) then
                           lo=val
                        end if
                     end do
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*block+1)=(hi-lo)/(2*step)
               end do
               else if (l .le. nat6) then
!
!    The d(g(i,0))/dq terms
!
                  do k=1,nat6
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2d0*step
                        end if
                        call deltaq(all_int(:,i),point_int(:,j),q_qts)
                        d_p=dot_product(q_qts,q_qts)
                        expo=exp(-0.5d0*alph(j)*d_p)
                        val=q_qts(l)*expo
                        if (m .eq. 1) then
                           hi=val
                        else if (m .eq. 2) then
                           lo=val
                        end if
                     end do
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+1)=(hi-lo)/(2*step)
                     all_int(k,i)=all_int(k,i)+step
                  end do
               else 
!
!    The d(g(i,j))/dq terms
!    Make a distinction between on and off diagonal with kronecker delta
!              
                  col=0
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
                           col=col+1
                           do k=1,nat6
                              do m=1,2
                                 if (m .eq. 1) then
                                    all_int(k,i)=all_int(k,i)+step
                                 else if (m .eq. 2) then
                                    all_int(k,i)=all_int(k,i)-2d0*step
                                 end if
                                 call deltaq(all_int(:,i),point_int(:,j),q_qts)
                                 d_p=dot_product(q_qts,q_qts)
                                 expo=exp(-0.5d0*alph(j+add_alph)*d_p)
                                 val=0.5d0*q_qts(indii)*q_qts(indj)*expo
                                 if (m .eq. 1) then
                                    hi=val
                                 else if (m .eq. 2) then
                                    lo=val
                                 end if
                              end do
                              d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                   (hi-lo)/(2*step)
                              all_int(k,i)=all_int(k,i)+step
                           end do
                        else 
                           col=col+1 
                           do k=1,nat6
                              do m=1,2
                                 if (m .eq. 1) then
                                    all_int(k,i)=all_int(k,i)+step
                                 else if (m .eq. 2) then
                                    all_int(k,i)=all_int(k,i)-2d0*step
                                 end if
                                 call deltaq(all_int(:,i),point_int(:,j),q_qts)
                                 d_p=dot_product(q_qts,q_qts)
                                 expo=exp(-0.5d0*alph(j+add_alph)*d_p)
                                 val=q_qts(indii)*q_qts(indj)*expo
                                 if (m .eq. 1) then
                                    hi=val
                                 else if (m .eq. 2) then
                                    lo=val
                                 end if
                              end do
                              d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                    (hi-lo)/(2*step)
                              all_int(k,i)=all_int(k,i)+step
                           end do
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
!
!   Last but not least, fill the lowest block with second derivatives
!
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            do l=0,nat6+1
               if (l .eq. 0) then
!
!    The d^2(g(0,0))/dq^2 terms
! 
                  do k=1,nat6  
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2.d0*step
                        end if
                        do p=1,nat6
                           do r=1,2
                              if (r .eq. 1) then
                                 all_int(p,i)=all_int(p,i)+step2
                              else if (r .eq. 2) then
                                 all_int(p,i)=all_int(p,i)-2.d0*step2
                              end if 
                              call deltaq(all_int(:,i),point_int(:,j),q_qts)
                              d_p=dot_product(q_qts,q_qts)
                              expo=exp(-0.5d0*alph(j)*d_p)
                              val=(1+0.5d0*alph(j)*d_p)*expo
                              if (r .eq. 1) then
                                 hi=val
                              else if (r .eq. 2) then
                                 lo=val
                              end if
                           end do
                           all_int(p,i)=all_int(p,i)+step2
                           difftmp(k,p,m)=(hi-lo)/(2.d0*step2)
                        end do
                      end do
                      all_int(k,i)=all_int(k,i)+step
                  end do
                  do m=1,nat6
                     do k=m,nat6
                        numdiff(m,k)=(difftmp(m,k,1)-difftmp(m,k,2))/(2.d0*step)
                     end do
                  end do
                  lin=0  ! the actual line index to fill
                  do m=1,nat6
                     do p=m,nat6
                        lin=lin+1 
                        d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1)=&
                                      &numdiff(m,p)
                     end do
                  end do
               else if (l .le. nat6) then
!
!    The d^2(g(i,0))/dq^2 terms
!  
                  do k=1,nat6
                     do m=1,2
                        if (m .eq. 1) then
                           all_int(k,i)=all_int(k,i)+step
                        else if (m .eq. 2) then
                           all_int(k,i)=all_int(k,i)-2.d0*step
                        end if
                        do p=1,nat6
                           do r=1,2
                              if (r .eq. 1) then
                                 all_int(p,i)=all_int(p,i)+step2
                              else if (r .eq. 2) then
                                 all_int(p,i)=all_int(p,i)-2.d0*step2
                              end if
                              call deltaq(all_int(:,i),point_int(:,j),q_qts)
                              d_p=dot_product(q_qts,q_qts)
                              expo=exp(-0.5d0*alph(j)*d_p)
                              val=q_qts(l)*expo
                              if (r .eq. 1) then
                                 hi=val
                           else if (r .eq. 2) then
                                 lo=val
                              end if
                           end do
                           all_int(p,i)=all_int(p,i)+step2
                           difftmp(k,p,m)=(hi-lo)/(2.d0*step2)
                        end do
                      end do
                      all_int(k,i)=all_int(k,i)+step
                  end do
                  do m=1,nat6
                     do k=m,nat6
                        numdiff(m,k)=(difftmp(m,k,1)-difftmp(m,k,2))/(2.d0*step)
                     end do
                  end do
                  lin=0  ! the actual line index to fill
                  do m=1,nat6
                     do p=m,nat6
                        lin=lin+1
                        d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                                     &numdiff(m,p)
                     end do
                  end do
               else 
!
!    The d^2(g(i,j))/dq^2 terms
!    Make a distinction between on and off diagonal with kronecker delta
!  
                  col=0
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
                           col=col+1
                           do k=1,nat6
                              do m=1,2
                                 if (m .eq. 1) then
                                    all_int(k,i)=all_int(k,i)+step
                                 else if (m .eq. 2) then
                                    all_int(k,i)=all_int(k,i)-2d0*step
                                 end if
                                 do p=1,nat6
                                    do r=1,2
                                       if (r .eq. 1) then
                                          all_int(p,i)=all_int(p,i)+step2
                                       else if (r .eq. 2) then
                                          all_int(p,i)=all_int(p,i)-2.d0*step2
                                       end if
                                       call deltaq(all_int(:,i),point_int(:,j),q_qts)
                                       d_p=dot_product(q_qts,q_qts)
                                       expo=exp(-0.5d0*alph(j+add_alph)*d_p)
                                       val=0.5d0*q_qts(indii)*q_qts(indj)*expo
                                       if (r .eq. 1) then
                                          hi=val
                                       else if (r .eq. 2) then
                                          lo=val
                                       end if
                                    end do
                                    all_int(p,i)=all_int(p,i)+step2
                                    difftmp(k,p,m)=(hi-lo)/(2.d0*step2)
                                 end do
                              end do
                              all_int(k,i)=all_int(k,i)+step
                           end do
                           do m=1,nat6
                              do k=m,nat6
                                 numdiff(m,k)=(difftmp(m,k,1)-difftmp(m,k,2))/(2.d0*step)
                              end do
                           end do
                           lin=0  ! the actual line index to fill
                           do m=1,nat6
                              do p=m,nat6
                                 lin=lin+1
                                 d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                       &numdiff(m,p)
                              end do
                           end do
                        else 
                           col=col+1
                           do k=1,nat6
                              do m=1,2
                                 if (m .eq. 1) then
                                    all_int(k,i)=all_int(k,i)+step
                                 else if (m .eq. 2) then
                                    all_int(k,i)=all_int(k,i)-2d0*step
                                 end if
                                 do p=1,nat6
                                    do r=1,2
                                       if (r .eq. 1) then
                                          all_int(p,i)=all_int(p,i)+step2
                                       else if (r .eq. 2) then
                                          all_int(p,i)=all_int(p,i)-2.d0*step2
                                       end if
                                       call deltaq(all_int(:,i),point_int(:,j),q_qts)
                                       d_p=dot_product(q_qts,q_qts)
                                       expo=exp(-0.5d0*alph(j+add_alph)*d_p)
                                       val=q_qts(indii)*q_qts(indj)*expo
                                       if (r .eq. 1) then
                                          hi=val
                                       else if (r .eq. 2) then
                                          lo=val
                                       end if
                                    end do
                                    all_int(p,i)=all_int(p,i)+step2
                                    difftmp(k,p,m)=(hi-lo)/(2.d0*step2)
                                 end do
                              end do
                              all_int(k,i)=all_int(k,i)+step
                           end do
                           do m=1,nat6
                              do k=m,nat6
                                 numdiff(m,k)=(difftmp(m,k,1)-difftmp(m,k,2))/(2.d0*step)
                              end do
                           end do
                           lin=0  ! the actual line index to fill
                           do m=1,nat6
                              do p=m,nat6
                                 lin=lin+1
                                 d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                        &numdiff(m,p)
                              end do
                           end do
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
   end if
!
!     restore the all_int array
!
   all_int=all_int_store
!   end do
!   call cpu_time(endtime)
!   write(*,*) "1000 cycles needed time of (s):",endtime-starttime
!   stop "Hiihpih"

   d_mat_num=d_mat  ! test
!
!     VARIANT B: ANALYTICAL CALCULATION  
!   
else   
   
!   call cpu_time(starttime)
!   do zzz=1,1000 
   d_mat=0.d0 
!
!  ---- MODE 1: only energies ----
!
   if (dg_evb_mode .eq. 1) then
      do i=1,dg_evb_points
         do j=1,dg_evb_points
!
!     All gaussians: g(0,0) : (1+0.5*alph*dq^2)*exp(-0.5*alph*dq^2)
!
            call deltaq(all_int(:,i),point_int(:,j),q_qts)
            d_p=dot_product(q_qts,q_qts)
            expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
            if (expo .lt. g_thres) cycle

            d_mat(i,j)=(1+0.5d0*alph(j)*d_p)*expo
         end do
      end do
!
!  ---- MODE 2: energies and gradients ----
!
   else if (dg_evb_mode .eq. 2) then
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            call deltaq(all_int(:,i),point_int(:,j),q_qts)
            d_p=dot_product(q_qts,q_qts)
            expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
            if (expo .lt. g_thres) cycle

!
!     Functions itself 
!
            do l=0,nat6
               if (l .eq. 0) then
!
!     The g(0,0) terms :  (1+0.5*alph*dq^2)*exp(-0.5*alph*dq^2)
!   
                  d_mat(i,(j-1)*(nat6+1)+1)=(1+0.5d0*alph(j)*d_p)*expo
               else
!
!     The g(i,0) terms : dq_l*exp(-0.5*alph*dq^2)
!
                  d_mat(i,(j-1)*(nat6+1)+1+l)=q_qts(l)*expo
               end if
            end do
!
!     First derivatives (k'th coordinate)
!
            do l=0,nat6
               do k=1,nat6
                  if (l .eq. 0) then
!
!     The d(g(0,0))/dq terms : -0.5*alph^2*q^2*dq_k*exp(-0.5*alph*dq^2)
!                
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*(nat6+1)+1)=-0.5d0*alph(j)*alph(j)*d_p*q_qts(k)*expo
                  else
!
!     The d(g(i,0))/dq terms : more cases :
!                 
                     if (k .eq. l) then
!
!     - l=k (in prefactor) : -(alph*dq_i*dq_k-1)*exp(-0.5*alph*dq^2)
!
                        d_mat((i-1)*(nat6)+k+mp,(j-1)*(nat6+1)+l+1)=-(alph(j)*q_qts(k)*q_qts(l)-1.d0)*expo
                     else
!
!     - l/=k (in prefactor) : -alph*dq_i*dq_k*exp(-0.5*alph*dq^2)
!
                        d_mat((i-1)*(nat6)+k+mp,(j-1)*(nat6+1)+l+1)=-alph(j)*q_qts(k)*q_qts(l)*expo
                     end if
                  end if
               end do
            end do
         end do
      end do
!
!  ---- MODE 3: energies, gradients and hessians ----
!
   else if (dg_evb_mode .eq. 3) then
!
!     if the DOUBLE-ALPHA option is activated, take extra alphas for g(i,j) functions (hessians)
!     Store the rather long variable add_alph into a short name: aa
!
      aa=add_alph
      do i=1,dg_evb_points
         do j=1,dg_evb_points
            call deltaq(all_int(:,i),point_int(:,j),q_qts)
            d_p=dot_product(q_qts,q_qts)
            expo=exp(-0.5d0*alph(j)*d_p)
!
!     Go to next loop if exponent is too low (neglect calculations)
!
            if (expo .lt. g_thres) cycle

            if (double_alpha) then
               expo_aa=exp(-0.5d0*alph(j+aa)*d_p)
            else 
               expo_aa=expo
            end if
!
!     Functions itself 
!
            do l=0,nat6+1
               if (l .eq. 0) then
!
!     The g(0,0) terms :  (1+0.5*alph*dq^2)*exp(-0.5*alph*dq^2)  
!
                  d_mat(i,(j-1)*block+1)=(1+0.5d0*alph(j)*d_p)*expo
               else if (l .le. nat6) then
!
!     The g(i,0) terms : dq_l*exp(-0.5*alph*dq^2)
!
                  d_mat(i,(j-1)*block+1+l)=q_qts(l)*expo
               else 
!
!     The g(i,j) terms : two major cases : 
!               
                  col=0   ! The actual column index to fill
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
!
!     - i=j (diagonal) : 0.5*dq_i*dq_j*exp(-0.5*alph*dq^2)
!
                           col=col+1
                           d_mat(i,(j-1)*block+l+col)=0.5d0*q_qts(indii)*q_qts(indj)*expo_aa
                        else 
!
!     - i/=j (off-diagonal) : dq_i*dq_j*exp(-0.5*alph*dq^2)
!
                           col=col+1
                           d_mat(i,(j-1)*block+l+col)=q_qts(indii)*q_qts(indj)*expo_aa
                        end if
                     end do
                  end do
               end if
            end do
!
!     First derivatives (k'th coordinate)
!
            do l=0,nat6+1
               if (l .eq. 0) then
!
!     The d(g(0,0))/dq terms : -0.5*alph^2*q^2*dq_k*exp(-0.5*alph*dq^2)
!                 
                  do k=1,nat6 
                     d_mat((i-1)*(nat6)+k+mp,(j-1)*block+1)=-0.5d0*alph(j)*alph(j)*d_p*q_qts(k)*expo
                  end do
               else if (l .le. nat6) then
!
!     The d(g(i,0))/dq terms : more cases : 
!
                  do k=1,nat6
                     if (k.eq.l) then
!
!     - i=k (in prefactor) : -(alph*dq_i*dq_k-1)*exp(-0.5*alph*dq^2)
!
                        d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+1)=-(alph(j)*q_qts(k)*q_qts(l)-1.d0)*expo
                     else 
!
!     - i/=k (in prefactor) : -alph*dq_i*dq_k*exp(-0.5*alph*dq^2)
!
                        d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+1)=-alph(j)*q_qts(k)*q_qts(l)*expo
                     end if
                  end do
               else 
!
!      The d(g(i,j))/dq terms : two major cases : 
!              
                  col=0
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
! 
!      - i=j (diagonal) : two subcases
!
                           col=col+1
                           do k=1,nat6
                              if (k .eq. indii) then
!
!      - - i=j=k : -0.5*dq_i*(alph*dq_i*dq_i-2)*exp(-0.5*alph*dq^2) 
!
                                 d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                      & -0.5d0*q_qts(indii)*(alph(j+aa)*q_qts(indii)*q_qts(indii)-2.d0)*expo_aa
                              else 
!
!      - - i=j/=k : -0.5*alph*dq_i*dq_j*dq_k*exp(-0.5*alph*dq^2)
!
                                 d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                      & -0.5d0*alph(j+aa)*q_qts(indii)*q_qts(indii)*q_qts(k)*expo_aa
                              end if
                           end do
                        else 
!
!      - i/=j (off-diagonal) : three subcases
!
                           col=col+1 
                           do k=1,nat6
                              if (k .eq. indii) then
!
!      - - i=k, j/=k : -dq_j*(alph*dq_i*dq_i-1)*exp(-0.5*alph*dq^2)
!
                                 d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                      & -q_qts(indj)*(alph(j+aa)*q_qts(indii)*q_qts(indii)-1.d0)*expo_aa
                              else if (k .eq. indj) then
!
!      - - i/=k, j=k : -dq_i*(alph*dq_j*dq_j-1)*exp(-0.5*alph*dq^2) 
!
                                 d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                      & -q_qts(indii)*(alph(j+aa)*q_qts(indj)*q_qts(indj)-1.d0)*expo_aa
                              else 
!
!      - - i/=k, j/=k : -alph*dq_i*dq_j*dq_k*exp(-0.5*alph*dq^2) 
!
                                 d_mat((i-1)*(nat6)+k+mp,(j-1)*block+l+col)=&
                                      & -alph(j+aa)*q_qts(indii)*q_qts(indj)*q_qts(k)*expo_aa
                              end if
                           end do
                        end if
                     end do
                  end do
               end if
            end do
!
!     Second derivatives (k'th and p'th coordinate)
!

            do l=0,nat6+1
               if (l .eq. 0) then
!
!     The d^2(g(0,0))/dq^2 terms : two cases :
! 
                  lin=0
                  do k=1,nat6  
                     do p=k,nat6
                     lin=lin+1
                        if (k .eq. p) then
!
!     - k=p (two equal partials) : 0.5*alph^2*dq^2*(alph*dq_k*dq_k-3)*exp(-0.5*alph*dq^2)
!
                           d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1)=&
                               & 0.5d0*alph(j)*alph(j)*q_qts(k)*q_qts(k)*(alph(j)*d_p-3.d0)*expo
                        else 
!
!     - k/=p (two diff. partials) : 0.5*alph^2*dq_k*dq_p*(alph*dq^2-2)*exp(-0.5*alph*dq^2)
!
                           d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1)=&
                               & 0.5d0*alph(j)*alph(j)*q_qts(k)*q_qts(p)*(alph(j)*d_p-2.d0)*expo
                        end if 
                     end do
                  end do
               else if (l .le. nat6) then
!
!     The d^2(g(i,0))/dq^2 terms : four cases : 
!  
                  lin=0
                  do k=1,nat6
                     do p=k,nat6
                        lin=lin+1
!
!     - k/=i,p/=i : two subcases
!
                        if ((k .ne. l) .and. (p .ne. l)) then
                           if (k .eq. p) then
!
!     - - k=p : alph*dq_i*(alph*dq_k*dq_k-1)*exp(-0.5*alph*dq^2)    
!
                              d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                                   & alph(j)*q_qts(l)*(alph(j)*q_qts(k)*q_qts(k)-1.d0)*expo
                           else 
!
!     - - k/=p : alph^2*dq_k*dq_p*dq_i*exp(-0.5*alph*dq^2)
!
                              d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                                   & alph(j)*alph(j)*q_qts(k)*q_qts(p)*q_qts(l)*expo
                           end if
!
!     - k/=i,p=i : alph*dq_p*(alph*dq_i*dq_i-1)*exp(-0.5*alph*dq^2)
!
                        else if ((k .ne. l) .and. (p .eq. l)) then
                           d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                               & alph(j)*(alph(j)*q_qts(l)*q_qts(l)-1.d0)*q_qts(k)*expo
!
!     - k=i,p/=i : alph*dq_k*(alph*dq_i*dq_i-1)*exp(-0.5*alph*dq^2)
!
                        else if ((k .eq. l) .and. (p .ne. l)) then
                           d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                               & alph(j)*(alph(j)*q_qts(l)*q_qts(l)-1.d0)*q_qts(p)*expo
!
!     - k=i,p=i : alph*dq_i*(alph*dq_i*dq_i-3)*exp(-0.5*alph*dq^2) 
!
                        else if ((k .eq. l) .and. (p .eq. l)) then
                           d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+1+l)=&
                               & alph(j)*q_qts(l)*(alph(j)*q_qts(l)*q_qts(l)-3.d0)*expo
                        end if
                     end do
                  end do
               else 
!
!     The d^2(g(i,j))/dq^2 terms : two major cases :
!  
                  col=0
                  do indii=1,nat6
                     do indj=indii,nat6
                        if (indii .eq. indj) then
!
!     - i=j (diagonal) : four subcases
!
                           col=col+1
                           lin=0  ! the actual line index to fill
                           do k=1,nat6
                              do p=k,nat6
                                 lin=lin+1
!
!      - - i/=k,j/=k,i/=p,j/=p : two subsubcases
!
                                 if ((k .ne. indii) .and. (p .ne. indii)) then
!
!      - - - k=p : 0.5*alph*dq_i*dq_j*(alph*dq_k*dq_k-1)*exp(-0.5*alph*dq^2)
!
                                    if (k .eq. p) then
                                        d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                            & 0.5d0*alph(j+aa)*q_qts(indii)*q_qts(indii)*(alph(j+aa)*&
                                            & q_qts(k)*q_qts(k)-1.d0)*expo_aa
!
!      - - - k/=p : 0.5*alph^2*dq_i*dq_j*dq_k*dq_p*exp(-0.5*alph*dq^2)
!
                                    else 
                                       d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                            & 0.5d0*alph(j+aa)*alph(j+aa)*q_qts(indii)*q_qts(indj)*&
                                            & q_qts(k)*q_qts(p)*expo_aa
                                    end if
!
!      - - i=j=k,p/=k : 0.5*alph*dq_p*dq_i*(alph*dq_i*dq_i-2)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indii) .and. (p .ne.indii)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                         & 0.5d0*alph(j+aa)*q_qts(p)*q_qts(indii)*(alph(j+aa)*&
                                         & q_qts(indii)*q_qts(indii)-2.d0)*expo_aa
!
!      - - i=j=p,p/=k : 0.5*alph*dq_k*dq_i*(alph*dq_i*dq_i-2)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .ne. indii) .and. (p .eq.indii)) then     
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                         & 0.5d0*alph(j+aa)*q_qts(k)*q_qts(indii)*(alph(j+aa)*&
                                         & q_qts(indii)*q_qts(indii)-2.d0)*expo_aa
!
!      - - i=j=k=p : 0.5*(alph^2*dq_i^4-5*alph*dq_i^2+2)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indii) .and. (p .eq. indii)) then 
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin,(j-1)*block+l+col)=&
                                         & 0.5d0*(alph(j+aa)*alph(j+aa)*q_qts(indii)*q_qts(indii)*&
                                         & q_qts(indii)*q_qts(indii)-5.d0*alph(j+aa)*q_qts(indii)*&
                                         & q_qts(indii)+2.d0)*expo_aa
                                 end if
                              end do
                           end do 
!
!      - i/=j (off-diagonal) : eight subcases
!
                        else 
                           col=col+1
                           lin=0
                           do k=1,nat6
                              do p=k,nat6
                                 lin=lin+1 
!
!      - - i/=k,j/=k,i/=p,j/=p : two subsubcases
!
                                 if ((k .ne. indii) .and. (k .ne. indj) .and. (p .ne. indii) &
                                       & .and. (p .ne. indj)) then
!
!      - - - k=p : alph*dq_i*dq_j*(alph*dq_k*dq_k-1)*exp(-0.5*alph*dq^2)
!
                                    if (k .eq. p) then
                                       d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                          & alph(j+aa)*q_qts(indii)*q_qts(indj)*(alph(j+aa)*q_qts(k)*&
                                          & q_qts(p)-1.d0)*expo_aa
!
!      - - - k/=p : alph^2*dq_i*dq_j*dq_k*dq_p*exp(-0.5*alph*dq^2) 
!
                                    else 
                                       d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                          & alph(j+aa)*alph(j+aa)*q_qts(i)*q_qts(j)*q_qts(indii)*&
                                          & q_qts(indj)*expo_aa 
                                    end if
!
!      - - k=i,p/=i,p/=j : alph*dq_j*(alph*dq_i*dq_i-1)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indii) .and. (p .ne. indii) .and. (p .ne. indj)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & alph(j+aa)*q_qts(indj)*(alph(j+aa)*q_qts(indii)*q_qts(indii)-1.d0)*&
                                       & q_qts(p)*expo_aa 
!
!      - - k=j,p/=i,p/=j : alph*dq_i*(alph*dq_j*dq_j-1)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indj) .and. (p .ne. indii) .and. (p .ne. indj)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & alph(j+aa)*q_qts(indii)*(alph(j+aa)*q_qts(indj)*q_qts(indj)-1.d0)*&
                                       & q_qts(p)*expo_aa
!
!      - - p=i,k/=i,k/=j : alph*dq_j*dq_k*(alph*dq_i*dq_i-1)*exp(-0.5*alph*dq^2)
!
                                 else if ((p .eq. indii) .and. (k .ne. indii) .and. (k .ne. indj)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & (alph(j+aa)*q_qts(p)*q_qts(p)-1.d0)*alph(j+aa)*q_qts(indj)*&
                                       & q_qts(k)*expo_aa
!
!      - - p=j,k/=i,k/=j : alph*dq_i*dq_k*(alph*dq_j*dq_j-1)*exp(-0.5*alph*dq^2)
!
                                 else if ((p .eq. indj) .and. (k .ne. indii) .and. (k .ne. indj)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & (alph(j+aa)*q_qts(p)*q_qts(p)-1.d0)*alph(j+aa)*q_qts(indii)*&
                                       & q_qts(k)*expo_aa
!
!      - - k=i,p=j or k=j,p=i : (alph*dq_i*dq_i-1)*(alph*dq_j*dq_j-1)*exp(-0.5*alph*dq^2)
!
                                 else if (((k .eq. indii) .and. (p .eq. indj)) .or. ((k .eq. indj) & 
                                           & .and. (p .eq. indii))) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & (alph(j+aa)*q_qts(k)*q_qts(k)-1.d0)*(alph(j+aa)*q_qts(p)*&
                                       & q_qts(p)-1.d0)*expo_aa
!
!      - - k=i,p=i : alph*dq_j*dq_i*(alph*dq_i*dq_i-3)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indii) .and. (p .eq. indii)) then
                                    d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & alph(j+aa)*q_qts(indj)*q_qts(indii)*(alph(j+aa)*q_qts(indii)*&
                                       & q_qts(indii)-3.d0)*expo_aa
!
!      - - k=j,p=j : alph*dq_i*dq_j*(alph*dq_j*dq_j-3)*exp(-0.5*alph*dq^2)
!
                                 else if ((k .eq. indj) .and. (p .eq. indj)) then
                                       d_mat(inc+(i-1)*nat6*(nat6+1)/2+lin ,(j-1)*block+l+col)=&
                                       & alph(j+aa)*q_qts(indii)*q_qts(indj)*(alph(j+aa)*q_qts(indj)*&
                                       & q_qts(indj)-3.d0)*expo_aa
                                 end if
                              end do
                           end do
                        end if
                     end do
                  end do
               end if
            end do
         end do
      end do
   end if
!   end do
!   call cpu_time(endtime)
!   write(*,*) "1000 cycles needed time of (s):",endtime-starttime
!   stop "Hiihpih"
end if
   
   
!   write(*,*) "numerical , analytical matrix:"
!   do i=1,mat_size
!      do j=1,mat_size
!      write(*,*) "i,j",i,j,d_mat_num(i,j),d_mat(i,j)
!   end do
!end do

!stop "in build_dmat!"
return
end subroutine build_dmat


