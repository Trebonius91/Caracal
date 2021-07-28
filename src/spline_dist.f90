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
!     The subroutine spline_dist calculates the parameter values s
!     of a certain spline polynomial for which the distance between it 
!     and a given point in coordinate space (here: geometries) is minimal 
!
!     part of EVB
!
subroutine spline_dist(dim,s_act,inter,inter_out,act_bead)
use evb_mod
implicit none
integer::dim   ! number of internal coordinates
integer::act_bead  ! the actual bead number/index to calculate
integer::n  ! number of points om the path
integer::i,j,k  ! loop indices
integer::j1,klo,khi  ! for interpolation routine
real(kind=8)::h,a,b    ! for interpolation routine
real(kind=8)::inter(dim)  ! the internal coordinates of the actual point
real(kind=8)::distvec(dim)  ! the returned distance vector of internals
real(kind=8)::s_act  ! the actual s value that shall be interpolated
real(kind=8)::ref_dist   ! current distance to an IRC ref point
integer::i_best,i_second  ! the positions of the nearest reference points
real(kind=8)::d_best,d_second  ! the distances to the nearest and second point
real(kind=8)::d_right,d_left  ! distances to the two neighbors of the best point
real(kind=8)::s_lo,s_hi,s_mid   ! borders and middle of the actual interval in s
real(kind=8)::path_hi(dim),path_lo(dim),inter_out(dim)  ! other coordinate arrays
real(kind=8)::step_diff(dim)  ! comparison tool to previous timesteps for rpmd.x
integer::i_start,i_end    ! lower and upper boundary for nearest structure search
logical::search_global  ! indicate if the whole IRC has been tested for nearest
real(kind=8)::dist_first,dist_last  ! for scan problem: distances to first or last frame
!
!     First, determine the interval, in which the needed interpolation
!     function is defined. It will be approximated, that the interval between the 
!     nearest reference point coordinates (to the given structure) and the 
!     neighbor reference point with the lower distance to the given structure 
!     is the one in which the nearest point is located
!     In practice: if an internal structure from run before is known, check, if 
!     the change between them is small: then check only +/- 5 intervals 
!     to the one before 
!
call int_diff(inter,inter_old(act_bead,:),step_diff)
!write(299,*) dot_product(step_diff,step_diff)

!
!     If the distance to the previous structure is too big, scan the whole IRC
!     for the nearest reference structure 
!     Scan also the whole IRC if we are in the generation step of the RP-EVB,
!     i.e. if the Xi values of the RP-EVB points shall be determined
!
if ((dot_product(step_diff,step_diff) .gt. change_tol) .or. (spline_all)) then
   i_start=1
   i_end=rp_irc_points
   search_global=.true.
else
!
!     If the distance to the previous structure is small, scan only the last
!     and plus/minus 5 intervals for the nearest structure 
!
   i_start=i_best_old(act_bead)-irc_local
   if (i_start .lt. 1) then
      i_start=1
   end if

   i_end=i_best_old(act_bead)+irc_local
   if (i_end .gt. rp_irc_points) then
      i_end=rp_irc_points
   end if
   search_global=.false.
end if
!
!    start the rough search for the nearest IRC structure 
!   
do i=i_start,i_end
   ref_dist=0.d0
   call int_diff(rp_spl_ref(1:dim,i),inter,distvec)
   ref_dist=dot_product(distvec,distvec)
   if (i .eq. i_start) then
      i_best=i_start
      d_best=ref_dist
   else
      if (ref_dist .lt. d_best) then
         i_best=i
         d_best=ref_dist
      end if
   end if
end do
!
!    Edit 29.11.2018: If for a global search the nearest structure was more than 
!       glob_seach_limit^2 away from the actual one, assume that 
!       we have a structure outside the path! 
!       Then, give s=0 or s=1 directly and return to the caller routine!
!      
!
if (search_global .and. (.not. recross_calc)) then
   if (d_best .gt. path_dist_limit**2) then
      write(*,*) "NOTE: The search for the nearest IRC structure on the path gave no"
      write(*,*) " useful result! Therefore s=0 or s=1 will be assumed!"
      write(*,*) " The distance was:",sqrt(d_best),", the limit is:",path_dist_limit,"."
      write(*,*) " Add/alter the keyword PATH_DIST_LIMIT if the warning occurs too often!"
      call int_diff(rp_spl_ref(1:dim,1),inter,distvec)
      dist_first=dot_product(distvec,distvec)          
      call int_diff(rp_spl_ref(1:dim,rp_irc_points),inter,distvec)
      dist_last=dot_product(distvec,distvec)
      inter_old(act_bead,:)=inter(:)
      if (dist_first .lt. dist_last) then
         s_act=0.d0
         i_best_old=1
         return
      else if (dist_first .gt. dist_last) then
         s_act=1.d0
         i_best_old=rp_irc_points
         return
      else 
         write(*,*) "ERROR! The distances to the first and last point on the path are equal!"
         call fatal
      end if
   end if
end if
!
!
!    store procedural global arrays for next step etc
!
i_best_old(act_bead)=i_best
inter_old(act_bead,:)=inter(:)

!
!     Which neighbor of the best IRC point is nearer to the coordinate point?
!
d_left=0.d0
d_right=0.d0
do j=1,dim
   d_left=d_left+(rp_spl_ref(j,i_best-1)-inter(j))**2
   d_right=d_right+(rp_spl_ref(j,i_best+1)-inter(j))**2
end do
!
!     Define the needed fixpoints for the interval (khi,klo)
!
if (d_left .le. d_right) then
   i_second=i_best-1
   klo=i_second
   khi=i_best
else 
   i_second=i_best+1
   klo=i_best
   khi=i_second
end if

!
!     Now solve the equation of the extremum of the distance between 
!     the curve and the given coordinate point
!     --> This cannot be done analytically, because the formula changes 
!     with each additional dimension and the scaling of complexity seems 
!     to be exponential!
!     
!     Strategy: Half the interval in each step, determine, which end is 
!     nearer to the point and set the new interval to this end to the 
!     middle of the old interval and so on...
!
s_lo=rp_spl_s(klo)
s_hi=rp_spl_s(khi)
h=rp_spl_s(khi)-rp_spl_s(klo)
do 
!
!     Calculate coordinates for the actual borders of the interval
!
   s_act=s_lo
   a=(rp_spl_s(khi)-s_act)/h
   b=(s_act-rp_spl_s(klo))/h 
   do j=1,dim
      path_lo(j)=a*rp_spl_ref(j,klo)+b*rp_spl_ref(j,khi)+((a*a*a-a)*&
           & rp_spl_d2(klo,j)+(b*b*b-b)*rp_spl_d2(khi,j))*(h*h)/6.d0
   end do
 
   s_act=s_hi
   a=(rp_spl_s(khi)-s_act)/h
   b=(s_act-rp_spl_s(klo))/h
   do j=1,dim
      path_hi(j)=a*rp_spl_ref(j,klo)+b*rp_spl_ref(j,khi)+((a*a*a-a)*&
           & rp_spl_d2(klo,j)+(b*b*b-b)*rp_spl_d2(khi,j))*(h*h)/6.d0
   end do
   d_left=0.d0
   d_right=0.d0
!
!     Calculate the distances of the interval borders to the actual structure
!
   call int_diff(path_hi,inter,distvec)
   d_left=dot_product(distvec,distvec)
   call int_diff(path_lo,inter,distvec)
   d_right=dot_product(distvec,distvec)
!
!     If interval is smaller than given condition, leave cycle
!
   if (s_hi-s_lo .lt. interp_tol) exit
!
!     Divide the interval in the middle and decide, which part of it shall 
!     remain
!
!   write(*,*) "d_left,d_right",d_left,d_right
   if (d_left .ge. d_right) then
      s_hi=s_lo+0.5d0*(s_hi-s_lo)
   else 
      s_lo=s_hi-0.5d0*(s_hi-s_lo)
   end if
!   write(*,*) s_hi,s_lo
end do
!
!     If cycle is converged, determine the middle of the last interval as 
!     the final s value
!
s_act=s_lo+0.5d0*(s_hi-s_lo)

a=(rp_spl_s(khi)-s_act)/h
b=(s_act-rp_spl_s(klo))/h
do j=1,dim
   inter_out(j)=a*rp_spl_ref(j,klo)+b*rp_spl_ref(j,khi)+((a*a*a-a)*&
        & rp_spl_d2(klo,j)+(b*b*b-b)*rp_spl_d2(khi,j))*(h*h)/6.d0
end do
!stop "splinedist!"

return
end subroutine spline_dist
