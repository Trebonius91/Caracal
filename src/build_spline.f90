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
!     The subroutine build_spline parameterizes the second derivatives 
!     which are needed to setup a cubic spline interpolation of a 1D
!     curve in multidimensional space. This is needed for the representation
!     of the reaction path in the case of the reaction path (RP)-EVB 
!     coupling term.
!     The routine is mostly taken from: Numerical Reciepes, 3rd ed., p. 148-149
!
!     part of EVB
!
subroutine build_spline(intdim,infdim,points,pts,ptsin,y2_nd,s,s_tot) 
implicit none 
integer::intdim,infdim   ! number of internal coordinates and number of E,G,H informations
integer::dim  ! the total spline dimension
integer::points,n  ! number of points om the path
integer::i,j,k  ! loop indices
integer::ii,im   ! loop indices
real(kind=8)::sig,p,yp1,ypn,qn,un  ! further parameters for interpolation
real(kind=8)::rad,db,de  ! for nD algorithm
real(kind=8)::rad2  ! TEST
real(kind=8)::all_rad(points) ! test for ana grad
real(kind=8)::all_rad2(points)    ! TEST
real(kind=8)::s_tot   ! the total length of the IRC for rescaling
real(kind=8)::ss  ! total arclength (for normalization)
real(kind=8)::ss2  ! TEST
real(kind=8)::fprime ! function of extrapolation for first/last points derivative
real(kind=8)::ptsin(points,intdim+infdim)  ! input array for nD data to interpolate
real(kind=8)::pts(intdim+infdim,points)   ! inverse nD data points
real(kind=8)::s(points)  ! the parameter to model the curve
real(kind=8)::s2(points) ! TEST
real(kind=8),allocatable::yv(:),xv(:)  ! arrays with s-values and coordinates 
                                       ! for the respective dimension
real(kind=8),allocatable::u(:),y2(:)  ! y2: second derivative (one dimension)
real(kind=8)::y2_nd(points,intdim+infdim)  ! all second derivatives (nD)

!
!     Define the total number of points to be used (number of structures on the 
!     reaction path)
!     Define the dimensionality of the interpolated data (internal coordinates and 
!     energies = nat6+1)
!     and fill the final data arrays
!
!
dim=intdim+infdim
n=points
allocate(yv(n),xv(n))  !for 1D calculations
allocate(u(n-1),y2(n))   ! for 1D calculations
s(1)=0.d0  ! start of the parameter range
!
!     Determine the topology of the path points and translate them into arclengths
!     of the s-values
!
open(unit=18,file="test.plot",status="unknown")
do i=1,points
   write(18,*) ptsin(i,:)
end do
close(18)

do i=1,n
!
!     different to NR: no modulo, only open curves
!
   ii=i
   im=ii-1
!
!     Calculate the euclidian distance between two reference points
!     and add it to the total parameter value
!  
   if (i .gt. 1) then
      rad=0.d0
      rad2=0.d0   ! TEST
      do j=1,intdim
         rad=rad+(ptsin(ii,j)-ptsin(im,j))*(ptsin(ii,j)-ptsin(im,j))
      !   rad=sqrt(rad)  !  square root leads to nonsymmetric paths!
      end do
      rad2=rad2+abs((ptsin(ii,intdim+2)-ptsin(im,intdim+2)))   ! TEST

!
!     Build square root of length parameter outside the inner loop!
!
    
      rad=sqrt(rad)
      all_rad(i)=rad
      all_rad2(i)=rad2       ! TEST
      s(i)=s(i-1) + rad
      s2(i)=s2(i-1) + rad2   ! TEST
!
!     Cancel with error if two consecutive points are identical
!
      if (s(i) .eq. s(i-1)) then 
         write(*,*) "Error in spline curve interpolation!"
         write(*,*) "Are there two identical structures in the structure file?"
         call fatal
      end if 
   end if
end do

!
!     rescale the parameter s to the interval [0,1]
!
ss=s(n)-s(1)
ss2=s2(n)-s2(1)   ! TEST
do i=1,n
   s(i)=s(i)/ss
   s2(i)=s2(i)/ss2
end do
!
!     Store the total length of the path for later use
!
s_tot=ss
!
!     Construct the splines: for each dimension, a pseudo 1D spline will
!     be constructed, endpoint derivatives will be used to asure convergence 
!     at the interval borders
!
do j=1,dim
   u=0.d0
   y2=0.d0
!
!     start- and endpoints have 4 point interpolation, if less than 4 points 
!     are given, set them to high default value
!
   if (n .lt. 4) then
      db=1D50
      de=1D50
   else
      db=fprime(n,s,pts(j,:),1)
      de=fprime(n,s,pts(j,:),-1)
   end if
!
!     Now do the 1D spline interpolation procedure of 1d_spline.f90 for each 
!     dimension separately!
!     s plays the role of the x coordinate array for each dimension!   
!     ---> first, calculate all second derivatives! 
!
   yp1=db
   ypn=de
   xv=s
   yv=pts(j,:)
!
!     set the lower boundary condition: either natural or specified
!       first derivative
!
   if (yp1 .gt. 0.99D50) then
      y2(1)=0.0
      u(1)=0.0
   else
      y2(1)=-0.5
      u(1)=(3.0/(xv(2)-xv(1)))*((yv(2)-yv(1))/(xv(2)-xv(1))-yp1)
   end if
!
!     decomposition loop of the tridiagonal algorithm
!
   do i=2,n-2
      sig=(xv(i)-xv(i-1))/(xv(i+1)-xv(i-1))
      p=sig*y2(i-1)+2.d0
      y2(i)=(sig-1.d0)/p
      u(i)=(yv(i+1)-yv(i))/(xv(i+1)-xv(i))-(yv(i)-yv(i-1))/(xv(i)-xv(i-1))
      u(i)=(6.d0*u(i)/(xv(i+1)-xv(i-1))-sig*u(i-1))/p
   end do
!
!     set the upper boundary condition: either natural or specified 
!       first derivative
!
   if (ypn .gt. 0.99D50) then
      qn=0.0d0
      un=0.0d0
   else
      qn=0.5
      un=(3.d0/(xv(n)-xv(n-1)))*(ypn-(yv(n-1)-yv(n-1))/(xv(n)-xv(n-1)))
   end if
!
!     backsubstitution loop of the tridiagonal algorithm
!     --> final values of the second derivatives
!
   y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
   do k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
   end do
!
!     fill the y2 values into global second derivative array
!  
   y2_nd(:,j)=y2
end do
return
end subroutine build_spline


!
!     Utility for estimating the derivatives at the endpoints
!     x and y point to the abscissa and ordinate for the endpoint.
!     If pm is +1 points to the right will be used (left endpoint),
!     if it is -1, points to the left will be used (right endpoint).
!
function fprime(n,x,y,pm)
implicit none
integer::n,pm,pt
real(kind=8)::fprime
real(kind=8)::x(n),y(n)
real(kind=8)::s1,s2,s3,s12,s23,s13
if (pm.eq.1) pt=1  ! left zone
if (pm.eq.-1) pt=n  ! right zone

s1=x(pt)-x(pt+pm*1)
s2=x(pt)-x(pt+pm*2)
s3=x(pt)-x(pt+pm*3)
s12=s1-s2
s13=s1-s3
s23=s2-s3
fprime=-(s1*s2/(s13*s23*s3))*y(pt+pm*3)+(s1*s3/(s12*s2*s23))*y(pt+pm*2)&
       & -(s2*s3/(s1*s12*s13))*y(pt+pm*1)+(1.d0/s1+1.d0/s2+1.d0/s3)*y(pt)

return
end function fprime

