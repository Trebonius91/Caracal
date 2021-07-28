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
!     The subroutine build_spline2 parameterizes the second derivatives 
!     which are needed to setup a cubic spline interpolation of a 1D
!     curve in multidimensional space. 
!     In contrast to build_spline, the s values of the curve are read 
!     in and gradients and frequencies will be interpolated there
!     The routine is mostly taken from: Numerical Reciepes, 3rd ed., p. 148-149
!
!     part of EVB
!
subroutine build_spline2(infdim,points,pts,ptsin,y2_nd,s) 
implicit none 
integer::intdim,infdim   ! number of internal coordinates and number of E,G,H informations
integer::dim  ! the total spline dimension
integer::points,n  ! number of points om the path
integer::i,j,k  ! loop indices
integer::ii,im   ! loop indices
real(kind=8)::sig,p,yp1,ypn,qn,un  ! further parameters for interpolation
real(kind=8)::rad,db,de  ! for nD algorithm
real(kind=8)::s_tot   ! the total length of the IRC for rescaling
real(kind=8)::ss  ! total arclength (for normalization)
real(kind=8)::fprime ! function of extrapolation for first/last points derivative
real(kind=8)::ptsin(points,infdim)  ! input array for nD data to interpolate
real(kind=8)::pts(infdim,points)   ! inverse nD data points
real(kind=8)::s(points)  ! the parameter to model the curve
real(kind=8),allocatable::yv(:),xv(:)  ! arrays with s-values and coordinates 
                                       ! for the respective dimension
real(kind=8),allocatable::u(:),y2(:)  ! y2: second derivative (one dimension)
real(kind=8)::y2_nd(points,infdim)  ! all second derivatives (nD)

!
!     Define the total number of points to be used (number of structures on the 
!     reaction path)
!     Define the dimensionality of the interpolated data (internal coordinates and 
!     energies = nat6+1)
!     and fill the final data arrays
!
!
dim=infdim
n=points
allocate(yv(n),xv(n))  !for 1D calculations
allocate(u(n-1),y2(n))   ! for 1D calculations
!
!     Determine the topology of the path points and translate them into arclengths
!     of the s-values
!
open(unit=18,file="test.plot",status="unknown")
do i=1,points
   write(18,*) ptsin(i,:)
end do
close(18)

!
!     The parameter values s(i) are not calculated here, but are 
!     taken as input parameters from the first spline routine
!     check only if two points have the same or nearly the 
!     same s value
!
do i=1,n-1
   if (s(i+1)-s(i) .lt. 1E-8) then
      write(*,*) "The parameter curve values of two reference points are"
      write(*,*) "(nearly) identical! This should not happen!"
      write(*,*) "Check the file ref.input and restart the calculation."
      call fatal
   end if
end do
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
end subroutine build_spline2

!
!     ---> fprime subroutine is contained in build_spline.f90!
!

