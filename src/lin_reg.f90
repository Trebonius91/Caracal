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
!     subroutine lin_reg: Performs a linear regression for 
!      a pair of data arrays (x and f(x) values).
!      It returns the slope, the y-intercept and the 
!      correlation coefficient  
!
!      part of EVB
!

subroutine lin_reg(x_data,y_data,npts,slope,y_inter,corre)

implicit none                                   
!     Number of points in the given data arrays
integer,intent(in) :: npts
!     The x coordinates of the given data
real(kind=8), intent(in) :: x_data(npts)
!     The y coordinates of the given data
real(kind=8), intent(in) :: y_data(npts)
!     The slope of the y data 
real(kind=8), intent(out) :: slope 
!     The y-intercept for the y data
real(kind=8), intent(out) :: y_inter
!     The correlation coefficient of the given data
real(kind=8), intent(out) :: corre
!     squared correlation coefficient 
real(kind=8) ::  r    
!     intermediate sums for calculations 
real(kind=8) :: sumx,sumy,sumx2,sumy2,sumxy                       
!     actual x and y coordinates 
real(kind=8) :: x,y
!     loop index
integer :: i
!
!     Set default values for variables 
!
slope=0.d0
y_inter=0.d0
corre=0.d0
sumx=0.d0
sumy=0.d0
sumx2=0.d0
sumy2=0.d0
sumxy=0.d0
!
!     Perform loop over all data pairs and calculate values 
!
do i=1,npts
!
!     read out actual data from array
!
   x=x_data(i)
   y=y_data(i)
!
!     compute sum of x
!
   sumx  = sumx + x  
!
!     compute sum of x**2
!
   sumx2 = sumx2 + x * x  
!
!     compute sum of x * y
!
   sumxy = sumxy + x * y
!
!     compute sum of y
!
   sumy  = sumy + y  
!
!     compute sum of y**2
!
   sumy2 = sumy2 + y * y
end do
!
!     compute slope
!
slope = (npts * sumxy  -  sumx * sumy) / (npts * sumx2 - sumx**2)   
!
!     compute y-intercept
!
y_inter = (sumy * sumx2  -  sumx * sumxy) / (npts * sumx2  -  sumx**2) 
!
!     compute correlation coefficient
!
corre = (sumxy-sumx*sumy/npts)/sqrt((sumx2-sumx**2/npts)*(sumy2-sumy**2/npts))

return
end subroutine lin_reg
