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
!     subroutine axis: calculate moment of inertia etc.
!
!     part of QMDFF
!    

subroutine axis(numat,nat,xyz,aa,bb,cc,avmom,sumw)        
implicit double precision (a-h,o-z)                        
dimension xyz(3,*), ams(86)                                 
integer::nat(*)
PARAMETER (BOHR=0.52917726)                                               
dimension t(6), rot(3), xyzmom(3), eig(3), evec(3,3)                      
dimension x(numat),y(numat),z(numat),coord(3,numat)
!                                                        
!     const1 =  10**40/(n*a*a)                                                  
!               n = avergadro's number                                          
!               a = cm in an angstrom                                           
!               10**40 is to allow units to be 10**(-40)gram-cm**2              
!                                                                               
const1 = 1.66053d0                                                        
!                                                                               
!     const2 = conversion factor from angstrom-amu to cm**(-1)                  
!                                                                               
!            = (planck's constant*n*10**16)/(8*pi*pi*c)                         
!            = 6.62618*10**(-27)[erg-sec]*6.02205*10**23*10**16/                
!              (8*(3.1415926535)**2*2.997925*10**10[cm/sec])                    
!                                                                               
const2=16.8576522d0 
     
t=0d0                                                 
!
!     table with atomic masses:
!
ams =(/ 1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0, &
   &  10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0, &
   &  20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0, &
   &  30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0, &
   &  40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0, &
   &  54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0, &
   &  65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0, &
   &  79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0, &
   &  91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0, &
   &  102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0, &
   &  118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0, &
   &  132.9054d0, 137.3300d0, &
   &  138.91d0, 140.12d0,140.91d0,144.24d0,147.00d0,150.36d0, &
   &  151.97d0, 157.25d0, 158.93d0,162.50d0,164.93d0,167.26d0, &
   &  168.93d0, 173.04d0,174.97d0, 178.4900d0, 180.9479d0, &
   &  183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0, &
   &  196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0, &
   &  18*0.000d0,   0.0000d0,  5*0.000d0/)
                                                      
sumw=1.d-20                                                               
sumwx=0.d0                                                                
sumwy=0.d0                                                                
sumwz=0.d0                                                                

coord(1:3,1:numat)=xyz(1:3,1:numat)*bohr
                                                                          
do i=1,numat                                                        
   atmass=ams(nat(i))                                                  
   sumw=sumw+atmass                                                    
   sumwx=sumwx+atmass*coord(1,i)                                       
   sumwy=sumwy+atmass*coord(2,i)
   sumwz=sumwz+atmass*coord(3,i)
end do                                                              
sumwx=sumwx/sumw                                                          
sumwy=sumwy/sumw                                                          
sumwz=sumwz/sumw                                                          
f=1.0d0/bohr        
do i=1,numat                                                           
   x(i)=coord(1,i)-sumwx                                                  
   y(i)=coord(2,i)-sumwy                                                  
   z(i)=coord(3,i)-sumwz  
end do                                                   
!
!    matrix for moments of inertia is of form                                   
!                                                                               
!           |   y**2+z**2                         |                             
!           |    -y*x       z**2+x**2             | -i =0                       
!           |    -z*x        -z*y       x**2+y**2 |                             
!                                                                               
do  i=1,6                                                               
   t(i)=dble(i)*1.0d-10  
end do                                                    
do i=1,numat                                                        
   atmass=ams(nat(i))                                                  
   t(1)=t(1)+atmass*(y(i)**2+z(i)**2)                                  
   t(2)=t(2)-atmass*x(i)*y(i)                                          
   t(3)=t(3)+atmass*(z(i)**2+x(i)**2)                                  
   t(4)=t(4)-atmass*z(i)*x(i)                                          
   t(5)=t(5)-atmass*y(i)*z(i)                                          
   t(6)=t(6)+atmass*(x(i)**2+y(i)**2)                                  
end do                                                             
call rsp(t,3,3,eig,evec)                                                  
do i=1,3                                                            
   if(eig(i).lt.3.d-4) then                                            
      eig(i)=0.d0                                                      
      rot(i)=0.d0                                                      
   else                                                                
      rot(i)=2.9979245d+4*const2/eig(i)  
   endif                                                               
   xyzmom(i)=eig(i)*const1 
end do                                               

aa=rot(3)/2.9979245d+4
bb=rot(2)/2.9979245d+4
cc=rot(1)/2.9979245d+4
avmom=1.d-47*(xyzmom(1)+xyzmom(2)+xyzmom(3))/3.
                                                                       
return
                                                                    
end subroutine axis                          
