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
!     prod: calculate product of two EHT gaussian basis functions
!       distinguish between difference combinations of angle functions
!
!     part of QMDFF
!
subroutine prod(c,d,s,iff1,iff2)                                          
implicit real(kind=8) (a-h,o-z)                                                  
dimension c(*),d(*),s(*)  
!     too large values                                                
if(iff1.gt.10.or.iff2.gt.10) goto 30      
!     s-s product                                
if (iff1.gt.4.or.iff2.gt.4) goto 20                                        
s( 1)=c( 1)*d( 1)                                                         
!     end of s - s      
!     s-p product                                               
if(iff1.eq.1.and.iff2.eq.1) return                                        
s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
!     end of s - p                                                     
if(iff1.eq.1.or.iff2.eq.1) return                                         
s( 5)=c( 2)*d( 2)                                                         
s( 6)=c( 3)*d( 3)                                                         
s( 7)=c( 4)*d( 4)                                                         
s( 8)=c( 2)*d( 3)+c( 3)*d( 2)                                             
s( 9)=c( 2)*d( 4)+c( 4)*d( 2)                                             
s(10)=c( 3)*d( 4)+c( 4)*d( 3)                                             
!     end of p - p                                                     
return                                                                    
20 continue                                                                  
s( 1)=c( 1)*d( 1)                                                         
s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)                                 
s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)                                 
s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)                                 
s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)                     
s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)                     
s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)                     
!     end of s - d                                                     
if(iff1.eq.1.or.iff2.eq.1) return                                         
s(11)=c( 2)*d( 5)+c( 5)*d( 2)                                             
s(12)=c( 3)*d( 6)+c( 6)*d( 3)                                             
s(13)=c( 4)*d( 7)+c( 7)*d( 4)                                             
s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)                     
s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)                     
s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)                     
s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)                     
s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)                     
s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)                     
s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9)  &                               
  &  +c( 9)*d( 3)+c( 4)*d( 8)+c( 8)*d( 4)                                 
!     end of p - d                                                     
if(iff1.lt.5.or.iff2.lt.5) return                                         
s(21)=c( 5)*d( 5)                                                         
s(22)=c( 6)*d( 6)                                                         
s(23)=c( 7)*d( 7)                                                         
s(24)=c( 5)*d( 8)+c( 8)*d( 5)                                             
s(25)=c( 5)*d( 9)+c( 9)*d( 5)                                             
s(26)=c( 6)*d( 8)+c( 8)*d( 6)                                             
s(27)=c( 6)*d(10)+c(10)*d( 6)                                             
s(28)=c( 7)*d( 9)+c( 9)*d( 7)                                             
s(29)=c( 7)*d(10)+c(10)*d( 7)                                             
s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)                                 
s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)                                 
s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)                                 
s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)                     
s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)                     
s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)                     
!     end of d - d                                                     
return                                                                    
30 continue                                                                  
s( 1)=c( 1)*d( 1)                                                         
s( 2)=c( 1)*d( 2)+c( 2)*d( 1)                                             
s( 3)=c( 1)*d( 3)+c( 3)*d( 1)                                             
s( 4)=c( 1)*d( 4)+c( 4)*d( 1)                                             
s( 5)=c( 1)*d( 5)+c( 5)*d( 1)+c( 2)*d( 2)                                 
s( 6)=c( 1)*d( 6)+c( 6)*d( 1)+c( 3)*d( 3)                                 
s( 7)=c( 1)*d( 7)+c( 7)*d( 1)+c( 4)*d( 4)                                 
s( 8)=c( 1)*d( 8)+c( 8)*d( 1)+c( 2)*d( 3)+c( 3)*d( 2)                     
s( 9)=c( 1)*d( 9)+c( 9)*d( 1)+c( 2)*d( 4)+c( 4)*d( 2)                     
s(10)=c( 1)*d(10)+c(10)*d( 1)+c( 3)*d( 4)+c( 4)*d( 3)                     
s(11)=c( 2)*d( 5)+c( 5)*d( 2)+c(1)*d(11)+c(11)*d(1)                       
s(12)=c( 3)*d( 6)+c( 6)*d( 3)+c(1)*d(12)+c(12)*d(1)                       
s(13)=c( 4)*d( 7)+c( 7)*d( 4)+c(1)*d(13)+c(13)*d(1)                       
s(14)=c( 2)*d( 8)+c( 8)*d( 2)+c( 3)*d( 5)+c( 5)*d( 3)+c(1)*d(14)+ &        
 &    c(14)*d(1)                                                          
s(15)=c( 2)*d( 9)+c( 9)*d( 2)+c( 4)*d( 5)+c( 5)*d( 4)+c(1)*d(15)+ &        
 &    c(15)*d(1)                                                          
s(16)=c( 2)*d( 6)+c( 6)*d( 2)+c( 3)*d( 8)+c( 8)*d( 3)+c(1)*d(16)+ &        
 &     c(16)*d(1)                                                          
s(17)=c( 3)*d(10)+c(10)*d( 3)+c( 4)*d( 6)+c( 6)*d( 4)+c(1)*d(17)+ &        
 &    c(17)*d(1)                                                          
s(18)=c( 2)*d( 7)+c( 7)*d( 2)+c( 4)*d( 9)+c( 9)*d( 4)+c(1)*d(18)+ &        
 &    c(18)*d(1)                                                          
s(19)=c( 3)*d( 7)+c( 7)*d( 3)+c( 4)*d(10)+c(10)*d( 4)+c(1)*d(19)+ &        
 &    c(19)*d(1)                                                          
s(20)=c( 2)*d(10)+c(10)*d( 2)+c( 3)*d( 9)+c(9)*d(3)+c(4)*d(8)+ &           
 &    c(8)*d(4)+c(1)*d(20)+c(20)*d(1)                                     
!      end of s - f                                                     
if(iff1.eq.1.or.iff2.eq.1) return                                         
s(21)=c( 5)*d( 5)+c(2)*d(11)+c(11)*d(2)                                   
s(22)=c( 6)*d( 6)+c(3)*d(12)+c(12)*d(3)                                   
s(23)=c( 7)*d( 7)+c(4)*d(13)+c(13)*d(4)                                   
s(24)=c( 5)*d( 8)+c( 8)*d( 5)+c(3)*d(11)+c(11)*d(3)+c(2)*d(14)+ &        
  &   c(14)*d(2)                                                          
s(25)=c( 5)*d( 9)+c( 9)*d( 5)+c(2)*d(15)+c(15)*d(2)+c(4)*d(11)+ &          
  &   c(11)*d(4)                                                          
s(26)=c( 6)*d( 8)+c( 8)*d( 6)+c(2)*d(12)+c(12)*d(2)+c(3)*d(16)+ &          
  &   c(16)*d(3)                                                          
s(27)=c( 6)*d(10)+c(10)*d( 6)+c(3)*d(17)+c(17)*d(3)+c(4)*d(12)+ &          
  &   c(12)*d(4)                                                          
s(28)=c( 7)*d( 9)+c( 9)*d( 7)+c(2)*d(13)+c(13)*d(2)+c(4)*d(18)+ &          
  &   c(18)*d(4)                                                          
s(29)=c( 7)*d(10)+c(10)*d( 7)+c(3)*d(13)+c(13)*d(3)+c(4)*d(19)+ &          
  &   c(19)*d(4)                                                          
s(30)=c( 5)*d( 6)+c( 6)*d( 5)+c( 8)*d( 8)+c(2)*d(16)+c(16)*d(2)+ &         
  &   c(3)*d(14)+c(14)*d(3)                                               
s(31)=c( 5)*d( 7)+c( 7)*d( 5)+c( 9)*d( 9)+c(2)*d(18)+c(18)*d(2)+ &         
  &   c(4)*d(15)+c(15)*d(4)                                               
s(32)=c( 6)*d( 7)+c( 7)*d( 6)+c(10)*d(10)+c(3)*d(19)+c(19)*d(3)+ &         
  &   c(4)*d(17)+c(17)*d(4)                                               
s(33)=c( 5)*d(10)+c(10)*d( 5)+c( 8)*d( 9)+c( 9)*d( 8)+c(3)*d(15)+ &        
  &   c(15)*d(3)+c(4)*d(14)+c(14)*d(4)+c(2)*d(20)+c(20)*d(2)              
s(34)=c( 6)*d( 9)+c( 9)*d( 6)+c( 8)*d(10)+c(10)*d( 8)+c(2)*d(17)+ &        
  &   d(2)*c(17)+c(3)*d(20)+c(20)*d(3)+c(4)*d(16)+c(16)*d(4)              
s(35)=c( 7)*d( 8)+c( 8)*d( 7)+c( 9)*d(10)+c(10)*d( 9)+c(2)*d(19)+ &        
  &   c(19)*d(2)+c(3)*d(18)+c(18)*d(3)+c(4)*d(20)+c(20)*d(4)              
!     end of p - f                                                     
if(iff1.eq.2.or.iff2.eq.2) return                                         
s(36)=c(5)*d(11)+c(11)*d(5)                                               
s(37)=c(6)*d(12)+c(12)*d(6)                                               
s(38)=c(7)*d(13)+c(13)*d(7)                                               
s(39)=c(6)*d(11)+c(11)*d(6)+c(5)*d(16)+c(16)*d(5)+c(8)*d(14)+ &            
  &   c(14)*d(8)                                                          
s(40)=c(7)*d(11)+c(11)*d(7)+c(5)*d(18)+c(18)*d(5)+c(9)*d(15)+ &            
  &   c(15)*d(9)                                                          
s(41)=c(5)*d(12)+c(12)*d(5)+c(6)*d(14)+c(14)*d(6)+c(8)*d(16)+ &            
  &   c(16)*d(8)                                                          
s(42)=c(5)*d(13)+c(13)*d(5)+c(7)*d(15)+c(15)*d(7)+c(9)*d(18)+ &            
  &   c(18)*d(9)                                                          
s(43)=c(7)*d(12)+c(12)*d(7)+c(6)*d(19)+c(19)*d(6)+c(10)*d(17)+ &           
  &   c(17)*d(10)                                                         
s(44)=c(6)*d(13)+c(13)*d(6)+c(7)*d(17)+c(17)*d(7)+c(10)*d(19)+ &           
  &   c(19)*d(10)                                                         
s(45)=c(8)*d(11)+c(11)*d(8)+c(5)*d(14)+c(14)*d(5)                         
s(46)=c(9)*d(11)+c(11)*d(9)+c(5)*d(15)+c(15)*d(5)                         
s(47)=c(8)*d(12)+c(12)*d(8)+c(6)*d(16)+c(16)*d(6)                         
s(48)=c(10)*d(12)+c(12)*d(10)+c(6)*d(17)+c(17)*d(6)                       
s(49)=c(10)*d(13)+c(13)*d(10)+c(7)*d(19)+c(19)*d(7)                       
s(50)=c(9)*d(13)+c(13)*d(9)+c(7)*d(18)+c(18)*d(7)                         
s(51)=c(8)*d(13)+c(13)*d(8)+c(7)*d(20)+c(20)*d(7)+c(9)*d(19)+ &            
  &   c(19)*d(9)+c(10)*d(18)+c(18)*d(10)                                  
s(52)=c(10)*d(11)+c(11)*d(10)+c(5)*d(20)+c(20)*d(5)+c(9)*d(14)+ &          
  &   c(14)*d(9)+c(8)*d(15)+c(15)*d(8)                                    
s(53)=c(9)*d(12)+c(12)*d(9)+c(6)*d(20)+c(20)*d(6)+c(10)*d(16)+ &           
  &   c(16)*d(10)+c(8)*d(17)+c(17)*d(8)                                   
s(54)=c(5)*d(17)+c(17)*d(5)+c(6)*d(15)+c(15)*d(6)+c(14)*d(10)+ &           
  &   d(14)*c(10)+c(9)*d(16)+c(16)*d(9)+c(8)*d(20)+c(20)*d(8)             
s(55)=c(5)*d(19)+c(19)*d(5)+c(7)*d(14)+c(14)*d(7)+c(10)*d(15)+ &           
  &   c(15)*d(10)+c(9)*d(20)+c(20)*d(9)+c(8)*d(18)+c(18)*d(8)             
s(56)=c(6)*d(18)+c(18)*d(6)+c(7)*d(16)+c(16)*d(7)+c(10)*d(20)+ &           
  &   c(20)*d(10)+c(9)*d(17)+c(17)*d(9)+c(8)*d(19)+c(19)*d(8)             
!     end of d - f                                                     
if(iff1.eq.3.or.iff2.eq.3) return                                         
s(57)=c(11)*d(11)                                                         
s(58)=c(12)*d(12)                                                         
s(59)=c(13)*d(13)                                                         
s(60)=c(11)*d(12)+c(12)*d(11)+c(14)*d(16)+c(16)*d(14)                     
s(61)=c(11)*d(13)+c(13)*d(11)+c(15)*d(18)+c(18)*d(15)                     
s(62)=c(12)*d(13)+c(13)*d(12)+c(17)*d(19)+c(19)*d(17)                     
s(63)=c(11)*d(14)+c(14)*d(11)                                             
s(64)=c(11)*d(15)+c(15)*d(11)                                             
s(65)=c(13)*d(18)+c(18)*d(13)                                             
s(66)=c(13)*d(19)+c(19)*d(13)                                             
s(67)=c(12)*d(17)+d(12)*c(17)                                             
s(68)=c(12)*d(16)+c(16)*d(12)                                             
s(69)=c(11)*d(16)+c(16)*d(11)+c(14)*d(14)                                 
s(70)=c(11)*d(18)+c(18)*d(11)+c(15)*d(15)                                 
s(71)=c(13)*d(15)+c(15)*d(13)+c(18)*d(18)                                 
s(72)=c(13)*d(17)+c(17)*d(13)+c(19)*d(19)                                 
s(73)=c(12)*d(14)+c(14)*d(12)+c(16)*d(16)                                 
s(74)=c(12)*d(19)+c(19)*d(12)+c(17)*d(17)                                 
s(75)=c(11)*d(17)+c(17)*d(11)+c(14)*d(20)+c(20)*d(14)+ &                   
  &   c(15)*d(16)+c(16)*d(15)                                             
s(76)=c(11)*d(19)+c(19)*d(11)+c(20)*d(15)+c(15)*d(20)+ &                   
  &   c(14)*d(18)+c(18)*d(14)                                             
s(77)=c(12)*d(18)+c(18)*d(12)+c(16)*d(19)+c(19)*d(16)+ &                   
  &   c(17)*d(20)+c(20)*d(17)                                             
s(78)=c(13)*d(14)+c(14)*d(13)+c(15)*d(19)+c(19)*d(15)+ &                   
  &   c(18)*d(20)+c(20)*d(18)                                             
s(79)=c(12)*d(15)+c(15)*d(12)+c(14)*d(17)+c(17)*d(14)+ &                   
  &   c(16)*d(20)+c(20)*d(16)                                             
s(80)=c(13)*d(16)+c(16)*d(13)+c(17)*d(18)+c(18)*d(17)+ &                   
  &   c(19)*d(20)+c(20)*d(19)                                             
s(81)=c(11)*d(20)+c(20)*d(11)+c(14)*d(15)+c(15)*d(14)                     
s(82)=c(12)*d(20)+c(20)*d(12)+c(16)*d(17)+c(17)*d(16)                     
s(83)=c(13)*d(20)+c(20)*d(13)+c(18)*d(19)+c(19)*d(18)                     
s(84)=c(14)*d(19)+c(19)*d(14)+c(15)*d(17)+c(17)*d(15)+  &               
  &   c(16)*d(18)+c(18)*d(16)+c(20)*d(20)                                 
!      end of f - f                                                     
return                                                                    
end subroutine prod 
