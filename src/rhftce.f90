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
!     subroutine rhftce: calculate cartesian prefactors for EHT gaussians
!
!     part of QMDFF
!
subroutine rhftce(cfs,a,e,iff)                                            
implicit real*8(a-h,o-z)                                                  
dimension cfs(*) ,a(*),e(*)                                               
data c2/2.0d0/,c3/3.0d0/   
!                                                 
!     e = center of product function, a = center of single gaussian     
!       
aex = e(1)-a(1)                                                           
aey = e(2)-a(2)                                                           
aez = e(3)-a(3)  
select case (iff)
   case(1)
   case(2)
      cfs(1)=aex*cfs(2) 
   case(3)
      cfs(1)=aey*cfs(3)
   case(4)
      cfs(1)=aez*cfs(4)
   case(5)
      cfs(1)=aex*aex*cfs(5)
      cfs(2)=c2*aex*cfs(5)
   case(6)
      cfs(1)=aey*aey*cfs(6)
      cfs(3)=c2*aey*cfs(6)
   case(7)
      cfs(1)=aez*aez*cfs(7)
      cfs(4)=c2*aez*cfs(7)
   case(8)
      cfs(1)=aex*aey*cfs(8)
      cfs(2)=aey*cfs(8)
      cfs(3)=aex*cfs(8)
   case(9)
      cfs(1)=aex*aez*cfs(9)
      cfs(2)=aez*cfs(9)
      cfs(4)=aex*cfs(9) 
   case(10)
      cfs(1)=aey*aez*cfs(10)
      cfs(3)=aez*cfs(10)
      cfs(4)=aey*cfs(10)
   case(11)
      cfs(1)=aex*aex*aex*cfs(11)
      cfs(2)=c3*aex*aex*cfs(11)
      cfs(5)=c3*aex*cfs(11)
   case(12)
      cfs(1)=aey*aey*aey*cfs(12)
      cfs(3)=c3*aey*aey*cfs(12)                 
      cfs(6)=c3*aey*cfs(12)      
   case(13)
      cfs(1)=aez*aez*aez*cfs(13)
      cfs(4)=c3*aez*aez*cfs(13)                
      cfs(7)=c3*aez*cfs(13) 
   case(14)
      cfs(1)=aex*aex*aey*cfs(14)            
      cfs(2)=c2*aex*aey*cfs(14)
      cfs(3)=aex*aex*cfs(14)                                           
      cfs(5)=aey*cfs(14)                                             
      cfs(8)=c2*aex*cfs(14)
   case(15)
      cfs(1)=aex*aex*aez*cfs(15)                               
      cfs(2)=c2*aex*aez*cfs(15)                   
      cfs(4)=aex*aex*cfs(15)                                        
      cfs(5)=aez*cfs(15)  
      cfs(9)=c2*aex*cfs(15)
   case(16)
      cfs(1)=aey*aey*aex*cfs(16)
      cfs(2)=aey*aey*cfs(16)
      cfs(3)=c2*aey*aex*cfs(16)
      cfs(6)=aex*cfs(16)
      cfs(8)=c2*aey*cfs(16)
   case(17)
      cfs(1)=aey*aey*aez*cfs(17)
      cfs(3)=c2*aey*aez*cfs(17)
      cfs(4)=aey*aey*cfs(17)
      cfs(6)=aez*cfs(17)
      cfs(10)=c2*aey*cfs(17)
   case(18)
      cfs(1)=aez*aez*aex*cfs(18)
      cfs(2)=aez*aez*cfs(18)
      cfs(4)=c2*aez*aex*cfs(18)
      cfs(7)=aex*cfs(18)
      cfs(9)=c2*aez*cfs(18)
   case(19)
      cfs(1)=aez*aez*aey*cfs(19)
      cfs(3)=aez*aez*cfs(19)
      cfs(4)=c2*aez*aey*cfs(19)
      cfs(7)=aey*cfs(19)
      cfs(10)=c2*aez*cfs(19)
   case(20)
      cfs(1)=aex*aey*aez*cfs(20)
      cfs(2)=aez*aey*cfs(20)
      cfs(3)=aex*aez*cfs(20)
      cfs(4)=aex*aey*cfs(20)
      cfs(8)=aez*cfs(20)
      cfs(9)=aey*cfs(20)
      cfs(10)=aex*cfs(20)
end select
return                                                       
end subroutine rhftce
