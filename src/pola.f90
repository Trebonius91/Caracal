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
!     subroutine pola: calculate product of two EHT gaussian basis functions
!
!     part of QMDFF
!
subroutine pola(a,b,ga,gb,gm,gm2,iff1,iff2,iall, &
   &   aa,bb,lmnexp,lmnfak,est,arg,va)
implicit none
!     aufpunkte,intarray
real*8  a(3),b(3),va,est,arg,lmnfak(84)
real*8  gm,gm2,ga,gb
integer iall,lmnexp(84),iff1,iff2
!     local
real*8 e(3),d(3),olap3,efact,val
real*8 aa(10),bb(10),dd(35)
integer i
!
!     if you want f-functions: dimension aa(20),bb(20),dd(84)
!
!     a,b are centres of gaussians                                              
!     ga,gb are their exponents                                                 
!     gama=ga+gb          
!     gm =1.0d0/gama
!     gm2=0.5*gm
!     apply product theorem                                                     
!     e is center of product gaussian with exponent gama    
!                  
do i=1,3
   e(i)=(ga*a(i)+gb*b(i))*gm
enddo
!
!     calculate cartesian prefactor for first gaussian             
!             
call rhftce(aa,a,e,iff1)
!                                      
!     calculate cartesian prefactor for second gaussian  
!                       
call rhftce(bb,b,e,iff2)
!                                                
!     form their product       
!     dd=0
!                                               
call prod(aa,bb,dd,iff1,iff2)
val = 0
do i=1,iall
   olap3=lmnfak(i)*arg*gm2**lmnexp(i)
   val=dd(i)*olap3+val
enddo
va=exp(-est)*val

return
end subroutine pola

