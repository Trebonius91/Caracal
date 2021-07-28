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
!     subroutine pfit: During the QMDFF setup, the torsional force 
!     constants are obtained by fitting the FF torsional potential 
!     to EHT reference
!
!     part of QMDFF
!
subroutine pfit(n,y,tangle,damp,s,o,p,mad,fx,np)       
implicit none
integer::n,np        
real(kind=8)::tangle,y(n),s(*),o(*),p(*),damp,fx(n)

integer::maxcycle, iiter, i, L, j,errorflag

real(kind=8)::rr(n),ll(n),aa(n,np),w(n)
real(kind=8)::dp(np),Nmat(np,np),tt(np,n)
real(kind=8)::step,inkre,thr,maxd,s0,s1,dsq,sumx,det,mad
!
!     start damping (=20%)
!
inkre=0.01
!
!     exit threshold (rel. change in RMSD)
!
thr=1.d-6

w = 1
!cccccccccccccccccc
!     ALGO STARTS HERE    
!     --> do no more than 50 iteration  
!cccccccccccccccccc

s1=0              
do iiter=1,50           

   call getxy(np,s,o,p,tangle,damp,n,fx)

   s0=s1
   s1=0
   do i=1,n
      s1=s1+(w(i)*(y(i)-fx(i)))**2
   end do
   dsq=s1-s0
  
   if(dsq.lt.0.and.iiter.gt.1) inkre=inkre/1.2
!     write(*,'('' iter'',i4,'' RMSD ='',F9.4,5x,''rel. change'',F9.4,
!    .          5x,''damping'',F9.4)')
!    .iiter,s1,dsq/s1,inkre

   if(abs(dsq/s1).lt.thr.or.s1.lt.thr) exit
!
!     parameter loop for num gradients
!
   do i=1,np
      step=abs(0.001*p(i))
      step=max(step,0.0001d0)

      p(i)=p(i)+step
      call getxy(np,s,o,p,tangle,damp,n,rr)

      p(i)=p(i)-2.*step
      call getxy(np,s,o,p,tangle,damp,n,ll)

      p(i)=p(i)+step

      do j=1,n
         aa(j,i)=(rr(j)-ll(j))/(2.*step)
      end do
  
   end do

   do i=1,n
      do j=1,np
         tt(j,i)=aa(i,j)
      end do
   end do
!
!     (* T^*A^ = Nmat *)
!
   do i = 1, np 
      do j = 1, np 
         sumx = 0
         do L = 1, n 
            sumx = sumx+tt(i, L)*aa(L, j)
         end do
         Nmat(i, j) = sumx
      end do
   end do
!
!      damping by diagonal shift
!
   do i=1,np
      Nmat(i,i)=Nmat(i,i)+inkre*Nmat(i,i)
   end do
!
!      (*  T^ * residuen = dp *)
!
   do i = 1, np 
      sumx = 0
      do L = 1, n 
         sumx = sumx+tt(i, L)*(y(L)-fx(L))*w(L)
      end do
      dp(i) = sumx
   end do
!
!      (* Matrix invert *)
!
   call dmatinv(Nmat,np,np,det)
!
!      (* N^-1 * (transp.*residuen) *)
!
   do i = 1, np 
      sumx = 0
      do L = 1, np 
         sumx = sumx+Nmat(i, L)*dp(L)
      end do
!
!      (* new Parameter *)
!
      p(i)=p(i)+sumx
   end do

end do

99   continue     

mad=0
do i=1,n
   sumx=sumx+abs(y(i))
   mad=mad+abs(y(i)-fx(i))
end do
mad=mad/sumx

return
end subroutine pfit

