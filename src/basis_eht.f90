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
!     basis_eht: second step of basis setup for EHT: calculate 
!     needed STO-3G functions and fill them into the Hamiltonian!
!     called subroutines:
!     - setsto3,setsto4
!
!     part of QMDFF
!
subroutine basis_eht(n,at,nbf,ok)    
use qmdff       
implicit none
integer::elem,n,nbf
integer::at(n)
logical::ok

integer::i,j,ibf,ipr,p,thisprim
real(kind=8)::a(6),c(6),zeta,k1,k2
real(kind=8)::parz(86,3)

ibf=0
ipr=0
ok=.true.
!
!     equalize arrays to them defined in setZETAandIP 
!     (setcommon.f90)!
!
parip=ip
parz=zet
do i=1,n
!     H-He
   if (at(i).le.2) then
!     1s
      ibf =ibf+1
      zeta=parz(1,1)
      call setsto3(thisprim,1,1,zeta,a,c)
      do p=1,thisprim  
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim 
      hdiag(ibf)=parip(at(i),1)
   end if

!     Li-F
   if (at(i).le.10.and.at(i).gt.2) then
!     2s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,2,1,zeta,a,c)
      do p=1,thisprim
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     2p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,2,2,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),2)
      end do
   end if

!     Na-Al
   if(at(i).le.13.and.at(i).gt.10)then
!     3s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,3,1,zeta,a,c)
      do p=1,thisprim 
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     3p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,3,2,zeta,a,c)
         do p=1,thisprim 
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim 
         hdiag(ibf)=parip(at(i),2)
      end do
   end if

!     Si-Ar
   if(at(i).le.18.and.at(i).gt.13)then
!     3s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,3,1,zeta,a,c)
      do p=1,thisprim 
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     3p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,3,2,zeta,a,c)
         do p=1,thisprim 
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim 
         hdiag(ibf)=parip(at(i),2)
      end do
   end if

!     K,Ca,Zn-Kr 
   if ((at(i).le.36.and.at(i).ge.30).or.&
    &  at(i).eq.19.or. at(i).eq.20) then
!     4s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,4,1,zeta,a,c)
      do p=1,thisprim
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     4p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,4,2,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),2)
      end do
   end if

!     Sc-Cu
   if (at(i).le.29.and.at(i).ge.21) then
!     4s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,4,1,zeta,a,c)
      do p=1,thisprim
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     4p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,4,2,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),2)
      end do
!     3d
      zeta=parz(at(i),3)
      do j=5,10
         ibf=ibf+1
         call setsto4(thisprim,3,3,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
!    dxy...are now correctly normalized
            if(j.gt.7)cont(ipr)=cont(ipr)*sqrt(3.0d0)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),3)
      end do
   end if

!     Y-ag
   if (at(i).le.47.and.at(i).ge.39) then
!     5s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,5,1,zeta,a,c)
      do p=1,thisprim
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      end do
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     5p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,5,2,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),2)
      end do
!     4d
      zeta=parz(at(i),3)
      do j=5,10
         ibf=ibf+1
         call setsto4(thisprim,4,3,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
!     dxy...are now correctly normalized
            if (j.gt.7) cont(ipr)=cont(ipr)*sqrt(3.0d0)
         end do
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),3)
      end do
   endif
!     In-Xe
   if (at(i).gt.47.and.at(i).le.54) then
!     5s
      ibf=ibf+1
      zeta=parz(at(i),1)
      call setsto4(thisprim,5,1,zeta,a,c)
      do p=1,thisprim
         ipr=ipr+1
         alp (ipr)=a(p)
         cont(ipr)=c(p)
      enddo
      fila(1,i) =ibf
      aoat (ibf)=i  
      lao  (ibf)=1  
      nprim(ibf)=thisprim
      hdiag(ibf)=parip(at(i),1)
!     5p
      zeta=parz(at(i),2)
      do j=2,4
         ibf=ibf+1
         call setsto4(thisprim,5,2,zeta,a,c)
         do p=1,thisprim
            ipr=ipr+1
            alp (ipr)=a(p)
            cont(ipr)=c(p)
         enddo
         aoat (ibf)=i  
         lao  (ibf)=j  
         nprim(ibf)=thisprim
         hdiag(ibf)=parip(at(i),2)
      end do
   end if

   fila(2,i)=ibf

enddo

do i=1,ibf
   if (hdiag(i).eq.0) ok=.false.
end do
do i=1,ipr
   if (alp(i).eq.0)  ok=.false.
end do

end subroutine basis_eht
