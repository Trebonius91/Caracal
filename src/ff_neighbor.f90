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
!     subroutine ff_neighbor: determine atom neighbors for the 
!        QMDFF reference structure
!
!     part of QMDFF
!
!
subroutine ff_neighbor(natoms,xyz,wbo,at,nb,hyb,mclust,wbocut_in)
use qmdff
implicit none  
integer::at(natoms),natoms,nb(20,natoms),hyb(natoms)
real(kind=8)::xyz(3,natoms),wbo(natoms,natoms),wbocut_in
logical::mclust

logical::da
character(len=40)::atmp
integer::iat,i,j,k,nn,ni,ii,jj,kk,ll
real(kind=8)::rb(20,natoms),valijkl,wbocut
real(kind=8)::dx,dy,dz,r,damp,xn,rr,rco,r2,f,aa1

write(10,*) 
write(10,*)'======================'
write(10,*)'topology analysis:'
write(10,*)'find nearest neighbors'
write(10,*)'======================'
write(10,*) 

nb=0
nn=min(natoms,2)-1
!
!write(10,*) "Topology analysis for molecule:"
!write(10,*) "1st  2nd  r_covalent  distance  bond-order bond-index"
do i=1,natoms
   f=1.3
   k=0
   100 continue
   do iat=1,natoms
      da=.false.
      do j=1,k
         if (nb(j,i).eq.iat) da=.true.
      end do
      if (iat.ne.i.and.(.not.da)) then
         dx=xyz(1,iat)-xyz(1,i)
         dy=xyz(2,iat)-xyz(2,i)
         dz=xyz(3,iat)-xyz(3,i)
         r2=dx*dx+dy*dy+dz*dz 
         r=sqrt(r2)*0.52917726
         rco=rad(at(i))+rad(at(iat))
!
!     important cut-off for small BE bonds     
!       
         wbocut=wbocut_in
         if (metal(at(i)).eq.1.or.metal(at(iat)).eq.1) then
            wbocut=-1
!
!     this makes high CN (coordination number?)    
!       
            if (mclust) f=1.4   
         end if
!
!     critical step: set Wiberg bond orders and coordinatiom numbers
!     together (only if they are lower than 19!)     
!      
         if (r.lt.f*rco.and.wbo(iat,i).gt.wbocut.and.k.lt.19) then
            k=k+1
            nb(k,i)=iat
            rb(k,i)=r  
         end if
!         write(10,'(2i4,3f10.5,i5)') i,iat,rco,r,wbo(i,iat),k
      end if
   end do
!
!     Rescale f factor to gain success and go back
!
   if (k.lt.1.and.f.lt.1.5) then
      f=f*1.1
      goto 100
   end if
   nb(20,i)=k
!
!     sort atoms and bonds to them (?)
!
   call sort(1,k,rb(1,i),nb(1,i))
end do
!
!     For lonely atoms (no coordination and WBO)
!
do i=1,natoms
   if(nb(20,i).eq.0)then
      write(atmp,'(''no bond partners for atom'',i4)')i
      call warn(atmp)
   endif
enddo
!
!     tag atoms in linear coordination 
!     
do i=1,natoms
   ni=nb(20,i)
   if(ni.eq.2)then
      k=nb(1,i)
      j=nb(2,i)
      call bangle(XYZ,k,i,j,aa1)          
      aa1=aa1*180./pi
      if(aa1.lt.athr.or.180-aa1.lt.athr) nb(19,i)=1
   endif
enddo
!
!     tag atoms bonded to metals    
!        
do i=1,natoms
   ni=nb(20,i)
   do j=1,ni
      k=nb(j,i)
      if (metal(at(k)).eq.1) nb(18,i)=1
   end do
end do
!
!     hybridization state for group 3-7 elements 
!     first main group: sp2 hybridization etc.
!     and for their heavier homologes depending on the 
!     number of bond partners
!
do i=1,natoms
   hyb(i)=0
!     B (Bor)         
   if (at(i).eq.5.or.at(i).eq.13.or.at(i) &
      &  .eq.31.or.at(i).eq.49.or.at(i).eq.81) then
      if (nb(20,i).eq.4) hyb(i)=3
      if (nb(20,i).eq.3) then       
         hyb(i)=3
         ll=nb(1,i)
         jj=nb(2,i)
         kk=nb(3,i)
         r=valijkl(natoms,xyz,ll,i,jj,kk)     
         if (abs(r*180./pi-180.).lt.athr) hyb(i)=2
      end if
   end if
!
!     C (Carbon)      
!    
   if(at(i).eq.6.or.at(i).eq.14.or.at(i) &
      &   .eq.32.or.at(i).eq.50.or.at(i).eq.82) then
      if (nb(20,i).eq.4) hyb(i)=3
      if (nb(20,i).eq.3) hyb(i)=2
      if (nb(20,i).eq.2) hyb(i)=1
      if (nb(20,i).eq.4.and.at(i).eq.6.and.nb(18,i).eq.1) hyb(i)=2
   end if
!
!     N (Nitrogen)  
!       
   if (at(i).eq.7.or.at(i).eq.15.or.at(i) &
      &   .eq.33.or.at(i).eq.51.or.at(i).eq.83) then
      if (nb(20,i).eq.4) hyb(i)=3   
      if (nb(20,i).eq.3) then       
         hyb(i)=3
         ll=nb(1,i)
         jj=nb(2,i)
         kk=nb(3,i)
         r=valijkl(natoms,xyz,ll,i,jj,kk)     
         if (abs(r*180./pi-180.).lt.athr) hyb(i)=2
      end if
      if (nb(20,i).eq.2) hyb(i)=2
      if (nb(20,i).eq.1) hyb(i)=1
   endif 
!
!     O (Oxygen)    
!     
   if(at(i).eq.8.or.at(i).eq.16.or.at(i) &
      &   .eq.34.or.at(i).eq.52.or.at(i).eq.84) then
      if (nb(20,i).eq.3) hyb(i)=3
      if (nb(20,i).eq.2) hyb(i)=3
      if (nb(20,i).eq.1) hyb(i)=2
   end if
end do

return
end subroutine ff_neighbor
