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
!     subroutine covbond: determine bond topology for EHT model system
!         all atoms outside torsion are artificial H atoms
!
!     part of QMDFF
! 
subroutine covbond(n,at,xyz,pair)
use qmdff
implicit none
integer::n,at(*),pair(*)
real(kind=8)::xyz(3,*)

real(kind=8)::r,rmin,rco,dum,aa,bb,f
integer::nn(n),list(1000,n),jj,nlist(1000,n),i2,k,iat
integer::i,ni,newi,ii,i1,ni1,iii,d,j,newatom,nj,m,pqn
integer::imem,jmem,nnn(n),l,tag,irow,jrow,nb(20,n)
logical::da,dai,daj

nb=0
do i=1,n
   f=1.3
   k=0
   100 continue
   do iat=1,n
      da=.false.
      do j=1,k
         if (nb(j,i).eq.iat) da=.true.
      end do
!
!     add bond if radius is small enough and no too many bonds exist
!
      if (iat.ne.i.and.(.not.da))then
         r=sqrt((xyz(1,iat)-xyz(1,i))**2 &
           &   +(xyz(2,iat)-xyz(2,i))**2 &
           &   +(xyz(3,iat)-xyz(3,i))**2) &
           &   *0.52917726
         rco=rad(at(i))+rad(at(iat))
         if (r.lt.f*rco.and.k.lt.19) then
            k=k+1
            nb(k,i)=iat
         end if
      end if
   end do
   if (k.lt.1) then
      f=f*1.1
      goto 100
   end if
   nb(20,i)=k
end do

nn(1:n)=nb(20,1:n)

pair(1:n*(n+1)/2)=0

list=0
do i=1,n
   ni=nn(i)
   list(1:ni,i)=nb(1:ni,i)
end do

nlist=list
!
!     one bond, tag=1      
!
tag=1
call pairs(n,nn,list,pair,tag)
!
!     determine up to 4 bonds in between      
! 

do d=1,3
!
!     loop over atoms    
!  
   do i=1,n
      ni=nn(i)
      newi=ni
!
!     all neighbors of atom i    
!     
      do ii=1,ni
         i1=list(ii,i)
         ni1=nb(20,i1)
!
!     all neighbors of neighbors of atom i   
!      
         do iii=1,ni1
            newatom=nb(iii,i1)
            da=.false.
            do j=1,newi
               if (newatom.eq.list(j,i)) da=.true.
            end do
            if (.not.da) then
               newi=newi+1
               nlist(newi,i)=newatom           
            end if
         end do
      end do
      nnn(i)=newi
   end do

   list=nlist
   nn=nnn
!
!     one bond more
!
   tag=tag+1
   call pairs(n,nn,list,pair,tag)

end do
!
!     5 tags five or more covalent bonds
!
do i=1,n*(n+1)/2
   if(pair(i).eq.0) pair(i)=5
enddo

return
end subroutine covbond
