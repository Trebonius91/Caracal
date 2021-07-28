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
!     subroutine ff_hb: energies and gradients of hydrogen/halogen
!     bond terms 
!
!     part of QMDFF
!

subroutine ff_hb(n,at,xyz,eh,g)
use qmdff
implicit none
integer::n,at(*)
real(kind=8)::xyz(3,n),g(3,n),eh
integer::i1,i2,i,k,j
real(kind=8)::c1,c2,r,er,el,step,edum,dum
integer::at_A,at_H  ! halogen bond partners
real(kind=8)::dum1,dum2
real(kind=8)::ri,rj,cpar
logical::halogen,hydrogen
real(kind=8)::rbond(3)

!write(*,*) "e before",eh
if (nhb.lt.1 .and. nmols.eq.0) return
!
!     loop over all hydrogen/halogen bonds:
!     first determine distance between donor/acceptor
!     count only if distance is smaller than 15 bohr
!
do k=1,nhb
   i1 =hb(1,k)
   i2 =hb(2,k)
   rbond(1)=xyz(1,i1)-xyz(1,i2)
   rbond(2)=xyz(2,i1)-xyz(2,i2)
   rbond(3)=xyz(3,i1)-xyz(3,i2)

   if (periodic) then
      call box_image(rbond)
   end if

   r=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)

   if (r.gt.15.0d0) cycle
   i=hb(3,k)
!
!     if an H-atom takes part, calculate hydrogen bond
!
   if (at(i).eq.1) then
      c1 =vhb(1,k)
      c2 =vhb(2,k)
      call eabhag(n,i1,i2,i,xyz,c1,c2,eh,g)
   else
      c1 =vhb(1,k)
      call eabxag(n,i1,i2,i,xyz,c1,eh,g)
   end if
end do
!
!    Add a second component for the inter-molecule interaction for solvent QMDFFs
!    Loop over all bonds in the system and check if 
!    a) One atom is a F, Cl, Br, I --> halogen bond
!    b) One atom is a H and the other is a N,O,F,S,Cl --> hydrogen bond
!
if (nmols .gt. 1) then
   do i=1,nbond
      i1 = bond(1,i)
      i2 = bond(2,i) 
!
!    Detect a Halogen bond and assign the roles of the halogen and acceptor atom
!    Only if the second atom is no hydrogen!
!

      if ((at(i1).eq.17) .or. (at(i1).eq.35) .or. (at(i1).eq.53) .or. & 
               &  (at(i1) .eq. 85)) then
         if (at(i2) .ne. 1) then
            at_H=i1
            at_A=i2
            halogen=.true.
         else 
            halogen=.false.
         end if
      else if ((at(i2).eq.17) .or. (at(i2).eq.35) .or. (at(i2).eq.53) .or. & 
               &  (at(i2) .eq. 85)) then
         if (at(i1) .ne. 1) then
            at_H=i2
            at_A=i1
            halogen=.true.
         else 
            halogen=.false.
         end if
      else 
         halogen=.false.
      end if 
!
!    Look for donor atoms (N or O) in all atoms of different molecules
! 
!    for all donor atoms: calculate the halogen bond energy/gradient
!
      if (halogen) then
         do j=1,n
            if (molnum(at_H) .ne. molnum(j)) then
               if ((at(j) .eq. 7) .or. (at(j) .eq. 8)) then          
!
!    Check if the atoms are close enough to each other
!  
                  rbond(1)=xyz(1,at_A)-xyz(1,at_H)
                  rbond(2)=xyz(2,at_A)-xyz(2,at_H)
                  rbond(3)=xyz(3,at_A)-xyz(3,at_H)
!
!     apply periodic boundaries, if needed
!
                  if (periodic) then
                     call box_image(rbond)
                  end if
             
                  ri=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)

                  rbond(1)=xyz(1,j)-xyz(1,at_H)
                  rbond(2)=xyz(2,j)-xyz(2,at_H)
                  rbond(3)=xyz(3,j)-xyz(3,at_H)
!
!     apply periodic boundaries, if needed
!
                  if (periodic) then
                     call box_image(rbond)
                  end if

                  rj=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)


                  dum1=1.2*(rad(at(at_A))+rad(at(at_H)))/0.52917726
                  dum2=1.2*(rad(at(j))+rad(at(at_H)))/0.52917726
                  if (ri.lt.dum1 .or. rj.lt.dum2) then
                     r=sqrt((xyz(1,at_A)-xyz(1,j))**2 &
                       &   +(xyz(2,at_A)-xyz(2,j))**2 &
                       &   +(xyz(3,at_A)-xyz(3,j))**2)
                     if (r.gt.15.d0) cycle
                     call hbpara(-6.5d0,1.0d0,q_glob(at_H),dum1)  !OK
                     c1=scalexb_glob(at(at_H))*dum1  !OK
                     call eabxag(n,at_A,j,at_H,xyz,c1,eh,g)  !OK
                  end if
               end if
            end if 
         end do 
      end if
      
!
!    Detect a Hydrogen bond and assign the roles of the hydrogen and the donor
!
!    Case 1) Atom i1 is hydrogen and atom i2 is the acceptor
!
      hydrogen=.false.
      if (at(i1).eq.1) then
         if ((at(i2).eq.7) .or. (at(i2).eq.8) .or. (at(i2).eq.9) .or. (at(i2).eq.16) &
            &  .or. (at(i2).eq.17)) then   
            at_H=i1
            at_A=i2
            hydrogen=.true.
         end if
      end if
!
!    Case 2) Atom i2 is hydrogen and atom i1 is the acceptor
!
      if (at(i2).eq.1) then
         if ((at(i1).eq.7) .or. (at(i1).eq.8) .or. (at(i1).eq.9) .or. (at(i1).eq.16) & 
            &  .or. (at(i1).eq.17)) then
            at_H=i2
            at_A=i1
            hydrogen=.true.
         end if
      end if

!
!    If a hydrogenbond-able bond was detected, look for donor atoms in the 
!    other molecules 
!
      if (hydrogen) then
         do j=1,n
            if (molnum(at_H) .ne. molnum(j)) then
!
!    Do a prescreening with the parameter array
!
               cpar=scalehb_glob(at(at_A))*scalehb_glob(at(j))
               if (cpar .gt. 1D-6) then
!
!    Check if the atoms are close enough to each other
!

                  rbond(1)=xyz(1,at_A)-xyz(1,at_H)
                  rbond(2)=xyz(2,at_A)-xyz(2,at_H)
                  rbond(3)=xyz(3,at_A)-xyz(3,at_H)
!
!     apply periodic boundaries, if needed
!
                  if (periodic) then
                     call box_image(rbond)
                  end if

                  ri=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)

                  rbond(1)=xyz(1,j)-xyz(1,at_H)
                  rbond(2)=xyz(2,j)-xyz(2,at_H)
                  rbond(3)=xyz(3,j)-xyz(3,at_H)
!
!     apply periodic boundaries, if needed
!
                  if (periodic) then
                     call box_image(rbond)
                  end if

                  rj=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)

                  dum1=1.3*(rad(at(at_A))+rad(1))/0.52917726
                  dum2=1.3*(rad(at(j))+rad(1))/0.52917726
                  if (ri.lt.dum1.or.rj.lt.dum2)then
                     rbond(1)=xyz(1,at_A)-xyz(1,j)
                     rbond(2)=xyz(2,at_A)-xyz(2,j)
                     rbond(3)=xyz(3,at_A)-xyz(3,j)
!
!     apply periodic boundaries, if needed
!

                     if (periodic) then
                        call box_image(rbond)
                     end if

                     r=sqrt(rbond(1)**2+rbond(2)**2+rbond(3)**2)
                     
                     if (r.gt.15.d0) cycle
                     call hbpara(10.0d0,5.0d0,q_glob(at_A),dum1)
                     c2=dum1*scalehb_glob(at(at_A))
                     call hbpara(10.0d0,5.0d0,q_glob(j),dum1)
                     c1=dum1*scalehb_glob(at(j))
                     call eabhag(n,j,at_A,at_H,xyz,c1,c2,eh,g)
                  end if
               end if
            end if
         end do

      end if

   end do

end if


return
end subroutine ff_hb
