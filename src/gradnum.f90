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
!     subroutine gradnum: calculates the potential EVB-QMDFF-energy
!     and numerical first derivatives with respect to Cartesian coordinates
!     Mainly for test-purposes (comparison with analytical gradient)
!
!     part of EVB
!
subroutine gradnum (xyz2,g_evb,deltaQ,deltaQ2)
use general
use evb_mod
implicit none
integer i,j,n2,k
real(kind=8)::xyz2(3,natoms),e_evb,g_evb(3,natoms) 
real(kind=8)::e1_shifted,e2_shifted,e3_shifted,e_two,gnorm_two
real(kind=8)::e_three,gnorm_three,gnorm_evb,step
real(kind=8)::ediff,offdiag,off4,root2,deldiscr,delsqrt
real(kind=8)::e,gnorm,e_qmdff1,e_qmdff2,e_qmdff3,deltaQ,deltaQ2
real(kind=8)::dE12,dE13,dE23,e_lower,e_upper
real(kind=8),dimension(3*natoms)::dq
real(kind=8),dimension(:,:),allocatable::U,Ut,g_mat,UtA
!
!     variables for lapack-dsyev-Subroutine
!
character(len=1)::JOBZ,UPLO
integer::INFO,LDA,LWORK,Nd
real(kind=8),dimension(:,:),allocatable::Aa,A1
real(kind=8),dimension(:),allocatable::W,WORK
real(kind=8),dimension(3,3)::evb_mat

JOBZ='N' !only eigenvalues
UPLO='U' !upper triangle of a
Nd=3
LDA=Nd
INFO=0
LWORK=Nd*Nd-1
allocate(A(Nd,Nd))
allocate(W(Nd))
allocate(WORK(LWORK))
if (nqmdff.eq.2) then
   step=num_grad_step
   do i=1,natoms
      do j=1,3
         do k=1,2
            if (k.eq.1) then
               xyz2(j,i)=xyz2(j,i)+step
            else if (k.eq.2) then
               xyz2(j,i)=xyz2(j,i)-2*step
            end if
            call ff_eg(n_one,at,xyz2,e,g_one)
            call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e,g_one)
            call ff_hb(n_one,at,xyz2,e,g_one)
            call ff_eg_two(n_one,at,xyz2,e_two,g_two)
            call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,r094_mod,sr42,&
                &c6xy_two,e_two,g_two)
            call ff_hb_two(n_one,at,xyz2,e_two,g_two)
            e1_shifted=e+E_zero1
            e2_shifted=e_two+E_zero2
            e_qmdff1=e1_shifted
            e_qmdff2=e2_shifted
            !  Maybe this is needed for EVB optimizsation: abs(...)
            ediff=e1_shifted-e2_shifted
            if (use_dq .eqv. .true.) then
               if (off_basis=="const") then
                  offdiag=offa
               else if (off_basis=="1g") then
                  offdiag=offa*exp(-offb*deltaQ*deltaQ)  !a*exp(-b*deltaE^2)
               else if (off_basis=="3g") then
                  offdiag=offa*exp(-offb*deltaQ*deltaQ)+offc*exp(-offd*(deltaQ+offm)&
                         &*(deltaQ+offm))+offe*exp(-offf*(deltaQ+offn)*(deltaQ+offn))
               else if (off_basis=="sd2") then
                  offdiag=offa*exp(-offb*deltaQ*deltaQ)+offc*deltaQ*deltaQ*&
                         &exp(-offd*deltaQ*deltaQ)+offe*deltaQ*deltaQ*&
                         &exp(-offf*deltaQ*deltaQ)
               end if
            else 
                if (off_basis=="const") then
                   offdiag=offa      !a
                else if (off_basis=="1g") then
                   offdiag=offa*exp(-offb*ediff*ediff)  !a*exp(-b*deltaE^2)
                else if (off_basis=="2g") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*exp(-offd*(ediff+offm)&
                           &*(ediff+offm))
                else if (off_basis=="3g") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*exp(-offd*(ediff+offm)&
                           &*(ediff+offm))+offe*exp(-offf*(ediff+offn)*(ediff+offn))
                else if (off_basis=="sp") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*exp(-offd*ediff*ediff)
                else if (off_basis=="sd") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*ediff*exp(-offd*ediff*ediff)
                else if (off_basis=="sd2") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*ediff*&
                           &exp(-offd*ediff*ediff)+offe*ediff*ediff*exp(-offf*ediff*ediff)
                else if (off_basis=="sp2d") then
                   offdiag=offa*exp(-offb*ediff*ediff)+offc*ediff*&
                           &exp(-offd*ediff*ediff)+offe*ediff*exp(-offf*ediff*ediff)+&
                           &ediff*ediff*offg*exp(-offh*ediff*ediff)
                end if
            end if
            off4=4.d0*offdiag*offdiag
            root2=sqrt(ediff*ediff+off4)
            e_evb=0.5d0*(e1_shifted+e2_shifted-root2)
            if (k.eq.1) then
               e_upper=E_evb
            else if (k.eq.2) then
               e_lower=E_evb
            end if
         end do
         g_evb(j,i)=(e_upper-e_lower)/(2*step)
         g_evb(j,i)=g_evb(j,i)*hartree*bohr
         xyz2(j,i)=xyz2(j,i)+step
      end do
   end do
   gnorm_evb=sqrt(sum(g_evb**2))
else if (nqmdff.eq.3) then
   step=num_grad_step
   do i=1,natoms
      do j=1,3
         do k=1,2
            if (k.eq.1) then
               xyz2(j,i)=xyz2(j,i)+step
            else if (k.eq.2) then
               xyz2(j,i)=xyz2(j,i)-2*step
            end if
            call ff_eg(n_one,at,xyz2,e,g_one)
            call ff_nonb(n_one,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e,g_one)
            call ff_hb(n_one,at,xyz2,e,g_one)
 
            call ff_eg_two(n_one,at,xyz2,e_two,g_two)
            call ff_nonb_two(n_one,at,xyz2,q_two,r0ab,zab,r094_mod,sr42,&
                    &c6xy_two,e_two,g_two)
            call ff_hb_two(n_one,at,xyz2,e_two,g_two)
  
            call ff_eg_three(n_one,at,xyz2,e_three,g_three)
            call ff_nonb_three(n_one,at,xyz2,q_three,r0ab,zab,r094_mod,sr42,&
                    &c6xy_three,e_three,g_three)
            call ff_hb_three(n_one,at,xyz2,e_three,g_three)

            e1_shifted=e+E_zero1
            e2_shifted=e_two+E_zero2
            e3_shifted=e_three+E_zero3
            e_qmdff1=e1_shifted
            e_qmdff2=e2_shifted
            e_qmdff3=e3_shifted
            evb_mat(1,1)=e_qmdff1
            evb_mat(2,2)=e_qmdff2
            evb_mat(3,3)=e_qmdff3
            if (use_dq) then
               if (off_basis=="1g") then
                  evb_mat(1,2)=offa*exp(-offb*deltaQ*deltaQ)
                  evb_mat(2,1)=evb_mat(1,2)
                  evb_mat(2,3)=offc*exp(-offd*deltaQ2*deltaQ2)
                  evb_mat(3,2)=evb_mat(2,3)
               else if (off_basis=="sd2") then
                  evb_mat(1,2)=offa*exp(-offb*deltaQ*deltaQ)+offc*deltaQ*deltaQ*&
                     &exp(-offd*deltaQ*deltaQ)+offe*deltaQ*deltaQ*exp(-offf*deltaQ*deltaQ)
                  evb_mat(2,1)=evb_mat(1,2)
                  evb_mat(2,3)=offg*exp(-offh*deltaQ2*deltaQ2)+offi*deltaQ2*deltaQ2*&
                     &exp(-offj*deltaQ2*deltaQ2)+offk*deltaQ2*deltaQ2*exp(-offl*&
                     &deltaQ2*deltaQ2)
                  evb_mat(3,2)=evb_mat(2,3)
               end if
                  evb_mat(1,3)=0
                  evb_mat(3,1)=0
            else
               dE12=abs(e_qmdff1-e_qmdff2)
               dE13=abs(e_qmdff1-e_qmdff3)
               dE23=abs(e_qmdff2-e_qmdff3)
               if (off_basis=="1g") then  
                  evb_mat(1,2)=offa*exp(-offb*dE12*dE12)
                  evb_mat(2,1)=evb_mat(1,2)
                  evb_mat(2,3)=offc*exp(-offd*dE23*dE23)
                  evb_mat(3,2)=evb_mat(2,3)
                  if (full .eqv. .true.) then
                     evb_mat(1,3)=offe*exp(-offf*dE13*dE13)
                     evb_mat(3,1)=evb_mat(1,3)
                  else
                     evb_mat(1,3)=0
                     evb_mat(3,1)=0
                  end if
               else if (off_basis=="sd2") then
                  evb_mat(1,2)=offa*exp(-offb*dE12*dE12)+offc*dE12*dE12*&
                      &exp(-offd*dE12*dE12)+offe*dE12*dE12*exp(-offf*dE12*dE12)
                  evb_mat(2,1)=evb_mat(1,2)
                  evb_mat(2,3)=offg*exp(-offh*dE23*dE23)+offi*dE23*dE23*&
                      &exp(-offj*dE12*dE23)+offk*dE23*dE23*exp(-offl*dE23*dE23)
                  evb_mat(3,2)=evb_mat(2,3)
                  if (full .eqv. .true.) then
                     evb_mat(1,3)=offm*exp(-offn*dE13*dE13)
                     evb_mat(3,1)=evb_mat(1,3)
                  else
                     evb_mat(1,3)=0
                     evb_mat(3,1)=0
                  end if
               end if
            end if
            Aa=evb_mat
            call DSYEV(JOBZ,UPLO,Nd,Aa,LDA,W,WORK,LWORK,INFO)
            E_evb=W(1)
            if (k.eq.1) then
               e_upper=E_evb
            else if (k.eq.2) then
               e_lower=E_evb
            end if
         end do
         g_evb(j,i)=(e_upper-e_lower)/(2*step)
         g_evb(j,i)=g_evb(j,i)*hartree*bohr
         xyz2(j,i)=xyz2(j,i)+step
      end do
   end do
   gnorm_evb=sqrt(sum(g_evb**2))
end if
return
end subroutine gradnum
