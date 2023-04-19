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
!     This subroutine calculates the fitness for all dE coupling terms
!     in a Levenberg Marquardt optimization
!  
!     part of EVB
!
subroutine lm_de_func(m_ind,n_ind,x_var,fvec,iflag)
use evb_mod
use lm_module
implicit none
integer::m_ind,n_ind
real(kind=8)::fvec(m_ind),fun
real(kind=8)::ep1,ep2,ep3,ep4,ep5,ep6,denominator
real(kind=8)::E1,E2,V12,deltaE
integer::i,iflag
real(kind=8)::x_var(n_ind),xdat(m_ind),ydat(m_ind)


!
!     The coupling parameters are now stored into the x(:) array
!     n=number of parameters, m=number of points 
!
! 
!    store the parameters in single real numbers
! 
offa=x_var(1)
offb=x_var(2)
if (off_basis .eq. "2g") then
   offc=x_var(3)
   offd=x_var(4)
   offm=x_var(5)
else if (off_basis .eq. "3g") then
   offc=x_var(3)
   offd=x_var(4)
   offe=x_var(5)
   offf=x_var(6)
   offm=x_var(7)
   offn=x_var(8)
else if (off_basis .eq. "sd" .or. off_basis .eq. "sp") then
   offc=x_var(3)
   offd=x_var(4)
else if (off_basis .eq. "sd2") then
   offc=x_var(3)
   offd=x_var(4)
   offe=x_var(5)
   offf=x_var(6)
else if (off_basis .eq. "sp2d") then
   offc=x_var(3)
   offd=x_var(4)
   offe=x_var(5)
   offf=x_var(6)
   offg=x_var(7)
   offh=x_var(8)
else if (off_basis .eq. "sp2d3") then
   offc=x_var(3)
   offd=x_var(4)
   offe=x_var(5)
   offf=x_var(6)
   offg=x_var(7)
   offh=x_var(8)
   offi=x_var(9)
   offj=x_var(10)
   offk=x_var(11)
   offl=x_var(12)
end if

do i=1,m_ind
   E1=ff_e1_lm(i)
   E2=ff_e2_lm(i)
   deltaE=abs(E1-E2)
   if (off_basis .eq. "1g") then
      ep1=exp(-offb*deltaE*deltaE)
      denominator=sqrt(deltaE*deltaE+4*offa*offa*ep1*ep1)
      fvec(i)=((energies_ref_lm(i)-0.5d0*((E1+E2-denominator)))*&
         (energies_ref_lm(i)-0.5d0*((E1+E2-denominator))))
   else if (off_basis .eq. "2g") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*(deltaE+offm)*(deltaE+offm))
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*&
           &deltaE+4*(offa*ep1+offc*ep2)*(offa*ep1+offc*ep2)))))*&
           &((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*&
           &deltaE+4*(offa*ep1+offc*ep2)*(offa*ep1+offc*ep2))))))
   else if (off_basis .eq. "3g") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*(deltaE+offm)*(deltaE+offm))
      ep3=exp(-offf*(deltaE+offn)*(deltaE+offn))
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*&
           &deltaE+4*(offa*ep1+offc*ep2+offe*ep3)*(offa*ep1+offc*&
           &ep2+offe*ep3)))))*((energies_ref_lm(i)-0.5*(E1+&
           &E2-sqrt(deltaE*deltaE+4*(offa*ep1+offc*ep2+offe*ep3)*&
           &(offa*ep1+offc*ep2+offe*ep3))))))
   else if (off_basis .eq. "sd") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*(deltaE*deltaE))
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+&
           &4*(offa*ep1+offc*deltaE*deltaE*ep2)*(offa*ep1+offc*deltaE*&
           &deltaE*ep2)))))*((energies_ref_lm(i)-0.5*(E1+E2-&
           &sqrt(deltaE*deltaE+4*(offa*ep1+offc*deltaE*deltaE*ep2)*&
           &(offa*ep1+offc*deltaE*deltaE*ep2))))))
   else if (off_basis .eq. "sp") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*(deltaE*deltaE))
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+&
           &4*(offa*ep1+offc*deltaE*ep2)*(offa*ep1+offc*deltaE*ep2)))))*&
           &((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*&
           &deltaE+4*(offa*ep1+offc*deltaE*ep2)*(offa*ep1+offc*deltaE*&
           &ep2))))))
   else if (off_basis .eq. "sd2") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*deltaE*deltaE)
      ep3=exp(-offf*deltaE*deltaE)
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+4*&
           &(offa*ep1+offc*deltaE*deltaE*ep2+offe*deltaE*deltaE*ep3)*&
           &(offa*ep1+offc*deltaE*deltaE*ep2+offe*deltaE*deltaE*ep3))))))*&
           &(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*&
           &deltaE+4*(offa*ep1+offc*deltaE*deltaE*ep2+offe*deltaE*deltaE*&
           &ep3)*(offa*ep1+offc*deltaE*deltaE*ep2+offe*deltaE*deltaE*&
           &ep3))))))
   else if (off_basis .eq. "sp2d") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*deltaE*deltaE)
      ep3=exp(-offf*deltaE*deltaE)
      ep4=exp(-offh*deltaE*deltaE)
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+&
           &4*(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*&
           &deltaE*ep4)*(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*&
           &deltaE*deltaE*ep4)))))*&
           &(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+&
           &4*(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*&
           &deltaE*ep4)*(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*&
           &deltaE*deltaE*ep4)))))))
   else if (off_basis .eq. "sp2d3") then
      ep1=exp(-offb*deltaE*deltaE)
      ep2=exp(-offd*deltaE*deltaE)
      ep3=exp(-offf*deltaE*deltaE)
      ep4=exp(-offh*deltaE*deltaE)
      ep5=exp(-offj*deltaE*deltaE)
      ep6=exp(-offl*deltaE*deltaE)
      fvec(i)=(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+&
           &4*(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*&
           &deltaE*ep4+offi*deltaE*deltaE*ep5+offk*deltaE*deltaE*ep6)*&
           &(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*&
           &deltaE*ep4+offi*deltaE*deltaE*ep5+offk*deltaE*deltaE*ep6))))))*&
           &(((energies_ref_lm(i)-0.5*(E1+E2-sqrt(deltaE*deltaE+4* &
           &(offa*ep1+offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*deltaE*&
           &ep4+offi*deltaE*deltaE*ep5+offk*deltaE*deltaE*ep6)*(offa*ep1+&
           &offc*deltaE*ep2+offe*deltaE*ep3+offg*deltaE*deltaE*ep4+offi*&
           &deltaE*deltaE*ep5+offk*deltaE*deltaE*ep6))))))
   end if
end do

return
end subroutine lm_de_func

