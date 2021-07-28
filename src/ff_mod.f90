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
!     subroutine ff_mod: correct torsions that are not calculated 
!          with Extended Hückel
!
!     part of QMDFF
!
subroutine ff_mod(echo,n,xyz,q,at,nb,wbo,hyb, &
    &    cring,ringsize)
use qmdff
implicit none
integer::at(n),n,nb(20,n),hyb(n)
integer::cring(8,n),ringsize(n)
real(kind=8)::xyz(3,n),q(n),wbo(n,n)
logical::echo
integer::i,k,j,ii,jj,kk
logical::ex,samering
real(kind=8)::scal

!
!     highly coordinated bend (fitted for SF6)
!
do i=1,nangl
   jj=angl(1,i)
   if (nb(18,jj).eq.1.or.metal(at(jj)).eq.1.or.nb(20,jj).gt.4) then
      vangl(2,i)=0.5*vangl(2,i)
   end if
end do
!
!     EMPIRICAL: 
!     torsion around single bonds which have not been 
!     computed by TB  
!    
do i=1,ntors
   ii=tors(2,i)
   jj=tors(3,i)
   if (tors(6,i).eq.1) then
      scal=1.0
!
!     only truly bonded pairs
!
      if(wbo(ii,jj).lt.0.1) cycle
!
!     cases like FLP where the model barrier is (eg due to B..P 
!     interaction) very high (FLP= Frustrated Lewis Pair)  
!     
!          
      if(hyb(ii).eq.3.and.hyb(jj).eq.3) scal=0.2
!
!     interpolated
!
      if(hyb(ii).eq.2.and.hyb(jj).eq.3) scal=0.5
      if(hyb(ii).eq.3.and.hyb(jj).eq.2) scal=0.5
!
!     fitted to PAH low-lying freqs
!
      if(hyb(ii).eq.2.and.hyb(jj).eq.2) scal=1.0
!
!     not used
!     small ring torsions (vib freq of cyclopentane and chexane 
!     conformer half twist
!        if(samering(n,ii,jj,cring,ringsize)) scal=0.3       
!     torsions for metals easier (eg benzene on pt)      
!        if(metal(at(ii)).eq.1.or.metal(at(jj)).eq.1) scal=1.0

      vtors(2,i)=vtors(2,i)*scal
!     else
!        scal=1.-0.1*abs(en(at(ii))-en(at(jj)))
!        vtors(2,i)=vtors(2,i)*scal
   end if
end do

return
end subroutine ff_mod
