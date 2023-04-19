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
!     subroutine d0cm5:
!     computes CM5 charges in ccm5(nat) from Hirshfeld charges
!     in array chir(nat)
!     nat      : # of atoms
!     iz(nat)  : ordinal numbers
!     echo     : logical for printout
!     q(3,nat) : Cartesian coord. in Bohr
!
!     Citation of CM5:
!     Marenich, A. V.; Jerome, S. V.; Cramer, C. J.; Truhlar, D. G.
!     J. Chem. Theor. Comput. 2012, 8, 527-541.
!
!     part of QMDFF
!
subroutine docm5(nat,iz,echo,q,chir,ccm5)
IMPLICIT REAL(kind=8) (A-H,O-Z)
DIMENSION Q(3,nat),IZ(nat),CHIR(nat),CCM5(nat)
logical echo
parameter (bohr=0.52917726)

q = q * Bohr  
!
!     call original CM5 code!
!
CALL CM5MOD(NAT,IZ,CHIR,Q,CCM5,DHIRX,DHIRY,DHIRZ, &
  & DCM5X,DCM5Y,DCM5Z)
DCM5=DSQRT(DCM5X**2+DCM5Y**2+DCM5Z**2)
DHIR=DSQRT(DHIRX**2+DHIRY**2+DHIRZ**2)

if(echo) then
   WRITE (10,'(/,A,/,/,A)') &
     & ' Charges (in a.u.) from CM5PAC (June 22, 2013):', &
     & ' -----------------------------------------------'
   WRITE (10,'(A,/,A)') &
     & ' Center     Atomic      CM5         Hirshfeld', &
     & ' Number     Number      Charge      Charge'
   WRITE (10,'(A)') &
     & ' -----------------------------------------------'
     DO I=1,NAT
        WRITE (10,'(I5,6X,I5,5X,F11.6,X,F11.6)') I,IZ(I),CCM5(I),CHIR(I)
     END DO
   WRITE (10,'(A)') &
     & ' -----------------------------------------------'

   WRITE (10,'(/,A,//,A)') &
     & ' Dipole moment (in Debye)', &
     & ' -----------------------------------------------'
   WRITE (10,'(A,/,A)') &
     & '                 X        Y        Z     Total', &
     & ' -----------------------------------------------'
   WRITE (10,'(A,4F9.4)') ' CM5       ',DCM5X,DCM5Y,DCM5Z,DCM5
   WRITE (10,'(A,4F9.4)') ' Hirshfeld ',DHIRX,DHIRY,DHIRZ,DHIR
   WRITE (10,'(A,/)') &
     & ' -----------------------------------------------'
end if
q = q / Bohr

return
end subroutine docm5
