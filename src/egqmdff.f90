!   
!     subroutine egqmdff: for an 2-dimensional evb-qmdff-forcefield, 
!     the energies and gradients are calculated for the 2 QMDFFs
!
!     part of EVB
!
subroutine egqmdff (xyz2,e_qmdff1,e_qmdff2,g1,g2)
use evb_mod
use general
implicit none
integer::i,j,n2
real(kind=8)::xyz2(3,natoms),e_evb,g_evb(3,natoms)
real(kind=8)::g1(3,natoms),g2(3,natoms)
real(kind=8)::e1_shifted,e2_shifted,e1,gnorm_two
real(kind=8)::ediff,offdiag,off4,root2,deldiscr,delsqrt
real(kind=8)::e,gnorm,e_qmdff1,e_qmdff2,e2
n=natoms
!
!
!     First QMDFF      
!
call ff_eg(n,at,xyz2,e1,g1)
call ff_nonb(n,at,xyz2,q,r0ab,zab,r094_mod,sr42,c6xy,e1,g1)
call ff_hb(n,at,xyz2,e1,g1)
!
!     Second QMDFF
!
call ff_eg_two(n,at,xyz2,e2,g2)
call ff_nonb_two(n,at,xyz2,q_two,r0ab,zab,r094_mod,sr42, &
  &             c6xy_two,e2,g2)
call ff_hb_two(n,at,xyz2,e2,g2)
!
!     Shift the energies
!
e1_shifted=e1+E_zero1  !E(qmdff1)
e2_shifted=e2+E_zero2  !E(qmdff2)
e_qmdff1=e1_shifted
e_qmdff2=e2_shifted

return
end subroutine egqmdff

