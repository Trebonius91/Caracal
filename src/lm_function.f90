!
!     subroutine lm_function: calculate Levenberg-Marquardt fitness
!     for a DG-EVB calculation
!

subroutine lm_function(m_ind,n_ind,x_var,fvec,iflag)
use evb_mod
use lm_module
implicit none
integer::m_ind,n_ind
real(kind=8)::fvec(m_ind),fun
real(kind=8)::dg_func,E1,E2,V12,deltaE,root
integer::i,iflag
real(kind=8)::x_var(n_ind),xdat(m_ind),ydat(m_ind)
!
!     the alpha values are now stored into the x(:) array 
!     n=number of alphas, m=number of path points
!
call build_dmat(dg_evb_mode_lm,mat_size_lm,x_var)
call mat_diag(mat_size_lm)

do i=1,m_ind
   V12=0
   E1=ff_e1_lm(i)
   E2=ff_e2_lm(i)
   deltaE=abs(E1-E2)
   call sum_v12(dg_evb_mode_lm,mat_size_lm,x_var,path_int_lm(:,i),V12)
   root=(0.5d0*(E1-E2))*(0.5d0*(E1-E2))+V12
   if (root .le. 0) then
      dg_func=0.5*(E1+E2)
   else
      dg_func=0.5*(E1+E2)-sqrt(root)
   end if
   fvec(i)=(dg_func-energies_ref_lm(i))*(dg_func-energies_ref_lm(i))
end do
return
end subroutine lm_function


