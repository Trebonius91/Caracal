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
!     subroutine andersen: Use the simple Andersen
!       thermostat to apply a temperature on the used system
!       Small additional routines for choice of random variables 
!       are also included 
!
!     part of EVB
!
subroutine andersen
use general
use evb_mod
implicit none
integer::i,j,k
real(kind=8)::dp(natoms)
real(kind=8)::beta_n
real(kind=8)::atom_rand
real(kind=8)::atom_frac
real(kind=8)::rand ! the current variable to be randomized
!
!     scale the inverse temperature with the number of beads
!
beta_n=beta/nbeads!/nbeads 
!
!     Calculate the fraction of atoms whose velocities shall be 
!     rescaled, based on the anderson step size (inverse)
!
!atom_frac=1.d0/(andersen_step)
dp=sqrt(mass(1:natoms)/beta_n)
do i=1,3
   do j=1,natoms
!      call random(atom_rand)
!
!     If an atom was chosen to be active, rescale its velocity
!
!      if (atom_rand .le. atom_frac) then
         do k=1,nbeads
             rand=p_i(i,j,k)
             call randomn(rand)
             p_i(i,j,k)=rand*dp(j)
         end do
!      end if
   end do
end do
!stop "Guog"
return
end subroutine andersen

!
!     subroutine randomn:
!     Compute a pseudo-random number weighted by a standard normal distribution.
!     The Marsaglia polar method is used to convert from a uniform distribution to
!     a normal distribution.
!     Returns:
!     rn - The pseudo-random number

subroutine randomn(rn)
implicit none
real(kind=8)::rn
integer::iset
real(kind=8)::gset,u,v,S,fac
save iset,gset

!
!     Do the Marsaglia polar method
!
if (iset .eq. 0) then
    S = 1.d0
    do while (S .ge. 1.d0 .or. S .eq. 0.d0)
        call random(u)
        call random(v)
        u = 2.d0 * u - 1.d0
        v = 2.d0 * v - 1.d0
        S = u * u + v * v
    end do
    fac = sqrt(-2 * log(S) / S)
    gset = u * fac
    rn = v * fac
    iset = 1
else
    rn = gset
    iset = 0
end if

return
end subroutine randomn
!
!     subroutine random:
!     Compute a pseudo-random number uniformly distributed in [0,1].
!     use the included fortran routine!
!     Returns:
!     rn - The pseudo-random number
!
subroutine random(rn)
implicit none
real(kind=8):: rn
call random_number(rn)
return
end subroutine random
!
!     subroutine random_init: initializes the random number generator 
!     needed for the anderson sampling
!
subroutine random_init_local(rank)
    implicit none
    integer :: i, n, clock
    integer:: rank
    integer, dimension(:), allocatable :: seed

    call random_seed(size=n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)+rank*166763
    call random_seed(put = seed)

    deallocate(seed)

end subroutine random_init_local

!
!     Seed the pseudo-random number generator using a particular seed value.
!     Useful for replicating runs.
!
subroutine random_init_seed(value)
    implicit none
    integer, intent(in) :: value
    integer, dimension(:), allocatable :: seed
    integer :: n, i

    call random_seed(size=n)
    allocate(seed(n))

    do i = 1, n
        seed(i) = abs(value) + (i - 1)
    end do
    call random_seed(put = seed)

    deallocate(seed)

end subroutine random_init_seed


