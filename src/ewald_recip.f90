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
!     ###################################################
!     ##  COPYRIGHT (C)  1999  by  Jay William Ponder  ##
!     ##              All Rights Reserved              ##
!     ###################################################
!
!     subroutine ewald_recip: Calculate the reciprocal part of the 
!       particle-Mesh Ewald method via discrete Fourier transform
!
!     part of QMDFF
!
subroutine  ewald_recip (n,xyz,q,energy,grad)
use qmdff
implicit none
integer::i,j,k,m  ! loop indices
integer::inorm,n
real(kind=8)::xyz(3,n),q(n)
real(kind=8)::recip(3,3) ! reciprocal lattice vectors as matrix columns
!  variables for finding of PME B-spline coefficients
integer::ifr
real(kind=8)::xi,yi,zi
real(kind=8)::w,fr,eps
real(kind=8),allocatable::igrid(:,:)  ! interpolation grid to fill
! ! B-spline coefficients along the a/b/c axes
real(kind=8),allocatable::thetai1(:,:,:),thetai2(:,:,:),thetai3(:,:,:)  
!   For filling of array with contributions from electrostatic sites
integer::cid(3),nearpt(3),abound(6),cbound(6)
integer::nlpts,nrpts,grdoff
integer::nchk1,nchk2,nchk3
integer::ngrd1,ngrd2,ngrd3
logical::negx,negy,negz
logical::posx,posy,posz
logical::midx,midy,midz
integer::nchunk

integer::pmetable(n,1)  ! PME grid spatial regions involved for each site
!   For placing of atomic charges onto PME grid
integer::ii,jj,kk
integer::ichk,isite,iatm
integer::offsetx,offsety,offsetz
real(kind=8)::v0,u0,t0
real(kind=8)::term
real(kind=8)::de1,de2,de3
real(kind=8)::denom
real(kind=8)::dn,dt1,dt2,dt3,e
real(kind=8)::energy,grad(3,n)
real(kind=8)::expterm,h1,h2,h3,hsq
real(kind=8)::fi,f
real(kind=8)::dielec
integer::npoint
integer::i0,igrd0,it1,it2,it3
real(kind=8)::ar1,ar2,ar3,br1,br2,br3,cr1,cr2,cr3
integer::j0,jgrd0,k0,k1,k2,k3,kgrd0
integer::m1,m2,m3,nf,nff
real(kind=8)::r1,r2,r3,t1,t2,t3
real(kind=8)::struc2,volterm,vterm
real(kind=8)::pterm
real(kind=8)::bohr



bohr=0.52917721092d0


!
!     compute and store reciprocal lattice vectors as column
!
ar1=box_len
ar2=0.d0
ar3=0.d0
br1=0.d0
br2=box_len
br3=0.d0
cr1=0.d0
cr2=0.d0
cr3=box_len

recip(1,1) = (br2*cr3 - cr2*br3) / volbox
recip(2,1) = (br3*cr1 - cr3*br1) / volbox
recip(3,1) = (br1*cr2 - cr1*br2) / volbox
recip(1,2) = (cr2*ar3 - ar2*cr3) / volbox
recip(2,2) = (cr3*ar1 - ar3*cr1) / volbox
recip(3,2) = (cr1*ar2 - ar1*cr2) / volbox
recip(1,3) = (ar2*br3 - br2*ar3) / volbox
recip(2,3) = (ar3*br1 - br3*ar1) / volbox
recip(3,3) = (ar1*br2 - br1*ar2) / volbox


!
!     zero out the particle mesh Ewald charge grid
!

do k = 1, nfft
   do j = 1, nfft
      do i = 1, nfft
         qgrid(1,i,j,k) = 0.0d0
         qgrid(2,i,j,k) = 0.0d0
      end do
   end do
end do

!
!     perform dynamic allocation of some global arrays
!
!if (first) then
!   first = .false.
   if (.not. allocated(igrid))  allocate (igrid(3,n))
!end if
!
!     offset used to shift sites off exact lattice bounds
!
eps = 1.0d-8
allocate(thetai1(4,bsorder,n))
allocate(thetai2(4,bsorder,n))
allocate(thetai3(4,bsorder,n))

!
!     get the B-spline coefficients for each atomic site
!
 do i = 1, n
    xi = xyz(1,i)
    yi = xyz(2,i)
    zi = xyz(3,i)
    w = xi*recip(1,1) + yi*recip(2,1) + zi*recip(3,1)
    fr = dble(nfft) * (w-dble(anint(w))+0.5d0)
    ifr = int(fr-eps)
    w = fr - dble(ifr)
    igrid(1,i) = ifr - bsorder
    call bsplgen (w,thetai1(1,1,i),bsorder)
    w = xi*recip(1,2) + yi*recip(2,2) + zi*recip(3,2)
    fr = dble(nfft) * (w-dble(anint(w))+0.5d0)
    ifr = int(fr-eps)
    w = fr - dble(ifr)
    igrid(2,i) = ifr - bsorder
    call bsplgen (w,thetai2(1,1,i),bsorder)
    w = xi*recip(1,3) + yi*recip(2,3) + zi*recip(3,3)
    fr = dble(nfft) * (w-dble(anint(w))+0.5d0)
    ifr = int(fr-eps)
    w = fr - dble(ifr)
    igrid(3,i) = ifr - bsorder
    call bsplgen (w,thetai3(1,1,i),bsorder)
end do


!
!      Set some parameters regarding the spline generation     
!
nchunk = 1
nchk1 = 1
nchk2 = 1
nchk3 = 1
nlpts =  2
nrpts =  2
grdoff = 4
ngrd1=nfft
ngrd2=nfft
ngrd3=nfft

!
!     Construct an array which stores the spatial regions of the 
!     particle mesh Ewald grid with contributions from each 
!     electrostatic site
!
    
!
!     zero out the PME table marking chunks per site
!
do k = 1, nchunk
   do i = 1, n
      pmetable(i,k) = 0
   end do
end do
!
!     loop over sites to find the spatial chunks for each
!
do i = 1, n
   nearpt(1) = igrid(1,i) + grdoff
   nearpt(2) = igrid(2,i) + grdoff
   nearpt(3) = igrid(3,i) + grdoff
   if (nearpt(1) .lt. 1) then
      nearpt(1) = mod(nearpt(1),nfft) + nfft
   else if (nearpt(1) .gt. nfft) then
      nearpt(1) = mod(nearpt(1),nfft)
   end if
   if (nearpt(2) .lt. 1) then
      nearpt(2) = mod(nearpt(2),nfft) + nfft
   else if (nearpt(2) .gt. nfft) then
      nearpt(2) = mod(nearpt(2),nfft)
   end if
   if (nearpt(3) .lt. 1) then
      nearpt(3) = mod(nearpt(3),nfft) + nfft
   else if (nearpt(3) .gt. nfft) then
      nearpt(3) = mod(nearpt(3),nfft)
   end if
   abound(1) = nearpt(1) - nlpts
   abound(2) = nearpt(1) + nrpts
   abound(3) = nearpt(2) - nlpts
   abound(4) = nearpt(2) + nrpts
   abound(5) = nearpt(3) - nlpts
   abound(6) = nearpt(3) + nrpts
   cid(1) = (nearpt(1)-1)/ngrd1 + 1
   cid(2) = (nearpt(2)-1)/ngrd2 + 1
   cid(3) = (nearpt(3)-1)/ngrd3 + 1
   cbound(1) = (cid(1)-1)*ngrd1 + 1
   cbound(2) = cbound(1) + ngrd1 - 1
   cbound(3) = (cid(2)-1)*ngrd2 + 1
   cbound(4) = cbound(3) + ngrd2 - 1
   cbound(5) = (cid(3)-1)*ngrd3 + 1
   cbound(6) = cbound(5) + ngrd3 - 1
!
!     set and store central chunk where the site is located
!
   k = 1 !(cid(3)-1)*nchk1*nchk2 + (cid(2)-1)*nchk1 + cid(1)
   pmetable(i,k) = 1
!
!     flags for atom bounds to left or right of central chunk
!
   negx = (abound(1) .lt. cbound(1))
   negy = (abound(3) .lt. cbound(3))
   negz = (abound(5) .lt. cbound(5))
   posx = (abound(2) .gt. cbound(2))
   posy = (abound(4) .gt. cbound(4))
   posz = (abound(6) .gt. cbound(6))
!
!     flags for atom bounds fully inside the central chunk
!
   midx = (.not.negx .and. .not.posx)
   midy = (.not.negy .and. .not.posy)
   midz = (.not.negz .and. .not.posz)
   if (midx .and. midy .and. midz)  cycle
!
!     flags for atom bounds that overlap the central chunk
!
   midx = (.not.negx .or. .not.posx)
   midy = (.not.negy .or. .not.posy)
   midz = (.not.negz .or. .not.posz)
!
!     check for overlap with any of the neighboring chunks
!
   if (midx .and. midy .and. negz)  call setchunk (i,cid,0,0,-1,pmetable)
   if (midx .and. midy .and. posz)  call setchunk (i,cid,0,0,1,pmetable)
   if (midx .and. negy .and. midz)  call setchunk (i,cid,0,-1,0,pmetable)
   if (midx .and. posy .and. midz)  call setchunk (i,cid,0,1,0,pmetable)
   if (negx .and. midy .and. midz)  call setchunk (i,cid,-1,0,0,pmetable)
   if (posx .and. midy .and. midz)  call setchunk (i,cid,1,0,0,pmetable)
   if (midx .and. negy .and. negz)  call setchunk (i,cid,0,-1,-1,pmetable)
   if (midx .and. negy .and. posz)  call setchunk (i,cid,0,-1,1,pmetable)
   if (midx .and. posy .and. negz)  call setchunk (i,cid,0,1,-1,pmetable)
   if (midx .and. posy .and. posz)  call setchunk (i,cid,0,1,1,pmetable)
   if (negx .and. midy .and. negz)  call setchunk (i,cid,-1,0,-1,pmetable)
   if (negx .and. midy .and. posz)  call setchunk (i,cid,-1,0,1,pmetable)
   if (posx .and. midy .and. negz)  call setchunk (i,cid,1,0,-1,pmetable)
   if (posx .and. midy .and. posz)  call setchunk (i,cid,1,0,1,pmetable)
   if (negx .and. negy .and. midz)  call setchunk (i,cid,-1,-1,0,pmetable)
   if (negx .and. posy .and. midz)  call setchunk (i,cid,-1,1,0,pmetable)
   if (posx .and. negy .and. midz)  call setchunk (i,cid,1,-1,0,pmetable)
   if (posx .and. posy .and. midz)  call setchunk (i,cid,1,1,0,pmetable)
   if (negx .and. negy .and. negz)  call setchunk (i,cid,-1,-1,-1,pmetable)
   if (negx .and. negy .and. posz)  call setchunk (i,cid,-1,-1,1,pmetable)
   if (negx .and. posy .and. negz)  call setchunk (i,cid,-1,1,-1,pmetable)
   if (posx .and. negy .and. negz)  call setchunk (i,cid,1,-1,-1,pmetable)
   if (negx .and. posy .and. posz)  call setchunk (i,cid,-1,1,1,pmetable)
   if (posx .and. negy .and. posz)  call setchunk (i,cid,1,-1,1,pmetable)
   if (posx .and. posy .and. negz)  call setchunk (i,cid,1,1,-1,pmetable)
   if (posx .and. posy .and. posz)  call setchunk (i,cid,1,1,1,pmetable)
end do

!
!     Place the atomic charges onto the particle mesh Ewald grid
!

do k = 1, nfft
   do j = 1, nfft
      do i = 1, nfft
         qgrid(1,i,j,k) = 0.0d0
         qgrid(2,i,j,k) = 0.0d0
      end do
   end do
end do

!
!     put the permanent multipole moments onto the grid
!
offsetx=0
offsety=0
offsetz=0
do ichk = 1, nchunk
   cid(1) = mod(ichk-1,nchk1)
   cid(2) = mod(((ichk-1-cid(1))/nchk1),nchk2)
   cid(3) = mod((ichk-1)/(nchk1*nchk2),nchk3)
   cbound(1) = cid(1)*ngrd1 + 1
   cbound(2) = cbound(1) + ngrd1 - 1
   cbound(3) = cid(2)*ngrd2 + 1
   cbound(4) = cbound(3) + ngrd2 - 1
   cbound(5) = cid(3)*ngrd3 + 1
   cbound(6) = cbound(5) + ngrd3 - 1
   do iatm = 1, n
      if (pmetable(iatm,ichk) .eq. 1) then
         nearpt(1) = igrid(1,iatm) + grdoff
         nearpt(2) = igrid(2,iatm) + grdoff
         nearpt(3) = igrid(3,iatm) + grdoff
         abound(1) = nearpt(1) - nlpts
         abound(2) = nearpt(1) + nrpts
         abound(3) = nearpt(2) - nlpts
         abound(4) = nearpt(2) + nrpts
         abound(5) = nearpt(3) - nlpts
         abound(6) = nearpt(3) + nrpts
         call ewald_adjust (offsetx,nfft,nchk1,abound(1), &
              &       abound(2),cbound(1),cbound(2)) 
         call ewald_adjust (offsety,nfft,nchk2,abound(3), &
              &       abound(4),cbound(3),cbound(4))
         call ewald_adjust (offsetz,nfft,nchk3,abound(5), & 
              &       abound(6),cbound(5),cbound(6))
         do kk = abound(5), abound(6)
            k = kk
            m = k + offsetz
            if (k .lt. 1)  k = k + nfft
            v0 = thetai3(1,m,iatm) * q(iatm)
            do jj = abound(3), abound(4)
               j = jj
               m = j + offsety
               if (j .lt. 1)  j = j + nfft
               u0 = thetai2(1,m,iatm)
               term = v0 * u0
               do ii = abound(1), abound(2)
                  i = ii
                  m = i + offsetx
                  if (i .lt. 1)  i = i + nfft
                  t0 = thetai1(1,m,iatm)
                  qgrid(1,i,j,k) = qgrid(1,i,j,k) + term*t0
               end do
            end do
         end do
      end if
   end do
end do

!
!     Perform the 3-D FFT forward transform using FFTW
!
call dfftw_execute_dft (planf,qgrid,qgrid)

!
!     use scalar sum to get reciprocal space energy and virial
!
f = 0.5d0
npoint = nfft**3
pterm = (pi/a_ewald)**2
energy=0.d0
volterm = pi * volbox
nff = nfft * nfft
nf = (nfft+1) / 2
do i = 1, npoint-1
   k3 = i/nff + 1
   j = i - (k3-1)*nff
   k2 = j/nfft + 1
   k1 = j - (k2-1)*nfft + 1
   m1 = k1 - 1
   m2 = k2 - 1
   m3 = k3 - 1
   if (k1 .gt. nf)  m1 = m1 - nfft
   if (k2 .gt. nf)  m2 = m2 - nfft
   if (k3 .gt. nf)  m3 = m3 - nfft
   r1 = dble(m1)
   r2 = dble(m2)
   r3 = dble(m3)
   h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
   h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
   h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
   hsq = h1*h1 + h2*h2 + h3*h3
   term = -pterm * hsq
   expterm = 0.0d0
   if (term .gt. -50.0d0) then
      denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
      expterm = exp(term) / denom
      struc2 = qgrid(1,k1,k2,k3)**2 + qgrid(2,k1,k2,k3)**2
      e = f * expterm * struc2
      energy = energy + e
      vterm = (2.0d0/hsq) * (1.0d0-term) * e
   end if
   qgrid(1,k1,k2,k3) = expterm * qgrid(1,k1,k2,k3)
   qgrid(2,k1,k2,k3) = expterm * qgrid(2,k1,k2,k3)
end do

!
!     Perform the 3-D FFT backward transform using FFTW
!

call dfftw_execute_dft (planb,qgrid,qgrid)

!
!     get first derivatives of the reciprocal space energy
!
f = 1.d0
dn = dble(nfft)
grad=0.d0
do isite = 1, n
   iatm = isite
   igrd0 = igrid(1,iatm)
   jgrd0 = igrid(2,iatm)
   kgrd0 = igrid(3,iatm)
   fi = f * q(isite)
   de1 = 0.0d0
   de2 = 0.0d0
   de3 = 0.0d0
   k0 = kgrd0
   do it3 = 1, bsorder
      k0 = k0 + 1
      k = k0 + 1 + (nfft-sign(nfft,k0))/2
      t3 = thetai3(1,it3,iatm)
      dt3 = dn * thetai3(2,it3,iatm)
      j0 = jgrd0
      do it2 = 1, bsorder
         j0 = j0 + 1
         j = j0 + 1 + (nfft-sign(nfft,j0))/2
         t2 = thetai2(1,it2,iatm)
         dt2 = dn * thetai2(2,it2,iatm)
         i0 = igrd0
         do it1 = 1, bsorder
            i0 = i0 + 1
            i = i0 + 1 + (nfft-sign(nfft,i0))/2
            t1 = thetai1(1,it1,iatm)
            dt1 = dn * thetai1(2,it1,iatm)
            term = qgrid(1,i,j,k)
            de1 = de1 + term*dt1*t2*t3
            de2 = de2 + term*dt2*t1*t3
            de3 = de3 + term*dt3*t1*t2
         end do
      end do
   end do
   grad(1,iatm) = grad(1,iatm) + fi*(recip(1,1)*de1+recip(1,2)*de2 &
            &                          +recip(1,3)*de3)
   grad(2,iatm) = grad(2,iatm) + fi*(recip(2,1)*de1+recip(2,2)*de2 &
            &                          +recip(2,3)*de3) 
   grad(3,iatm) = grad(3,iatm) + fi*(recip(3,1)*de1+recip(3,2)*de2 &
            &                          +recip(3,3)*de3)
end do


end subroutine ewald_recip


