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
!     subroutine hess: hess: calculate the QMDFF hessian and compare it 
!           with reference data
!
!     part of QMDFF
!
subroutine hess(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy, &
       &   outval,echo,hname,h,ffonly,freq,dist, &
       &   calc,num,software,length,fname_pre)

implicit none
integer::n
integer::mode
integer::at(n)
integer::finished ! for gaussian hessian read in: no formchk here!
real(kind=8)::xyz(3,n)
real(kind=8)::q(n)
real(kind=8)::h(3*n,3*n)
real(kind=8)::freq(3*n)
real(kind=8)::outval
logical::echo,ffonly,orca,dist,calc,num
character(len=*)::hname
real(kind=8)::r0ab(94,94),zab(94,94),r094(94,94),sr42(94,94),c6xy(n,n)
real(kind=8)::ams(86)  ! maybe here should be more... (107 elements)
real(kind=8)::step,amu2au,au2cm,zpve,e,dumi,dumj,zpve2,vmad,hstep
integer::n3,i,j,k,ic,jc,ia,ja,ii,jj,info,lwork
integer::mincol,maxcol,lin,bvib,kk
real(kind=8),allocatable::h2(:,:)
real(kind=8),allocatable::u (:,:)
real(kind=8),allocatable::modesff (:,:)
real(kind=8),allocatable::modesqm (:,:)
real(kind=8),allocatable::hs(:)
real(kind=8),allocatable::freq2(:)
real(kind=8),allocatable::aux (:)
real(kind=8),allocatable::isqm(:)
real(kind=8),allocatable::gr  (:,:)
real(kind=8),allocatable::gl  (:,:)
real(kind=8),allocatable::xsave(:,:)
character(len=3)::sym   
character(len=10)::a10  
character(len=80)::fname
character(len=1)::software
!
!     for molden test output
!
integer::length
character(len=length)::fname_pre

parameter (amu2au=1.66053886E-27/9.10938215E-31)
parameter (au2cm =219474.63067d0)                      
!
!     table with atomic masses:
!
ams =(/ 1.00790d0,  4.00260d0,  6.94000d0,  9.01218d0, &
   &  10.81000d0, 12.01100d0, 14.00670d0, 15.99940d0, 18.99840d0, &
   &  20.17900d0, 22.98977d0, 24.30500d0, 26.98154d0, 28.08550d0, &
   &  30.97376d0, 32.06000d0, 35.45300d0, 39.94800d0, 39.09830d0, &
   &  40.08000d0, 44.95590d0, 47.90000d0, 50.94150d0, 51.99600d0, &
   &  54.93800d0, 55.84700d0, 58.93320d0, 58.71000d0, 63.54600d0, & 
   &  65.38000d0, 69.73500d0, 72.59000d0, 74.92160d0, 78.96000d0, &
   &  79.90400d0, 83.80000d0, 85.46780d0, 87.62000d0, 88.90590d0, &
   &  91.22000d0, 92.90640d0, 95.94000d0, 98.90620d0, 101.0700d0, &
   &  102.9055d0, 106.4000d0, 107.8680d0, 112.4100d0, 114.8200d0, &
   &  118.6900d0, 121.7500d0, 127.6000d0, 126.9045d0, 131.3000d0, &
   &  132.9054d0, 137.3300d0, &
   &  138.91d0, 140.12d0,140.91d0,144.24d0,147.00d0,150.36d0, &
   &  151.97d0, 157.25d0, 158.93d0,162.50d0,164.93d0,167.26d0, &
   &  168.93d0, 173.04d0,174.97d0, 178.4900d0, 180.9479d0, &
   &  183.8500d0, 186.2070d0, 190.2000d0, 192.2200d0, 195.0900d0, &
   &  196.9665d0, 200.5900d0, 204.3700d0, 207.2000d0, 208.9804d0, &
   &  18*0.000d0,   0.0000d0,  5*0.000d0/) 

n3=3*n
dist=.false.
!
!     get QM Hessian     
! 
if (.not.calc) then
!   if (hname.eq.'hessian'.or.hname.eq.'hessian_driver') then
!      call rdhess(n3,h,hname)
!   else if (index(fname,'.hess').ne.0) then
!      call rdohess(n3,h,hname)
!   else if (index(fname,'.out').ne.0) then
!      call rdghess(n3,h,hname)
!   end if
   if (software .eq. "O") then
      call rdohess(n*3,h,hname)
   else if (software .eq. "C") then
      call rdchess(n,n*3,h,hname)
   end if

   return
end if

allocate(hs(n3*(n3+1)/2),u(n3,n3),freq2(n3),aux(3*n3))
allocate(gr(3,n),gl(3,n),isqm(n3),h2(n3,n3),xsave(3,n))
allocate(modesff(n3,n3),modesqm(n3,n3))

if (echo) then
   write(10,*)
   write(10,*)'============='
   write(10,*)'Hessian check'
   write(10,*)'============='
   write(10,*) 
endif
!
!     define the mass weighted matrix
!
do ia = 1, n
   do ic = 1, 3
      ii = (ia-1)*3+ic
      isqm(ii)=1./sqrt(ams(at(ia)))
   end do
end do

!     if(num) then
xsave = xyz
!     due to unknown reasons, it is necessary to rotate
!     the molecule by a tiny amount to avoid some artifical
!     mixing of rotations/translations with the modes
!     the same problem appeared in old EHT versions
!     step=0.1
!     call rotmol(n,xyz,step,2.*step,3.*step)
!
!     step length
step=0.01
!
!     compute FF Hessian
!
do ia = 1, n
   do ic = 1, 3
      ii = (ia-1)*3+ic

      xyz(ic,ia)=xyz(ic,ia)+step

      call ff_eg(n,at,xyz,e,gr)
      call ff_nonb(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,e,gr)
      call ff_hb(n,at,xyz,e,gr)

      xyz(ic,ia)=xyz(ic,ia)-2.*step

      call ff_eg(n,at,xyz,e,gl)
      call ff_nonb(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,e,gl)
      call ff_hb(n,at,xyz,e,gl)

      xyz(ic,ia)=xyz(ic,ia)+step

      do ja = 1, n
         do jc = 1, 3
            jj = (ja-1)*3 + jc
            h (ii,jj) =(gr(jc,ja) - gl(jc,ja)) / (2.*step) 
         end do
      end do 
   end do 
end do 
!     done
      xyz = xsave 
!     call prmat(6,h,3*n,3*n,'H num')
!     else

!     analytical (not fully implemented yet)         
!     call ffh(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,h)
!     call prmat(6,h,3*n,3*n,'H ana')
!     endif

k=0
do i=1,n3
   do j=1,i 
      k=k+1
      hs (k)=(h(i,j)+h(j,i))*0.5
      if (abs(h(i,j)-h(j,i)).gt.1) then
         call warn('FF not fully consistent')
         write(10,*)'Hessian ',i,j,' not symmetric ',h(i,j),h(j,i)
         write(10,*)'a bend or inversion is close to linearity'
         write(10,*)'or planarity. Change the potential functions by'
         write(10,*)'the -athr <real> option. Try values'
         write(10,*)'of 2-5. If it persists call S. Grimme!'
      end if
   end do
end do

call trproj(n,n3,xyz,hs,.false.)
!
!     include masses
!
k=0
do i=1,n3
   do j=1,i 
      k=k+1
      u(j,i)=hs(k)*isqm(i)*isqm(j)
      h(j,i)=hs(k)
      h(i,j)=hs(k)
   enddo
enddo

!
!     calculate the frequencies of the QMDFF at the minimum!
!
lwork=3*n3
!
!     the diagonalization routine
!
call dsyev ('V','U',n3,u,n3,freq,aux,lwork,info)

if (echo) then
   write(*,*)
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
   write(10,*) 'vibrational frequencies from projected FF Hessian'
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
end if
!
!     save normal modes for sorting 
!     --> in g98 fake output: visible with gaussview!
!
modesff=u

zpve=0
k   =0
do i=1,n3
   freq(i)=au2cm*sign(sqrt(abs(freq(i))),freq(i))/sqrt(amu2au)
   if (freq(i).gt.1.d-6) zpve=zpve+freq(i)
   if (freq(i).lt.-10.) k=k+1                
end do

fname=fname_pre//'_molden.out'
if (ffonly) fname=fname_pre//'_molden.out2'
if (echo) then
   call preig2(6,freq,n3)
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
   write(10,*) 
   call g98fake(fname,n,at,xyz,freq,u,h2)
   if (ffonly) write(*,'(''lowest frequency = '',f10.2)')freq(1)
   if (freq(1).lt.-500) then
      call warn('FF not fully consistent')
   end if
end if
if (freq(1).lt.-10.and.ffonly) then
   call warn2('distorting structure along lowest eigenmode', &
      &     'and repeating optimization for true minimum')    
   k=0
   do i=1,n
      do j=1,3
         k=k+1
         xyz(j,i)=xyz(j,i)+0.1*u(k,1)
      end do
   end do
   dist=.true.
end if

if (k.gt.0.and.(.not.ffonly)) then
   write(10,*) 
   write(10,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(10,*)'   incorrect curvatures of FF Hessian detected'     
   write(10,*)' but dont worry because the FF PES is not exactly'
   write(10,*)'the same as the input one and re-opt may be ncessary'
   write(10,*)'in order to remove the imaginary modes. If they are '
   write(10,*)'< 200-300 cm-1 the FF seems to be ok.'
   write(10,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end if
!
!     if ffonly is set no real analysis is made!
!
if (ffonly) return
write(10,*)
write(10,*) 'reading <',trim(hname),'>'

if (echo) then
   write(10,*)
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
   write(10,*) 'vibrational frequencies from projected QM Hessian'
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
end if
!
!     read Hessian of the reference method   
!
finished=1
if(.not.ffonly) then
   if (software .eq. "O") then
      call rdohess(n3,h2,hname)
      call hpack1(n3,h2,hs)
      call trproj(n,n3,xyz,hs,.false.)
      call hpack2(n3,hs,h2)
   else if (software .eq. "C") then
      call rdchess(n,n3,h2,hname)
      call hpack1(n3,h2,hs)
      call trproj(n,n3,xyz,hs,.false.)
      call hpack2(n3,hs,h2)
   else if (software .eq. "G") then
      call rdghess(n3,h2,hname,fname_pre,finished)
      call hpack1(n3,h2,hs)
      call trproj(n,n3,xyz,hs,.false.)
      call hpack2(n3,hs,h2)
   end if
end if
!
!     include masses
!
do i=1,n3
   do j=1,i 
      u(j,i)=(h2(i,j)+h2(j,i))*0.5*isqm(i)*isqm(j)
   end do
end do
!
!     diagonalitation of the reference hessian
!      
call dsyev ('V','U',n3,u,n3,freq2,aux,lwork,info)
!
!     save modes for sorting
!
modesqm=u
!
!     call frequency sort/alignment before calculating the mad
!       call sortfreq(n,modesqm,modesff,freq2)
!     --> only if different orders exist..
!

vmad=0
zpve2=0
k=0
kk=0
do i=1,n3
   freq2(i)=au2cm*sign(sqrt(abs(freq2(i))),freq2(i))/sqrt(amu2au)
   if (freq2(i).gt.1.d-6) zpve2=zpve2+freq2(i)
   if (freq2(i).lt.-10 )  kk=kk+1                
   if (freq2(i).gt. 10. ) then
      vmad=vmad+abs(freq2(i)-freq(i))
      k=k+1
   end if
end do

if(echo)then
   call preig2(6,freq2,n3)
   write(10,*) "---------------------------------------------------", &
      &  "--------------------"
   write(10,*) 
   write(10,*) "###############################################################"
   write(10,'('' ZPVE comparison (FF/true, kcal) :'',2F12.3)') &
       &  627.51*0.5*zpve/au2cm,627.51*0.5*zpve2/au2cm 
   write(10,'('' MAD of frequencies (cm-1)       :'',F12.3)') &
       &  vmad/k                 
   write(*,*) "Results of the optimization:"
   write(*,'('' Zero point vibr. energy comparison (FF/true, kcal)  :'',2F12.3)') &
       &  627.51*0.5*zpve/au2cm,627.51*0.5*zpve2/au2cm
   write(*,'('' Mean absolute deviation (MAD) of frequencies (cm-1) :'',F12.3)') &
       &  vmad/k
   write(*,*) 
   
   write(10,*) "###############################################################"
end if
if (kk.gt.0) then
    write(10,*)
    write(10,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(10,*)' incorrect curvatures of true Hessian detected'
    write(10,*)'if the imags are >100cm-1 this is a reason to cry!'
    write(10,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
end if


outval=vmad/k

return
end subroutine hess
