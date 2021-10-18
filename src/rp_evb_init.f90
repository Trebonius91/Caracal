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
!     The subroutine rp_evb_init sets up a new reaction path (RP)-EVB
!     coupling term. Based on the energies and structures along the 
!     reaction path stored into ref.dat and struc.xyz, a 1D parameterized 
!     and interpolated path is constructed via cubic spline interpolation.
!     A shape function Zeta will be optimized that archives a smooth 
!     asymptotic behavior of the coupling strength for s=0 and s=1 (where 
!     s is the parameter that determines the 1D reaction path).
!
!     part of EVB
!
subroutine rp_evb_init(filegeo,fileenergy,points,rank)
use evb_mod
use general
implicit none
logical::has_next   ! for geometry read in
integer::points   ! number of reference points
integer::nat,nat3   ! number of atoms (short)
integer::readstat ! file IO status
integer::i,j,k  ! loop indices
character(len=60)::filegeo,fileenergy ! names of files with energies and structures
real(kind=8)::e1,e2,e_diff   ! energies of the QMDFFs and their difference
real(kind=8)::str_e1,str_e2  ! the same ??
real(kind=8),dimension(3*natoms-6)::int_buf  ! buffer for gradient conversion
real(kind=8),dimension(3,natoms)::xyz2,coord  ! xyz structures
real(kind=8),dimension(3*natoms)::g1,g2,gref  ! gradient vectors
real(kind=8),dimension(3,natoms)::gqmdff1,gqmdff2  ! QMDFF gradients
real(kind=8),dimension(3*natoms,3*natoms)::hess1,hess2  ! QMDFF hessians
real(kind=8),dimension(points)::energies_ref,e_path,ff_e1,ff_e2  ! energies
real(kind=8)::diff1,diff2,e_one,e_two  ! correction of QMDFF energies
real(kind=8),dimension(:,:,:),allocatable::xyz_path  ! xyz structures of the reactionpath
real(kind=8),dimension(:),allocatable::int_path  ! internal coordinates for a single struc
real(kind=8),dimension(:,:),allocatable::path_int ! internal coordinates for all structures
real(kind=8),dimension(:),allocatable::v12  !  the coupling strength on the whole path
integer::intdim,infdim  ! spline interpolation: number of internal coordinates and informations
real(kind=8),dimension(:,:),allocatable::ptsin,y2_nd ! spline interpolation: data and 2nd derivates
real(kind=8),dimension(:,:),allocatable::pts   ! spline interpolation: transposed data input array
real(kind=8),dimension(:),allocatable::path_s ! spline interpolation: resulting coordinates for s_act
real(kind=8),dimension(:),allocatable::int_dum  ! dummy for output of interpolation routine
real(kind=8),dimension(points)::s  ! spline interpolation: values of parameter
real(kind=8)::en_ts  ! highest energy along the IRC
real(kind=8)::tot_left,tot_right  ! nonsymmetric borders of energy interpolation region
real(kind=8)::s_act    ! the actual value of the path parameter
real(kind=8)::s_tot_tmp  ! the total IRC length for rescaling
real(kind=8)::corr_stren  ! strengh of the QMDFF correction force constant
real(kind=8)::step  ! TEST
real(kind=8)::pi,coup  ! TEST
real(kind=8),dimension(3*natoms)::freqs  ! TEST for hessian frequencies
integer::int_mode ! for internal coordinate definition
integer::klo,khi  ! borders of the actual spline polynom
integer::ind_ts  ! position of the TS (which structure it is)
integer::rank  ! MPI dummy parameter

pi=3.1415926535897932384626433832795029d0
!
!     first, define the used set of internal coordinates!
!
int_mode=3   ! for RP-EVB, choose always mode 3 as default!
if (read_coord) then 
   int_mode=1  ! read in from file coord_def.inp
else if (dist_matrix) then 
   int_mode=4   ! generation of simple but expensive distance matrix coordinates
end if
call init_int(filegeo,points,rank,int_mode)
!     short forms for variables
nat=natoms
nat3=3*nat
nats=natoms
!
!    Arrays for all needed input informations: energies (all_ens),
!    gradients (all_grads) and hessians (all_hess)
!    of the reference and the single QMDFFs (f1,f2)
!
allocate(xyz_path(3,nat,points))
allocate(all_ens(rp_evb_points))
allocate(allf1_ens(rp_evb_points))
allocate(allf2_ens(rp_evb_points))
allocate(all_xyz(nat*3,rp_evb_points))
allocate(all_int(nat6,rp_evb_points))
allocate(path_int(nat6,points))
allocate(int_path(nat6))
allocate(point_int(nat6,rp_evb_points))
allocate(all_grad(nat6,rp_evb_points))
allocate(allf1_grad(nat6,rp_evb_points))
allocate(allf2_grad(nat6,rp_evb_points))
allocate(all_hess(nat6,nat6,rp_evb_points))
allocate(allf1_hess(nat6,nat6,rp_evb_points))
allocate(allf2_hess(nat6,nat6,rp_evb_points))
!
!    Arrays for the coupling V12, as well as its first and second 
!    derivatives
!
allocate(v12(points))   ! for all given structures on the path
allocate(dv12(rp_evb_points,nat6))   ! only for the treq reference points
allocate(d2v12(rp_evb_points,nat6,nat6))
!
!    Arrays for the cubic spline interpolation of the 1D reaction path
!
allocate(ptsin(points,nat6+2))
allocate(pts(nat6+2,points))
allocate(y2_nd(points,nat6+2))
allocate(path_s(nat6+2))
allocate(int_dum(nat6+2))
allocate(rp_point_s(rp_evb_points))
allocate(rp_point_v12(rp_evb_points))
!
!    TEST: Arrays for frequency corrections of hessians
!
!allocate(indi(natoms),indi3(natoms*3))
!do i=1,natoms
!   indi(i)=i
!   mass(i)=1.0079d0
!end do
!
!    open the ref.input file with informations about the reference points
!    read in the lines and convert geometries, gradients and hessians to 
!    internal coordinates  --> dummy DG-EVB mode = 3
!
call read2int(3)
open (unit=126,file=fileenergy,status="old")
open (unit=127,file=filegeo,status='old')
e_diff=0
!
!  calculate energies for the reaction path and convert all strucutres 
!  to internal coordinates
!
en_ts=-1.d20
do i=1,points
   read(126,*,iostat=readstat) energies_ref(i)
!
!   estimate TS position: if new highest energy was found, take its index
!
   if (energies_ref(i) .ge. en_ts) then
      en_ts=energies_ref(i)
      ind_ts=i
   end if
   if (readstat .ne. 0) then
      write(*,*) "The reference-energies file contains too few lines!"
      call fatal
   end if
   call next_geo(coord,natoms,127,has_next)
   call eqmdff(coord,str_e1,str_e2)
 !  stop "HUphpu"
   xyz_path(:,:,i)=coord
   ff_e1(i)=str_e1
   ff_e2(i)=str_e2
   call xyz_2int(coord,int_path,nat)
   if (.not.has_next) exit
   do j=1,nat6
      path_int(j,i)=int_path(j)
   end do
end do
close(126)

!
!    Correct energies of QMDFF in order to reproduce end points of the IRC
!    exactly   --> can be deactivated with shift_manual
!
if (.not. shift_man) then
   diff1=abs(energies_ref(1)-ff_e1(1))
   diff2=abs(energies_ref(points)-ff_e1(points))
!
!     if the first QMDFF describes the first minimum
!
   if (diff1 .lt. diff2) then
      e_one=ff_e1(1)-E_zero1
      E_zero1=energies_ref(1)-e_one

      e_two=ff_e2(points)-E_zero2
      E_zero2=energies_ref(points)-e_two
!
!    if the first QMDFF describes the second minimum
! 
   else
      e_two=ff_e2(1)-E_zero2
      E_zero2=energies_ref(1)-e_two

      e_one=ff_e1(points)-E_zero1
      E_zero1=energies_ref(points)-e_one
   end if
!
!    refill global QMDFF arrays
!
   do i=1,points
      call eqmdff(xyz_path(:,:,i)*bohr,str_e1,str_e2)
      ff_e1(i)=str_e1
      ff_e2(i)=str_e2
   end do

!
!     print info message that correction was sucessful
!
   if (rank .eq. 0) then
      write(*,*) 
      write(*,*) "QMDFF energies were corrected in order to reproduce energies of"
      write(*,*) "IRC endpoints exactly. The new QMDFF shift parameters are:"
      write(*,*) E_zero1,E_zero2
      write(*,*) 
   end if
end if
!
!     For usage with evb_kt_driver.f90: write single QMDFF energies as well as reference 
!     energies to file for a later plot
!
open(unit=87,file="qmdff_ref.dat",status="unknown")
do i=1,points 
   write(87,*) energies_ref(i),ff_e1(i),ff_e2(i)
end do
close(87)



do i=1,rp_evb_points
!
!   Read in the reference structures for every fixpoint with enhanced 
!   degree of reference information (gradients and hessians)
!

   do j=1,natoms
      do k=1,3
         xyz2(k,j)=all_xyz((j-1)*3+k,i)!*bohr
      end do
   end do
   call xyz_2int(xyz2,int_path,nat)
   do j=1,nat6
      all_int(j,i)=int_path(j)
   end do
   xyz2=xyz2*bohr
!
!   Calculate QMDFF energies of all reference-strucutures
!
   call eqmdff(xyz2,e1,e2)
   allf1_ens(i)=e1
   allf2_ens(i)=e2

!
!   Calculate QMDFF gradients and convert them to internal coordinates
!
   s_transfer=0.5d0
   call egqmdff(xyz2,e1,e2,gqmdff1,gqmdff2)
   do j=1,nat
      do k=1,3
         g1((j-1)*3+k)=gqmdff1(k,j)
         g2((j-1)*3+k)=gqmdff2(k,j)
      end do
   end do
!
!   both internal and cartesian coordinates needed for gradient conversion!
!
!   write(*,*) "int_path",int_path
   call grad2int(xyz2,int_path,allf1_grad(:,i),g1)
   call grad2int(xyz2,int_path,allf2_grad(:,i),g2)
!   write(*,*) "f1_g",i,allf1_grad(:,i)
!   write(*,*) "f2_g",i,allf2_grad(:,i)

!
!   Calculate QMDFF hessians and convert them to internal coordinates
!
   if (rp_evb_mode .eq. 3) then
      call hessqmdff(xyz2,hess1,hess2)
      call hess2int(xyz2,int_path,allf1_hess(:,:,i),hess1,allf1_grad(:,i))
      call hess2int(xyz2,int_path,allf2_hess(:,:,i),hess2,allf2_grad(:,i))
   end if
!   call project_hess(allf1_hess(:,:,i))
!   call project_hess(allf2_hess(:,:,i))
!   write(*,*) "f1_h",i,allf1_hess(:,:,i)
!   write(*,*) "f2_h",i,allf2_hess(:,:,i)
end do
point_int=all_int
!
!    Determine the needed coupling strength at all given points on the 
!    path as well as their derivatives om the RP-EVB points with gradient 
!    and hessian informations
!
!    1.: The coupling strengths
!
do i=1,points
   v12(i)=(ff_e1(i)-energies_ref(i))*(ff_e2(i)-energies_ref(i))
end do
!
!    2.: The coupling gradients
!
do i=1,rp_evb_points
   do j=1,nat6
      dv12(i,j)=((allf1_grad(j,i)-all_grad(j,i))*&
                  &(allf2_ens(i)-all_ens(i)))+((allf1_ens(i)-all_ens(i))*&
                  &(allf2_grad(j,i)-all_grad(j,i)))
!      write(*,*) dv12(i,j),allf1_grad(j,i),all_grad(j,i),allf2_grad(j,i)
   end do
!    write(*,*) all_grad(:,i)
!   write(*,*) "grad_read",all_grad(:,i)
!    write(*,*) "v12_grad_new",dv12(i,:)
end do
!
!    3.: The coupling hessians
!
do i=1,rp_evb_points
   do j=1,nat6
      do k=j,nat6
         d2v12(i,j,k)=(allf2_grad(k,i)-all_grad(k,i))*(allf1_grad(j,i)-all_grad(j,i))+&
                     &(allf1_grad(k,i)-all_grad(k,i))*(allf2_grad(j,i)-all_grad(j,i))+&
                     &(allf2_ens(i)-all_ens(i))*(allf1_hess(j,k,i)-all_hess(j,k,i))+&
                     &(allf1_ens(i)-all_ens(i))*(allf2_hess(j,k,i)-all_hess(j,k,i))
!         write(*,*) i,j,k,allf1_hess(j,k,i),allf2_hess(j,k,i)
         d2v12(i,k,j)=d2v12(i,j,k)
      end do
   end do
end do
!
!    Call the subroutine for the cubic spline interpolation of the reference 
!    reaction path and determine the second derivatives of the polynomina
!    first, fill the coordinate array ptsin
!    Now, fill also the reference energies of the path for the TS region
!
intdim=nat6
ptsin=1.d0
do i=1,points
   ptsin(i,1:nat6)=path_int(:,i)
   ptsin(i,nat6+1)=v12(i)
   ptsin(i,nat6+2)=energies_ref(i)
!   ptsin(i,nat6+2)=
end do
!
!    Convert the input array to the needed transposed values
!
do i=1,points
   do j=1,nat6+2
      pts(j,i)=ptsin(i,j)
   end do
end do

!
!    Call spline builder for energies/coupling strenghts on full path
!
call build_spline(intdim,2,points,pts,ptsin,y2_nd,s,s_tot_tmp)
!
!    Store the optimitzed parameters into global arrays
!
allocate(rp_spl_ref(nat6+2,points))
allocate(rp_spl_d2(points,nat6+2))
allocate(rp_spl_s(points))
rp_spl_ref=pts ! reference information for spline interpolation
rp_spl_d2=y2_nd ! second derivatives for spline interpolation
rp_spl_s=s   !  the parameter reference values for spline interp.
rp_irc_points=points  ! the number of structures along the IRC
s_tot=s_tot_tmp  ! the total length of the IRC path for rescaling
!do i=1,points
!   write(78,*) 2*(i-1),rp_spl_s(i)
!end do
!
!    For acceleration of gradient calculations during dynamics etc: store
!    the previous internal structure/path position index globally
!
allocate(inter_old(1,nat6))
allocate(i_best_old(1))
inter_old=0.d0
i_best_old=0

!
!    The tolerance factor for internal coordinate changes above it 
!    the whole reference IRC will be scanned for the nearest point!
!    Else, only the 5+5=10 neighbored structures will be scanned 
!    It will be scaled with the length of the total IRC devided thorugh
!    the number of RP-energy points and the number of internal coordinates
!
change_tol=s_tot/rp_irc_points*nat6/100d0

!
!    determine position of the TS along the IRC (s-parameter)
!
s_ts=rp_spl_s(ind_ts)
!
!    Next, determine the positions of the gradient/hessian reference points 
!    on the IRC in terms of the parameter s 
!    Store also the V12 values of these points
!    In order to avoid bad results, the whole IRC will be scanned by using 
!    the spline_all option
!
!write(*,*) path_int(:,3)
spline_all=.true.
do i=1,rp_evb_points 
   int_path=all_int(:,i)!/bohr
   call spline_dist(nat6,s_act,int_path,int_dum,1)
  
   call interp_spline(nat6+1,s_act,path_s,klo,khi)
   rp_point_v12(i)=path_s(nat6+1)
   rp_point_s(i)=s_act
  ! write(*,*) "i sact",i,s_act,int_path(1:3)
end do
spline_all=.false.
!
!    Check if the first and last structure for the gradient/hessian reference 
!    is also the total first/last structure on the whole path
!
if (rp_point_s(1) .gt. 1D-5) then
   write(*,*) "OOPS, the first gradient/hessian reference structure for RP-EVB"
   write(*,*) "has to be also the first structure on the whole reference IRC!"
   write(*,*) "Please change the ref.input file!"
   call fatal
else if (rp_point_s(rp_evb_points) .lt. 0.9999d0) then
   write(*,*) "OOPS, the last gradient/hessian reference structure for RP-EVB"
   write(*,*) "has to be also the last structure on the whole reference IRC!"
   write(*,*) "Please change the ref.input file!"
   call fatal
end if
!
!    To avoid numerical problems: set s of the first point to 0 and s of the last point to 1
!
rp_point_s(1)=0.d0
rp_point_s(rp_evb_points)=1.d0  
!
!     For TS that are not located in the centre, also determine the 
!     borders nonsymmetric
!
tot_left=s_ts*rp_mid_tot
tot_right=(1d0-s_ts)*rp_mid_tot
!
!     then, determine borders of outer RP_EVB regions; look, where the next
!     reference structure along the path is located!
!     Take always the largest possible interval, i.e. for left_lo take the klo
!     index and for right_hi the khi index
!     Only if no manual read in for these borders is activated!
!
!     if the NO_EVB option is activated, use the same values for the determination
!     of pure QMDFF/ transition/ pure direct interpolation regions
!
if (.not. no_evb) then
   if (.not. s_bord_man) then
      trans_l_lo=s_ts-tot_left
      trans_l_hi=s_ts-tot_left+rp_mid_trans
      trans_r_lo=s_ts+tot_right-rp_mid_trans
      trans_r_hi=s_ts+tot_right
   end if
else 
   iq_l_lo=s_ts-tot_left
   iq_l_hi=s_ts-tot_left+rp_mid_trans
   iq_r_lo=s_ts+tot_right-rp_mid_trans
   iq_r_hi=s_ts+tot_right
!   write(*,*) "iq",iq_l_lo,iq_l_hi,iq_r_lo,iq_r_hi
end if
!write(*,*) "trans_l_lo",trans_l_lo
!write(*,*) "trans_l_hi",trans_l_hi
!write(*,*) "trans_r_lo",trans_r_lo
!write(*,*) "trans_r_hi",trans_r_hi

!
!     For better interpolation behavior: shift the borders so that they coincide
!     with borders of single spline interpolation parts
!
if (.not. no_evb) then
   call interp_spline(nat6+1,trans_l_lo,path_s,klo,khi)
   trans_l_lo=rp_spl_s(klo)

   call interp_spline(nat6+1,trans_l_hi,path_s,klo,khi)
   trans_l_hi=rp_spl_s(khi)

   call interp_spline(nat6+1,trans_r_lo,path_s,klo,khi)
   trans_r_lo=rp_spl_s(khi)

   call interp_spline(nat6+1,trans_r_hi,path_s,klo,khi)
   trans_r_hi=rp_spl_s(klo)
end if
!write(*,*) "trans_l_lo",trans_l_lo
!write(*,*) "trans_l_hi",trans_l_hi
!write(*,*) "trans_r_lo",trans_r_lo
!write(*,*) "trans_r_hi",trans_r_hi
!write(*,*) "iq_l_lo",iq_l_lo
!write(*,*) "iq_l_hi",iq_l_hi
!write(*,*) "iq_r_lo",iq_r_lo
!write(*,*) "iq_r_hi",iq_r_hi


inter_old=0.d0
i_best_old=0


return
end subroutine rp_evb_init
