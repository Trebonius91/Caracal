program caracal_quality
implicit none 
integer::i,j,k,l
integer::readstat
integer::nframes,natoms
character(len=50)::cdum
character(len=150)::arg
real(kind=8)::adum
real(kind=8)::bohr,evolt
real(kind=8),allocatable::ener_ref(:),ener_mace(:)
real(kind=8),allocatable::grad_ref(:,:,:),grad_mace(:,:,:)
real(kind=8),allocatable::xi_vals(:)
real(kind=8)::ener_mae,grad_mae,gradnorm_mae
real(kind=8)::ener_rmse,grad_rmse
real(kind=8)::nhisto_2d_abs_range,nhisto_2d_angle_range
real(kind=8)::gradnorm_ref,gradnorm_mace,cos_theta
real(kind=8)::angle
real(kind=8)::xi_min,xi_max,xi_range,ener_min
real(kind=8),allocatable::ener_xi_min(:)
integer,allocatable::xi_bin_loc(:)
integer::nhisto_2d_abs,nhisto_2d_angle
integer::nhisto_xi
integer,allocatable::histo_grad_2d(:,:)
integer::ngrads
integer::max_histo
parameter (bohr=0.52917721092d0)
parameter (evolt=27.21138503d0)

write(*,*) 
write(*,*) "Program caracal_quality: "
write(*,*) "This program evaluates the energies and gradients"
write(*,*) "of MD frames collected from a calc_rate structure"
write(*,*) "generation calculation by comparing them to a "
write(*,*) "reference method of interest."
write(*,*) "For this, the MD energies/gradients must be written"
write(*,*) "by calc_rate with the print_gen / print_train mace"
write(*,*) "commands, both traj_gen.xyz and traj_gen_xi.dat"
write(*,*) "must be placed in the current folder."
write(*,*) "Further, the files gradients.dat and energies.dat"
write(*,*) " in the format of explore.x (job egrad) must be there,"
write(*,*) " containing the reference method energies/gradients"
write(*,*) "Further commands for plots:"
write(*,*) " -histo_2d_abs_range=[value]: The x axis range for"
write(*,*) "    2D gradient histogram (absolute value)"
write(*,*) " -histo_2d_angle_range=[value]: The y axis range for "
write(*,*) "    2D gradient histogram (gradient angle deviation)"
write(*,*) " -nhisto_2d_abs=[number]: Number of x bins "
write(*,*) " -nhisto_2d_angle=[number]: Number of y bins "
write(*,*) 
!
!     Ranges for 2D plots of gradient angles and absolute values
!
nhisto_2d_abs_range=5.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:20))  .eq. "-histo_2d_abs_range=") then
      read(arg(21:),*,iostat=readstat) nhisto_2d_abs_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_2d_abs_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_angle_range=10d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:22))  .eq. "-histo_2d_angle_range=") then
      read(arg(23:),*,iostat=readstat) nhisto_2d_angle_range
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -histo_2d_angle_range=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_abs=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:15))  .eq. "-nhisto_2d_abs=") then
      read(arg(16:),*,iostat=readstat) nhisto_2d_abs
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_2d_abs=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

nhisto_2d_angle=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:17))  .eq. "-nhisto_2d_angle=") then
      read(arg(18:),*,iostat=readstat) nhisto_2d_abs
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -nhisto_2d_angle=..., something went wrong!"
         write(*,*)
      end if
   end if
end do

!
!    First, determine the number of frames and atoms
!
open(unit=17,file="energies.dat",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'energies.dat' does not exist!"
   stop
end if
nframes=0
do 
   read(17,*,iostat=readstat) adum
   if (readstat .ne. 0) exit
   nframes=nframes+1
end do
close(17)

open(unit=18,file="traj_gen.xyz",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'traj_gen.xyz' does not exist!"
   stop
end if
read(18,*) natoms
close(18)


allocate(ener_mace(nframes))
allocate(ener_ref(nframes))
allocate(xi_vals(nframes))
allocate(grad_mace(3,natoms,nframes))
allocate(grad_ref(3,natoms,nframes))
allocate(histo_grad_2d(nhisto_2d_abs,nhisto_2d_angle))
!
!    Then, read in the MACE reference energies and gradients
!
open(unit=18,file="traj_gen.xyz",status="old")
do i=1,nframes
   read(18,*) 
   read(18,*) cdum,cdum,ener_mace(i)
   do j=1,natoms
      read(18,*) cdum,adum,adum,adum,grad_mace(:,j,i)
   end do
end do
close(18)

!
!    Then read in the analytical energies and gradients
!
open(unit=19,file="energies.dat",status="old")
do i=1,nframes
   read(19,*) ener_ref(i)
end do
close(19)

open(unit=20,file="gradients.dat",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'gradients.dat' does not exist!"
   stop
end if
do i=1,nframes
   read(20,*)
   do j=1,natoms
      read(20,*) grad_ref(:,j,i)
   end do
end do
!
!     Then read in the reaction progress/xi values
!
open(unit=21,file="traj_gen_xi.dat",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'traj_gen_xi.dat' does not exist!"
   stop
end if
do i=1,nframes
   read(21,*) xi_vals(i)
end do
close(21)

!
!     Convert reference to same unit 
!
ener_ref=ener_ref*evolt
grad_ref=-grad_ref*evolt/bohr

!
!     Calculate the mean absolute errors (MAEs) of energies and gradients
!

ener_mae=0.d0
ener_rmse=0.d0
do i=1,nframes
   ener_mae=ener_mae+abs(ener_ref(i)-ener_mace(i))
   ener_rmse=ener_rmse+(ener_ref(i)-ener_mace(i))**2
end do
ener_mae=ener_mae/nframes
ener_rmse=sqrt(ener_rmse)/nframes

write(*,*) "The MAE of energies is (in meV):",ener_mae*1000.d0
!write(*,*) "The RMSE of energies is (in meV):",ener_rmse*1000.d0


grad_mae=0.d0
grad_rmse=0.d0
do i=1,nframes
   do j=1,natoms
      do k=1,3
         grad_mae=grad_mae+abs(grad_ref(k,j,i)-grad_mace(k,j,i))
         grad_rmse=grad_rmse+(grad_ref(k,j,i)-grad_mace(k,j,i))**2
      end do
   end do
end do
grad_mae=grad_mae/(nframes*natoms*3)
grad_rmse=sqrt(grad_rmse)/(nframes*natoms*3)

write(*,*) "The MAE of gradients is (in meV/A):",grad_mae*1000.d0
!write(*,*) "The RMSE of gradients is (in meV/A):",grad_rmse*1000.d0
!
!    Generate 2D histogram of gradient norm vs angle to show their accuracies
!
ngrads=0
histo_grad_2d=0
do i=1,nframes
   do j=1,natoms
      gradnorm_mace=sqrt(grad_mace(1,j,i)**2+grad_mace(2,j,i)**2+grad_mace(3,j,i)**2)
      gradnorm_ref=sqrt(grad_ref(1,j,i)**2+grad_ref(2,j,i)**2+grad_ref(3,j,i)**2)
      cos_theta=dot_product(grad_ref(:,j,i),grad_mace(:,j,i))/gradnorm_mace/gradnorm_ref
      if (cos_theta .lt. -1.d0) cos_theta=-1.d0
      if (cos_theta .gt. 1.d0) cos_theta=1.d0
      angle=acos(cos_theta)
      angle=abs(angle*180.d0/acos(-1.d0))

      k=nint(angle/nhisto_2d_angle_range*nhisto_2d_angle)+1

      l=nint(gradnorm_ref/nhisto_2d_abs_range*nhisto_2d_abs)+1

      if (k .gt. nhisto_2d_angle .or. l .gt. nhisto_2d_abs) cycle
      histo_grad_2d(l,k)=histo_grad_2d(l,k)+1
      ngrads=ngrads+1
   end do
end do
max_histo=maxval(histo_grad_2d)

open(unit=60,file="gradnorm_vs_angle_histo.dat",status="replace")
do i=1,nhisto_2d_abs
   do j=1,nhisto_2d_angle
      write(60,*) real(i-1)/real(nhisto_2d_abs)*nhisto_2d_abs_range,real(j-1)/real(nhisto_2d_angle)* &
              & nhisto_2d_angle_range,real(histo_grad_2d(i,j))/max_histo
   end do
   write(60,*)
end do
close(60)
write(*,*) "2D gradient direction histogram written to 'gradnorm_vs_angle_histo.dat'"
write(*,*)

!
!     The 2D force direction error histogram
!
write(*,*) "Execute gnuplot file 'plot_2d_grad_histo.gnu'"
open(unit=61,file="plot_2d_grad_histo.gnu",status="replace")
write(61,*) "set encoding iso_8859_1"
write(61,*) "set terminal png lw 6.0 size 2400,2000 font 'Helvetica,56'"
write(61,*) "set output '2d_gradient_histo.png'"
write(61,*) "set xlabel 'absolute atomic force (eV/{\305})'"
write(61,*) "set ylabel 'error in direction (°)'"
write(61,*) "set zlabel 'relative frequency'"
write(61,*) "set grid"
write(61,*) "set lmargin at screen 0.18"
write(61,*) "set rmargin at screen 0.85"
write(61,*) "set bmargin at screen 0.15"
write(61,*) "set tmargin at screen 0.95"
write(61,*) "unset key"
write(61,*) "set pm3d map interpolate 5,5"
write(61,*) "set palette defined ( 0 'white',\"
write(61,*) "     0.01 'black',\"
write(61,*) "     0.15 'blue',\"
write(61,*) "     0.4 'green',\"
write(61,*) "     0.7 'yellow',\"
write(61,*) "     1 'red' )"
write(61,*) "splot 'gradnorm_vs_angle_histo.dat' u 1:2:3 with pm3d"
close(61)
call system("gnuplot plot_2d_grad_histo.gnu")
write(*,*) "Plot picture written to '2d_gradient_histo.png'!"

write(*,*)
write(*,*) "Program mlp_quality exited normally."
write(*,*)


!
!    Write file with energies/energy errors vs reaction path Xi value
!
ener_min=minval(ener_ref)
open(unit=28,file="xi_ener.dat",status="replace")
write(28,*) "# xi value    ref. energy (eV)    energy diff. (meV)"
do i=1,nframes
   write(28,*) xi_vals(i),ener_ref(i)-ener_min,abs(ener_ref(i)-ener_mace(i))*1000.d0
end do
close(28)
!
!    The energy/energy error vs Xi scattering plot
!
write(*,*) "Execute gnuplot file 'plot_xi_ener.gnu'"
open(unit=62,file="plot_xi_ener.gnu",status="replace")
write(62,*) "set terminal png truecolor lw 6.0 size 2400,2000 font 'Helvetica,56'"
write(62,*) "set output 'plot_xi_ener.png'"
write(62,*) "set encoding utf8"
write(62,*) "set xlabel 'ξ (a.u.)'"
write(62,*) "set ylabel 'total energy (eV)'"
write(62,*) "set cblabel 'energy error (meV)'"
write(62,*) "set style fill transparent solid 0.25 noborder"
write(62,*) "set style circle radius 0.013"
write(62,*) "unset key"
write(62,*) "set grid"
write(62,*) "set palette rgb 33,13,10"
write(62,*) "plot 'xi_ener.dat' u 1:2:3 with circles lc palette z"
close(62)
call system("gnuplot plot_xi_ener.gnu")
write(*,*) "Plot picture written to 'plot_xi_ener.png'!"


write(*,*)
write(*,*) "Program caracal_quality exited normally."
write(*,*)

end program caracal_quality
