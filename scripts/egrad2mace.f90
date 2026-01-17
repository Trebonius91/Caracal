
!
!    This program converts the results of a energy+gradient 
!     calculation with explore.x (from Caracal) into a MACE
!     training set file. 
!    For this, the files trajectory.xyz, energies.dat and 
!     gradients.dat need to be present
!

program egrad2mace
implicit none
integer::i,j 
integer::natoms,nframes
integer::readstat
real(kind=8),allocatable::energies(:),xyz(:,:,:),grad(:,:,:)
character(len=2),allocatable::names(:,:)
real(kind=8)::bohr,evolt

write(*,*)
write(*,*) "This program converts the results of a energy+gradient"
write(*,*) " calculation with explore.x (from Caracal) into a MACE"
write(*,*) " training set file."
write(*,*) "For this, the files trajectory.xyz, energies.dat and"
write(*,*) " gradients.dat need to be present"

bohr=0.52917721092d0
evolt=27.21138503d0

nframes=0
open(unit=27,status="old",file="energies.dat",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'energies.dat' is not there!"
   stop
end if
do
   read(27,*,iostat=readstat) 
   if (readstat .ne. 0) exit
   nframes=nframes+1
end do
close(27)

open(unit=28,status="old",file="trajectory.xyz",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'trajectory.xyz' is not there!"
   stop
end if

read(28,*) natoms
close(28)

allocate(energies(nframes))
allocate(xyz(3,natoms,nframes))
allocate(grad(3,natoms,nframes))
allocate(names(natoms,nframes))

open(unit=27,status="old",file="energies.dat")
do i=1,nframes
   read(27,*) energies(i)
end do
close(27)

open(unit=28,status="old",file="trajectory.xyz")
do i=1,nframes
   read(28,*)
   read(28,*) 
   do j=1,natoms
      read(28,*) names(j,i),xyz(:,j,i)
   end do
end do
close(28)

open(unit=29,status="old",file="gradients.dat",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'gradients.dat' is not there!"
   stop
end if

do i=1,nframes
   read(29,*)
   do j=1,natoms
      read(29,*) grad(:,j,i)
   end do
end do
close(29)


open(unit=51,status="replace",file="training_set.xyz")
do i=1,nframes
   write(51,*) natoms
   write(51,'(a,f20.10,a)') 'Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3 Nmols=1 &
                           &REF_energy=',energies(i)*evolt,' pbc="F F F"'
   do j=1,natoms
      write(51,'(a,a,3f14.8,a,3f14.8)') names(j,i),"  ",xyz(:,j,i),"  0  ",grad(:,j,i)*bohr*evolt
   end do
end do
close(51)
write(*,*)
write(*,*) "Execution of egrad2mace sucessful!"
write(*,*) "File 'training_set.xyz' has been written."
write(*,*)


end program egrad2mace
