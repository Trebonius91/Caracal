
!
!    This program converts the results of a energy+gradient 
!     calculation with explore.x (from Caracal) into a MACE
!     training set file. 
!    For this, the files trajectory.xyz, energies.dat and 
!     gradients.dat need to be present
!

program egrad2mace
implicit none
integer::i,j,k 
integer::natoms,nframes
integer::readstat
integer::at_ind
real(kind=8),allocatable::energies(:),xyz(:,:,:),grad(:,:,:)
character(len=2),allocatable::names(:,:)
real(kind=8)::bohr,evolt
character(len=200)::line,adum
character(len=4)::folder
integer::status

write(*,*)
write(*,*) "This program converts the results of energy+gradient"
write(*,*) " calculations with Orca (e.g., wB97X-V) into a MACE"
write(*,*) " training set file."
write(*,*) "For this, this program needs to be started in the "
write(*,*) " main folder, containing subfolders 0000, 0001, ..."
write(*,*) " with Orca calculations, basename 'ener'"
write(*,*) "An output file train_mace.xyz will be written."

bohr=0.52917721092d0
evolt=27.21138503d0

nframes=1

!
!    Number of frames: number of Orca folders
!
do
   write(folder, '(I4.4)') nframes
   call execute_command_line('test -d ' // folder, exitstat=status)
   if (status /= 0) exit
   nframes = nframes + 1
end do
if (nframes .lt. 2) then
   write(*,*) "No Orca result folders could be found!"
end if
nframes=nframes-1
write(*,*)
write(*,*) "Number of calculation folders/frames: ",nframes

call chdir("0001")
open(unit=28,file="ener.engrad",status="old")
read(28,*)
read(28,*)
read(28,*)
read(28,*) natoms
close(28)
call chdir("..")
write(*,*)
write(*,*) "Number of atoms: ",natoms


allocate(energies(nframes))
allocate(xyz(3,natoms,nframes))
allocate(grad(3,natoms,nframes))
allocate(names(natoms,nframes))

!
!    Loop through all folders and read in geometries, 
!    gradients and energies
!

do i=1,nframes
   write(folder, '(I4.4)') nframes
   open(unit=29,file=folder//"/ener.engrad",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file ",folder,"/ener.engrad is missing!"
      stop
   end if
   do j=1,7
      read(29,*) 
   end do
   read(29,*) energies(i)
   read(29,*)
   read(29,*)
   read(29,*)
   do j=1,natoms
      do k=1,3
         read(29,*) grad(k,j,i)
      end do
   end do
   read(29,*)
   read(29,*)
   read(29,*)
   do j=1,natoms
      read(29,*) at_ind,xyz(:,j,i)
      call atomname(at_ind,names(j,i))
   end do
   close(29)
end do



open(unit=51,status="replace",file="training_set.xyz")
do i=1,nframes
   write(51,*) natoms
   write(51,'(a,f20.10,a)') 'Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3 Nmols=1 &
                           &REF_energy=',energies(i)*evolt,' pbc="F F F"'
   do j=1,natoms
      write(51,'(a,a,3f14.8,a,3f14.8)') names(j,i),"  ",xyz(:,j,i)/bohr,"  0  ",grad(:,j,i)*bohr*evolt
   end do
end do
close(51)
write(*,*)
write(*,*) "Execution of orca2mace sucessful!"
write(*,*) "File 'training_set.xyz' has been written."
write(*,*)


end program egrad2mace

!
!     subroutine atomname: is an atomic charge is given, extract the 
!        corresponding element label 
!
!     part of EVB 
!

subroutine atomname(nat,key1)
implicit none
CHARACTER(len=3)::KEY1  ! the name of the element
CHARACTER(len=2)::ELEMNT(107)  ! array with all names 
integer::nat  ! atomic charge 

DATA ELEMNT/'H ','He', &
 & 'Li','Be','B ','C ','N ','O ','F ','Ne', &
 & 'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
 & 'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu', &
 & 'Zn','Ga','Ge','As','Se','Br','Kr', &
 & 'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag', &
 & 'Cd','In','Sn','Sb','Te','I ','Xe', &
 & 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
 & 'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt', &
 & 'Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
 & 'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/

!
!     Return the i'th entry of the elements array of i is the given
!     atomic charge of the atom 
!
key1=elemnt(nat)


return
end subroutine atomname

