

program call_external
implicit none 
real(kind=8),allocatable::xyz(:,:),grad(:,:)
integer::i,j,natoms
real(kind=8)::energy
character(len=3)::adum
integer::info

open(unit=37,file="coords.xyz",status="old")
read(37,*) natoms

allocate(xyz(3,natoms))
allocate(grad(3,natoms))

read(37,*)
do i=1,natoms
   read(37,*) adum, xyz(:,i) 
end do
close(37)
xyz=xyz/0.52917721092d0


call potential(xyz,natoms,energy,grad,info)

open(unit=38,file="egrad_out.dat",status="replace")
write(38,*) energy
do i=1,natoms
   write(38,*) grad(:,i)
end do



end program call_external
