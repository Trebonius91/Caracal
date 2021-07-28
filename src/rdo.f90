!
!     subroutine rdo: Read in coordinates and atom types from orca output
!
!     part of QMDFF
!
subroutine rdo(echo,fname,n,xyz,iat)
implicit real(kind=8) (a-h,o-z)
dimension::xyz(3,n),iat(n),xx(10)
character(len=128)::line
character(len=2)::a2
character(len=*)::fname
logical::echo

if (echo) then
   write(10,*) '========================='
   write(10,*) 'reading ... ',trim(fname)
   write(10,*) '========================='
end if

ich=142
open(unit=ich,file=fname)
do
   read(ich,'(a)') line
   if(index(line,'$atoms').ne.0) then
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) a2,dum,xyz(1:3,i)
         call elem(a2,iat(i))
      end do
      exit
   end if
end do

close(ich)
return
end subroutine rdo
