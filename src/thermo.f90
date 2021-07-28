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
!     subroutine thermo: Calculate (statistical) thermodynamics properties
!
!     part of QMDFF
!


subroutine thermo(nat,it,coord,freq,h298)
implicit none 
                       
integer::nat
integer::it(nat)
real(kind=8)::freq(3*nat),coord(3,nat)
character(len=128)::A,arg
character(len=30)::sym                                  
real(kind=8)::xx(10),sthr,temp,scale_factor,vibs(3*nat)
real(kind=8)::aa,bb,cc,vibthr,zp
real(kind=8)::g298,h298,escf,symnum,ts298,wt,avmom
integer::nn,nvib,nimag,i,n,nvib_theo
logical::lin,atom,da

write(10,*) 
write(10,*) '=================='  
write(10,*) 'thermo calculation'  
write(10,*) '=================='    
write(10,*)
!
!     frequencies read in are considered
!     as being real if .gt. this value in cm-1      
!     this threshold requires projected freqs.!
!
vibthr=1.0

sym='c1'
atom=.false.
lin=.false.
temp=298.15
sthr=100.
scale_factor=1.0
nvib=0
nimag=0
write(10,*) "Molecule properties:"

call axis(nat,it,coord,aa,bb,cc,avmom,wt) 

nvib_theo=3*nat-6
if (cc.lt.1.d-6) lin=.true.
if (lin) nvib_theo=3*nat-5
write(10,*)'linear? (True/False)   ',lin

if (aa+bb+cc.lt.1.d-6)then
   atom=.true.
   nvib=0
   nvib_theo=0            
   goto 300
end if

open(unit=1,file='control')
19 continue 
read(1,'(a)',end=190)a
if (index(a,'$symmetry').ne.0) sym=a
goto 19
190  close(1,status='delete')

300  symnum=1
!
!     distinct between different molecular symmetries
!
if (index(sym,'ci').ne.0) symnum=1
if (index(sym,'cs').ne.0) symnum=1
if (index(sym,'c2').ne.0) symnum=2
if (index(sym,'c3').ne.0) symnum=3
if (index(sym,'s6').ne.0) symnum=3
if (index(sym,'c4').ne.0) symnum=4
if (index(sym,'c5').ne.0) symnum=5
if (index(sym,'c6').ne.0) symnum=6
if (index(sym,'c7').ne.0) symnum=7
if (index(sym,'c8').ne.0) symnum=8
if (index(sym,'c9').ne.0) symnum=9
if (index(sym,'c10').ne.0) symnum=10
if (index(sym,'c11').ne.0) symnum=11
if (index(sym,'d2').ne.0) symnum=4
if (index(sym,'d3').ne.0) symnum=6
if (index(sym,'d4').ne.0) symnum=8
if (index(sym,'d5').ne.0) symnum=10
if (index(sym,'d6').ne.0) symnum=12
if (index(sym,'d7').ne.0) symnum=14
if (index(sym,'d8').ne.0) symnum=16
if (index(sym,'d9').ne.0) symnum=18
if (index(sym,'td').ne.0) symnum=12
if (index(sym,'oh').ne.0) symnum=24
if (index(sym,'ih').ne.0) symnum=60
if (index(sym,'c').ne.0.and.lin) symnum=1
if (index(sym,'d').ne.0.and.lin) symnum=2
 
do i=1,3*nat
   if (freq(i).gt.vibthr) then
      nvib=nvib+1
      vibs(nvib)=freq(i)
   end if
end do
write(10,'(A,I5)') ' Number of vibrations: ',nvib
write(10,*) 'Molecular symmetry:   ',trim(sym)                
write(10,'(A,F6.2)') ' Rotational symmetry number ',symnum
!
!     scale   
!      
vibs(1:nvib)=vibs(1:nvib)*scale_factor
!
!     do calc. with old mopac code
!    
escf=0
call thermocal(aa,bb,cc,avmom,lin,atom,symnum,wt,vibs,nvib,escf, &
    &  temp,sthr,h298,g298,ts298,zp)

return
end subroutine thermo
