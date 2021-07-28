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
!     subroutine rd: Read in coordinates and atom types from Turbomole output
!
!     part of QMDFF
!
subroutine rd(echo,fname,n,xyz,iat)
implicit real(kind=8) (a-h,o-z)
dimension::xyz(3,n),iat(n),xx(10)
character(len=128)::line
character(len=80)::s(3)
character(len=*)::fname  ! variable length for dummy parameter (..)
logical::echo
CHARACTER(len=2)::el1(98),el2(98)
DATA el1/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf'/
DATA el2/'H ','HE', &
 & 'LI','BE','B ','C ','N ','O ','F ','NE', &
 & 'NA','MG','AL','SI','P ','S ','CL','AR', &
 & 'K ','CA','SC','TI','V ','CR','MN','FE','CO','NI','CU', &
 & 'ZN','GA','GE','AS','SE','BR','KR', &
 & 'RB','SR','Y ','ZR','NB','MO','TC','RU','RH','PD','AG', &
 & 'CD','IN','SN','SB','TE','I ','XE', &
 & 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY', &
 & 'HO','ER','TM','YB','LU','HF','TA','W ','RE','OS','IR','PT', & 
 & 'AU','HG','TL','PB','BI','PO','AT','RN', &
 & 'FR','RA','AC','TH','PA','U ','NP','PU','AM','CM','BK','CF'/

ich=142
open(unit=ich,file=fname)
read(ich,'(a)') line
do i=1,n
   read(ich,*) xyz(1:3,i)
enddo

if (echo) then
   write(10,*) '========================='
   write(10,*) 'reading ... ',trim(fname)
   write(10,*) '========================='
endif

rewind ich
read(ich,'(a)') line
do i=1,n
   read(ich,'(a)')line
   call cutline(line,xx,s)
   do j=1,98
      k1=index(s(1),el1(j))
      k2=index(s(1),el2(j))
      if (k1.eq.1.or.k2.eq.1) then
         iat(i)=j
         exit
      end if
   end do
end do
close(ich)
return
end subroutine rd
