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
!     function asym: read in number of element and assign symbol
!
!     part of QMDFF
! 

function asym(I)
CHARACTER(len=2)::ASYM
CHARACTER(len=2)::ELEMNT(107), AS
DATA ELEMNT/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/ 
AS=ELEMNT(I)
!     convert the symbols to upper case
CALL UPPER(AS)
ASYM=AS
if (i.eq.103) asym='XX'
RETURN
end function asym
