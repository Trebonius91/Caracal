!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CARACAL - Ring polymer molecular dynamics and rate constant calculations
!             on black-box generated potential energy surfaces
!
!   Copyright (c) 2023 by Julien Steffen (mail@j-steffen.org)
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
!     subroutine main_gen: the main entry for any attempt to optimize a new 
!     QMDFF force field. It is called from the qmdffgen program
!
!     part of QMDFF
!
subroutine main_gen(fname_pre,length)
use qmdff  ! include global variables
implicit none
!     define all variables that are locally used!

integer::n  ! number of atoms in system
real(kind=8),allocatable::xyz (:,:)  ! atomic coordinates
real(kind=8),allocatable::xyz1(:,:)
real(kind=8),allocatable::xyz2(:,:)
real(kind=8),allocatable::g(:,:)  ! gradient
real(kind=8),allocatable::grep(:,:) 
real(kind=8),allocatable::wbo(:,:)  ! Wiberg-Mayer bond orders
real(kind=8),allocatable::q(:)  
real(kind=8),allocatable::chir(:)
real(kind=8),allocatable::cn  (:)
real(kind=8),allocatable::c6xy(:,:)
real(kind=8),allocatable::epair(:)
real(kind=8),allocatable::h(:,:)
real(kind=8),allocatable::freq(:)
integer,allocatable::at(:)
integer,allocatable::nb(:,:)
integer,allocatable::mlist(:)
integer,allocatable::molist(:)
integer,allocatable::hyb(:)
integer,allocatable::covpair(:)
integer,allocatable::cring(:,:)
integer,allocatable::ringsize(:)
!     loop indices:
integer::i,j,k,l,nvar,nst,izero,nn,i1,i2,iz1,iz2
integer::mol1,nbf,nel,m,mm,time,icode,mdmode,mtors
integer::length
character(len=80)::a,hname,rcfile,homedir,method,fname1,fname,fname2
character(len=80)::outname
character(len=80)::fname3
character(len=length+6)::solvname
character(len=length)::fname_pre
character(len=80)::arg(20)
real(kind=8)::e,temp,valijk,f1(3),f2(3),xx(10),one,e0input,chrg,sum1
real(kind=8)::gnorm,val,einf,dens,c6,r42,qscal,dipol2(3),dipol(3)
real(kind=8)::d1,d2,dum,eat,h298,e0,er,el,vz(94),zab(94,94),r0ab(94,94)
real(kind=8)::sr42(94,94),e0ff,t,enci(3),hneglect,sum2,wbocut_in
real(kind=8)::scalehb(94),scalexb(94)


logical::ex,bondcalc,tm,echo,domd,vdw,codetest,hex,hesscalc,dist
logical::rds,doopt,okbas,analyse,refit,restart,mclust,tstat,numhess

character(len=8)::date
character(len=10)::times
character(len=5)::zone
integer,dimension(8)::values
integer::finished ! if formchk should be executed

!
!     get the name of the hessian input file with reference data
!     --> depends on the used reference software
!

if (software .eq. "O") then
   fname=fname_pre//'.hess'
   outname=fname_pre//'.out'
else if (software .eq. "C") then
   fname=fname_pre//'.out'
else if (software .eq. "G") then
   fname=fname_pre//'.log'
else if (software .eq. "V") then
   fname=fname_pre//".OUTCAR"
end if 
!
!     set calculation settings default values
!
bondcalc=.true. ! calculate bond energies separately
codetest=.false.
doopt=.false.
domd=.false.
echo=.false.
rds =.false.
vdw =.false.
mclust=.false.
analyse=.false.  ! more datails written after optimization
hesscalc=.true.
numhess=.false.

nrot1=0
nrot2=0
dipol=0
dipol2=0
dens=0.0
chrg=0.0
!
!     global scale factor for all reference hessian elements! 
!
scalh=1.
!
!     additional torsions aroung metal centers?
!
mtors=0 
!
!     important: hessian elements smaller than the value are neglected!
!
hneglect=1.d-4
!
!     when is an inversion symmetrized to planarity (degree)
!     and when are angles (for torsion assignment only) linear      
!     a peptide and the tm_64 test are better with athr=2 than for athr=3
!     requires sometimes user input
!
athr=2.0d0
!
!     imporant: default CM5 scaling factor
!
qscal=1.15
!
!     when is a neighbor covalently bonded?
!     was bondorder .ge. 0.5 in the published version   
!
wbocut_in=0.5
!
!     H-bond(HB) parameters/scaling for interesting elements
!     HB as in qmsolv (D3 is in block.f)
!
scalehb=0
scalehb(7  )=0.8
scalehb(8  )=0.3
scalehb(9  )=0.1
scalehb(15 )=2.0
scalehb(16 )=2.0
scalehb(17 )=2.0
scalehb(34 )=2.0
scalehb(35 )=2.0
!
!     Halogenbond(XB) parameters/scaling for halogenes
!
scalexb=0
scalexb(17 )=0.30
scalexb(35 )=0.60
scalexb(53 )=0.80
scalexb(85 )=1.00
!
!     write out informations about HB and XB parameters
!
call date_and_time(date,times,zone,values)
write(10,*)
write(10,*) '==================================='
write(10,*) 'QMDFF GENERATION:'
write(10,*) 'Program written by S.Grimme'
write(10,*) 'Fortran90 translation by J. Steffen'
if (check_coord) then
   write(10,*) 'The CHECK_COORD option was activated.'
   write(10,*) 'Only coordinate analysis done, QMDFF'
   write(10,*) 'generation will be skipped!'
end if
write(10,*) '==================================='
write(10,'('' '',I2,''.'',I2,''.'',I4,'', '',I2,'':'',I2,'':'',I2)') & 
    & values(3),values(2),values(1),values(5),values(6),values(7)
write(10,*)
write(10,*) "General parameters:"
write(10,'('' Hessian neglect threshold (au)'',E10.3)')hneglect
write(10,'('' linear or plane cutoff (deg.) '',F10.3)')athr
write(10,*) "Non covalent and hydrogen/halogen bond parameters:"
write(10,'('' D3             '',3F18.4)')a1,a2,s8
write(10,'('' HB(N,O,F)      '',3F18.4)')scalehb(7:9)
write(10,'('' HB(P,S,Cl)     '',3F18.4)')scalehb(15:17)
write(10,'('' HB(As,Se,Br)   '',3F18.4)')scalehb(33:35)
write(10,'('' HX(Cl,Br,I)    '',3F18.4)')  &
         &   scalexb(17),scalexb(35),scalexb(53)
write(10,*)
write(10,*) "++++++++++++++++++ 1: Force field setup ++++++++++++++++++", &
          &  "+++++++++++++++++++++++++++++"
write(10,*) 
!
!     read in different input files and determine force field 
!     name as well as number of atoms in system
!     options: O = orca, G = gaussian, T = turbomole,C = CP2K
!
if (software .eq. "O") call rdo0  (fname,outname,n,check_coord)
if (software .eq. "T") call rd0   (fname,n)
if (software .eq. "G") call gaurd0(fname,n)
if (software .eq. "C") call rdc0 (fname,n)
if (software .eq. "V") call rdv0 (fname,n)
!
!     allocate all needed arrays (if n > 0)
!
if(n.lt.1) stop 'no atoms!'

!
!     For VASP calculations: compare geometry of second minimum to first one!
!
if (software .eq. "V") then
   if (.not. second_qmdff) then
      allocate(xyz_previous(3,n))
   end if
end if

allocate(at(n),xyz(3,n),g(3,n),nb(20,n),q(n),c6xy(n,n),h(3*n,3*n), &
   &      grep(3,n),cn(n),chir(n),xyz2(3,n),epair(n*(n+1)/2),  &
   &      cring(8,n),ringsize(n),mlist(n),covpair(n*(n+1)/2), &
   &      wbo(n,n),xyz1(3,n),hyb(n),molist(n),freq(3*n))
!
!     after allocation, read in coordinates, hessian, charges and 
!     Wiberg-Mayer bond orders from same files as above
!
if (software .eq. "O") call rdo(.true.,fname,outname,n,xyz,at,check_coord)
if (software .eq. "T") call rd (.true.,fname,n,xyz,at)
if (software .eq. "G") call gaurd(fname,n,at,h,xyz,chir,wbo,check_coord)
if (software .eq. "C") call rdc(.true.,fname,n,xyz,at)
if (software .eq. "V") call rdv(.true.,fname,n,xyz,at,second_qmdff,xyz_previous)
!write(*,*) "summmm",sum(wbo)
!
!     if the system contains more than 100000 (!) atoms abort
!
if (10*n.gt.ndim) stop 'system too large. recompile (ndim)'
!
!     if the system contains more than 300 atoms calculate no EHT
!     theory and bonds
!
if (n.gt.300) bondcalc=.false.
!
!     Not PRECOMPUTE parameters and others for NON-BONDED terms
!     E=Erep + Edisp + Ees
!
write(10,*)
write(10,*) 'Initialization of nonbonded interactions:'
write(10,*)
!
!     copy C6 parameters for D3 formula von hard coded array
!
call copyc6
!
!     set up the D3 part and nonbonded interactions in general
!     (not molecule but element sepecific)
!
call setnonb(scalehb,scalexb,vz,sr42,zab,r0ab)
!
!     set numbers of bond partners by calculating distance 
!     and include inverse damping function
!
call ncoord_qmdff(n,rcov,at,xyz,cn,5000.0d0)
!
!     set up parameters for EHT calculation
!
call setvalel
call setZETAandIP
!
!     interpolate C6 parameters according to D3
!     loop over all atoms in the system (i1) and their atomtypes (at(i1))
!
do i1=1,n
   iz1=at(i1)
   do i2=1,i1
      iz2=at(i2)
      call getc6(5,94,c6ab,maxci,iz1,iz2,cn(i1),cn(i2),c6)  
      c6xy(i2,i1)=c6
      c6xy(i1,i2)=c6
   end do
end do
!
!     get name of hessian file
!
if (index(fname,'.hess').ne.0) then
   hname=fname
   hex=.true.
else if (index(fname,'coord').ne.0) then
   hname='hessian'
   inquire (file=hname,exist=hex)
   if (.not.hex)then
      hname='hessian_driver'
      inquire (file=hname,exist=hex)
   end if
else if (index(fname,'.out').ne.0) then
   hname=fname
   hex=.true.
else if (software.eq."G") then
   hname=fname
   hex=.true.
else if (software .eq. "V") then
   hname=fname
   hex=.true.
end if
if (.not.codetest) call system('rm -f solvent')
!
!     Get the Hirshfeld charges from input and calculate 
!     the needed CM5 charges with them
!
!     (special treatment for HF molecule)
!

if (n.eq.2.and.at(1).eq.1.and.at(2).eq.9) qscal=1.6
if (n.eq.2.and.at(1).eq.9.and.at(2).eq.1) qscal=1.6

!
!     read Hirshfeld charges from charge file
!     only if Gaussian isn´t used!
!
!     Skip the whole process if the check_coord option is activated!
!     In that case, simply set charges to zero..
!
fname2=fname_pre//'.out'
fname3=fname_pre//'.log'
if (software .eq. "V") then
   fname3=fname_pre//'.charges'
end if
if (check_coord .and. software .ne. "G") then
   q=0.0
else 
   if (index(fname_pre,'.out' ).eq.0) then
      call gethirsh(n,chir,ex,fname2,fname3)  !modificated
      ! ex is set true if read in was successfull
   else
      ex=.true.
   end if
!
!     In the case of Gaussian input: charges were already read in
!
   if (software .eq. "G") then
      q=chir
      ex = .true.
   end if
!write(*,*) "qjdh",chir(1:3)
   if (ex) then
!
!     calculate the CM5 chagres from the Hirshfeld charges
!
      call docm5(n,at,.true.,xyz,chir,q)
   else
     ! stop 'no charges!'
   endif
   write(method,'(''CM5*'',F4.2)')qscal
   write(10,*)
   write(10,'('' SCALING CM5 charges by '',F5.2)')qscal
   write(10,*)'and as such written to FF file'
   write(10,*)
end if
!
!    Scale the obtained CM5 charges!
!

q = q * qscal + (1.-qscal) *chrg/dble(n)

!
!    Read the restrain potential (whoever that needs...)
!

call rdfix(n)

!
!    Read the Wiberg-Mayer-Bondorders (WBO) matrix from file
!    --> distinction between gaussian and orca!
!
! ---> TEST TEST
!if (index(fname,'.out' ).eq.0) then
!   call rdwbo(n,wbo,ex,fname2)
!else if (index(fname,".log") .eq. 0) then
!   stop "hhfgg"
!else
!   ex=.false.
!end if

ex=.true.
if (software .eq. "O") then
   call rdwbo(n,wbo,ex,fname2)
else if (software .eq. "G") then
!   call rdgwbo(n,wbo,ex,fname3) 
!   ex=.true. 
!else 
!   ex=.false.
end if
!
!     If no bond orders were calculated (i.e. all are zero), 
!     activate the EHT calculation of them!
!
if (sum(wbo) .le. 0.1d0 .and. software .ne. "C") then
   ex=.false.

end if
!
!     test 29.11.2018 : commented out wbo calculation No. 2 (was already done above)
!

!
!     status 15.03.17
!
!     If that bondorders can´t be found, calculate them based on the 
!     EHT calculation! (for example if no orca output is present)
!
!     The calculated values are similar (within 5-10% deviation) to the 
!     reference data!!!
!
!ex=.false.

if (.not.ex) then
   write(10,*)
   write(10,*) 'No Wiberg Mayer bond orders found in reference output!!!!'
   write(10,*) 'An Extended Hückel calculation will be done to determine them.'
   write(10,*) 'In Most cases these bond orders should be reliable.'
   write(10,*)
!
!     Define basis functions for the EHT hamiltonian
!
   call basis0_eht(n,at,nel,nbf)
   call basis_eht (n,at,nbf,okbas)



   if (.not.okbas) stop 'TB Hamiltonian incomplete'
   nel=nel-chrg
!
!     Do the EHT calculation 
! 
   call ehtfull(n,at,xyz,q,nel,nbf,el,wbo,1.0d0)
   write(10,*) "Wiberg-Mayer bond orders were calculated "
   write(10,*) "   based on EHT hamiltonian!"
end if
!
!     Determine the nearest neighbors of all atoms
!

call setmetal
call ff_neighbor(n,xyz,wbo,at,nb,hyb,mclust,wbocut_in)

!
!     Write informations about charges, momenta and neighbors
!

dipol2=0
write(10,*) "Summary for reference molecule:"
write(10,*)'---------------------------------------------------------', &
  &   '-----------------------------'
write(10,*)'  #   Z           coordinates           CN   q', &
  & '   hyb/lin/metbond   neighbors'
write(10,*)'---------------------------------------------------------', &
  &   '-----------------------------'
do i=1,n
   write(10,'(2i4,3f10.5,f6.2,f7.3,3x,3i2,3x,20i4)')  &
     &  i,at(i),xyz(1:3,i),cn(i),q(i), &
     &  hyb(i),nb(19,i),nb(18,i), &
     &  nb(1:nb(20,i),i)
   dipol2(1:3)=dipol2(1:3)+xyz(1:3,i)*q(i)
end do
write(10,*)'---------------------------------------------------------', &
  &   '-----------------------------'
write(10,*)

d2=sqrt(dipol2(1)**2+dipol2(2)**2+dipol2(3)**2)
write(10,'('' monopole moment for charges used (au/Debye): '', &
  &  2f10.5)') d2,d2*2.5418
!
!     Get hessian from file 
!
if (.not.hex) then
   write(*,*) 'no Hessian found'
   doopt=.false.
   domd =.false.
   nbond=0
   nangl=0
   ntors=0
   goto 1111
else
   write(10,*) 
   write(10,*) 'Contains the structure several fragments?'
!   if (hname.eq.'hessian'.or.hname.eq.'hessian_driver') then
!     The gaussian hessian
!      call rdhess(n*3,h,hname)
!   else if (index(hname,'.hess').ne.0) then
!     The orca hessian
!      call rdohess(n*3,h,hname)
!   end if
    if (software .eq. "O") then
       call rdohess(n*3,h,hname)
    else if (software .eq. "C") then
       call rdchess(n,n*3,h,hname)
!
!     for gaussian small change by JS: first, convert output 
!     to fchk-file!
!     now, a single calculation for geoopt and frequency etc. is possible!
!
    else if (software .eq. "G") then
       finished=0
       call rdghess(n*3,h,hname,fname_pre,finished)
!
!
!     Read in the Hessian from a VASP calculation: combine the Hessian from
!     the gradient vectors stated in the OUTCAR file numerically!
!
    else if (software .eq. "V") then
       call rdvhess(n*3,h,hname)
       if (vasp_hessian_sel) then
          write(*,*)
          write(*,*) "WARNING: you have used selective dynamics for the VASP reference!"
          write(*,*) " qmdffgen will nevertheless optimize all force constants."
          write(*,*) " This might lead to some additional errors!"
       end if
    end if    

end if

!
!     determine bond distance list, 5 or more covalent bonds must
!     be in between to atoms to include their NCI terms
!     otherwise they get a bond stretch term 
!     in fragmentation metal systems are treated specially      
!

call fragmentation(n,at,nb,molist)
call ff_bond(n,at,nb,xyz,covpair,q,wbo,cn,molist)

!
!     Define the QMDFF! (Bonds, angles, torsions etc.)
!
call ff_set(n,xyz,xyz2,at,nb,vdw,q,scalehb,scalexb,covpair, &
     &  r0ab,zab,r094,sr42,c6xy,echo,wbo,hyb, &
     & cring,ringsize,mlist,mtors)
!
!     Write force field information output
!

!if (echo) then
   write(10,*) '-----------------------------------------------------------'
   write(10,*) 'Force field definition (bonded part): '
   write(10,*) '-----------------------------------------------------------'
   do i=1,nbond
      write(10,'('' bond ij,R:'',2i4,2F8.4)') &
        &  bond(1,i),bond(2,i),vbond(1,i),vbond(3,i)
   end do
   write(10,*) '-----------------------------------------------------------'
   do i=1,nangl
      write(10,'('' bend  ijk,theta:'',3i4,f8.2)') &
        &  angl(1,i),angl(2,i),angl(3,i),vangl(1,i)*180./pi
   end do
   write(10,*) '-----------------------------------------------------------'
   do i=1,ntors
      write(10,'('' tors ijkl,type,phi  :'',6i4,3f8.2)') &
        &  tors(1:6,i),vtors(1,i)*180./pi,vtors(3,i)
      if (tors(5,i).gt.1) then
         mm=3
         do m=1,ntterm
            write(10,'(i4,10f12.6)') m,vtors(2,i),vtors(mm,i), &
              &  vtors(mm+1,i),vtors(mm+2,i)
            mm=mm+3
         end do
      end if
   end do
   write(10,*) '-----------------------------------------------------------'
!end if
!
!     Calculate total number of terms!
!

nvar=nbond+nangl+ntors
write(10,*) 
write(10,*) "Number of terms included into the force field:"
write(10,'(A,I4)') " * Bonded terms total: ",nvar
write(10,'(A,I4)') "  - Number of bonds:",nbond
write(10,'(A,I4)') "  - Number of angles:",nangl
write(10,'(A,I4)') "  - Number of torsions:",ntors
write(10,'(A,I4)') " * Non-covalent terms total: ",nnci
write(10,*) 

!
!     If the CHECK_COORD mode is activated, skip the Hessian read in and 
!     jump directy to the end!
!
if (check_coord) then
   goto 1111
end if

!     If no variable is defined (fatal...)
if (nvar.eq.0) goto 1999

!
!     Calculate fractions of noncovalent interaction (NCI) 
!     energies and write them out
!

e=0
grep=0
call ff_nonb(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,e,grep)
write(10,*) 'Intramolecular non-bonding energies:'
write(10,*) 'E(Non-covalent)           (Re)=',e
e=0
call ff_hb(n,at,xyz,e,grep)
call ff_nonbe(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,epair,enci)
write(10,*) 'E(hydrogen/halogen bond)  (Re)=',e
write(10,*) 'Enci(electrostatics)      (Re)=',enci(2)
write(10,*) 'Enci(dispersion)          (Re)=',enci(1)
write(10,*) 'Enci(repulsion)           (Re)=',enci(3)

if(echo) stop
!
!     Morst important part:
!     FIT OF THE QMDFF HESSIAN TO REFERENCE!!!
!
write(10,*) 
write(10,*) "++++++++++++++++++ 2: Force field parameter fit ++++++++++", &
          &  "+++++++++++++++++++++++++++++"
call hfit(n,nb,at,xyz,q,r0ab,zab,r094,sr42,c6xy,h,.false.,hneglect)

!
!     Correction: Remove small parts of the potential and 
!     fit again...
!
write(10,*) 
write(10,*) 'removing small force constants from force field ...'
write(10,*) 
call ff_clean(n,xyz,at,refit)
!     If some FF-parts were removed..
i=nbond+nangl+ntors
if (i.ne.nvar.or.refit) then
   write(10,*)
   write(10,*) "Do a second Levenberg-Marquardt fit with reduced hessian: "
   call hfit(n,nb,at,xyz,q,r0ab,zab, &
      &  r094,sr42,c6xy,h,.false.,hneglect)
   write(10,*) 'removing small force constants from force field ...'
   write(10,*) 
   call ff_clean(n,xyz,at,refit)
end if
! 
!     Modify FF for single bonds for which
!     Extended Hückel (also named TB) failed 
!     and scale also metal bendings
!

call ff_mod(echo,n,xyz,q,at,nb,wbo,hyb,cring,ringsize)

!     Test call (?)
!     calculate fraction of bonded energy
999 continue
call ff_eg(n,at,xyz,e,g)
write(10,*)
write(10,*) "====================="
write(10,*) "Energy and gradient"
write(10,*) "====================="
write(10,*) 

write(10,*) 'E(Re,bonded only)=',e
write(10,*) 'G(Re,bonded only)=',sum(abs(g))
!
!     Compute Force Field hessian 
!     and compare it with reference hessian!
!
if (vasp_hessian_sel) then
  write(10,*) 
  write(10,*) "WARNING: you have used selective dynamics for the VASP reference!"
  write(10,*) " qmdffgen will nevertheless optimize all force constants."
  write(10,*) " This might lead to some additional errors!"
end if

call hess(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy, &
  &   val,.true.,hname,h,.false.,freq,dist, &
  &   hesscalc,numhess,software,length,fname_pre)
!
!     Calculate energy and gradient of the new force field! 
!

1999  write(10,*)
write(10,*) "++++++++++++++++++ 3: Force field properties +++++++++++++", &
          &  "+++++++++++++++++++++++++++++"
write(10,*) 
!if (echo) then
   write(10,*)
   write(10,*)'======================='
   write(10,*)'Further investigations:'
   write(10,*)'======================='
!endif
write(10,*)
call ff_eg(n,at,xyz,e,g)
call ff_nonb(n,at,xyz,q,r0ab,zab,r094,sr42,c6xy,e,g)
call ff_hb(n,at,xyz,e,g)
call moment(n,xyz,g,f1,f2)
gnorm=sqrt(sum(g**2))
write(10,*) "Energy and gradient for the reference system:"
write(10,*) 'total E(Re,input)=',e
write(10,*) 'total G(Re,input)=',gnorm
write(10,*) 'total E(Re,input) per atom=',e/n
write(10,*) 'total G(Re,input) per atom=',gnorm/n
write(10,*) '(should be < 0.01)'
!     (Special format for writeout of array elements!)
write(10,50) f1(1:3),f2(1:3)
50 format(/' Fxyz = ',2(d12.3,','),d12.3, &
  & /, ' Mxyz = ',2(d12.3,','),d12.3,/)

!
!     Calculate dissociation energy (of all existing bonds):
!     Move parts away from each other and calculate energy
!
e0=e
e0input=e
xyz2 = xyz * 10000.
call ff_e(n,at,xyz2,q,r0ab,zab,r094,sr42,c6xy,einf)
eat=einf-e0

write(10,'('' Total De (dissociation energy) molecule -> atoms &
  &  (Eh/kcal)'')')
write(10,'(2F12.3)') eat,eat*627.509541
!
!     Calculate several specific thermochemical properties etc.
!
if (.not. doopt) then
   call thermo(n,at,xyz,freq,h298)
   call heat(n,at,eat,h298,freq)
end if

!
!     If a van-der-waals complex is present, calculate charges 
!     of both parts
!
if (vdw)then
   sum1=0
   sum2=0
!
!     make a geometry with separated monomers     
!    
   do i=1,n
      if (mlist(i).eq.1)then
         sum1=sum1+q(i)
         xyz2(1:3,i)=xyz(1:3,i)
      else
         sum2=sum2+q(i)
         xyz2(1:3,i)=xyz(1:3,i)+1.d+10
      end if
   end do
!   call ff_e(n,at,xyz2,q,r0ab,zab,r094,sr42,c6xy,e)
   write(10,*)
   write(10,*) 'separated monomers:'
   write(10,*) 
   write(10,*) 'charge on fragment 1 =',sum1
   write(10,*) 'charge on fragment 2 =',sum2
   write(10,*) 'De (vdW complex, kcal) =',(e-e0input)*627.509541
   write(10,*)
end if

!
!     Optional: analyze different energy parts..
!
if (bondcalc) then
   call bonde(n,at,xyz,covpair,nb,cring,ringsize)
end if
if (analyse)then
!   call rd(.true.,'coord1',n,xyz1,at)   
   call ff_eg(n,at,xyz1,e,g)
   call ff_nonb(n,at,xyz1,q,r0ab,zab,r094,sr42,c6xy,e,g)
   call ff_hb(n,at,xyz1,e,g)
   write(10,*) 'E (coord1) =',e
   write(10,*) 'G (coord1) =',sqrt(sum(g**2))
   call moment(n,xyz1,g,f1,f2)
   write(6,50) f1(1:3),f2(1:3)

   inquire(file='coord2',exist=ex)
   if (ex) then
      call rd(.true.,'coord2',n,xyz2,at)
   else
      call splitmol(n,at,xyz1,wbo,mlist,mol1,0,0)
      if(mol1.ne.n) then
         write(10,*)
         write(10,*) '================='                
         write(10,*) 'found vdW complex'                
         write(10,*) '================='      
         write(10,*)           
!     make a geometry with separated monomers         
         do i=1,n
            if (mlist(i).eq.1) then
               xyz2(1:3,i)=xyz1(1:3,i)
            else
               xyz2(1:3,i)=xyz1(1:3,i)+10000.
            end if
         end do
      end if
   end if
   call ff_eg(n,at,xyz2,e,g)
   call ff_nonb(n,at,xyz2,q,r0ab,zab,r094,sr42,c6xy,e,g)
   call ff_hb(n,at,xyz2,e,g)
   write(10,*) 'E (coord2) =',e
   write(10,*) 'G (coord2) =',sqrt(sum(g**2))
   call moment(n,xyz2,g,f1,f2)
   write(6,50) f1(1:3),f2(1:3)
   call ff_anal(n,at,xyz1,xyz2,q,r0ab,zab,r094,sr42,c6xy)
endif
!
!     At least: Write the QMDFF force field to file!
!
!     entry point if Hessian does not exist
!     continue: dummy-statement for jump numbers!
! 
1111 continue

if (.not.rds) then
   solvname=fname_pre // '.qmdff'
   write(10,*)
   write(10,*) 'writing <QMDFF> file (',fname_pre,'.qmdff).'
   open(unit=1,file=solvname)
   izero=0
   one=1.0
   write(1,'(i4,2F10.4)') n,one,dens
   write(1,'(a)') trim(method)
   do i=1,n
      write(1,'(i5,4F18.12,i5)')at(i),xyz(1:3,i),q(i),izero
   end do
   write(1,'(6i8)') nbond,nangl,ntors,nhb,nnci,nbond12
   do i=1,nbond
      write(1,'(2i5,3F18.12)')bond(1,i),bond(2,i),vbond(1:3,i)
   end do
   do i=1,nangl
      write(1,'(3i5,5F18.12)')angl(1,i),angl(2,i),angl(3,i),vangl(1:2,i)
   end do
   do i=1,ntors
      write(1,'(4i4,2i2,2F16.12,10(2F4.1,E16.8))') &
        &   tors(1:6,i),vtors(1:2,i),(vtors(3*(j-1)+3,i), &
        &   vtors(3*(j-1)+4,i)/pi,vtors(3*(j-1)+5,i),j=1,tors(5,i))
   end do
   if (nhb.gt.0) write(1,'(6(3i4,2x))')hb(1:3,1:nhb)
   write(1,'(8(3i4,2x))')nci(1:3,1:nnci)
   close(1)
   write(10,*)
   write(10,*) "QMDFF generation finished!!!"
   write(10,*) 
end if

return
end subroutine main_gen


