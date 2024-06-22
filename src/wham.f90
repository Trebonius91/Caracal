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
!     Subroutine wham: Calculates the free energy surface by using 
!     all sampled histograms from the umbrella samplings with the 
!     Weigted Histogram Analysis Method (WHAM)
!
!     part of EVB
!
!
!
!     REFERENCE: B. ROUX, COMPUTER PHYSICS COMMUNICATIONS VOL. 91, PP. 275-282 (1995)
!     THE CODE ASSUMES THAT THE TRAJECTORIES 1.7.xyz, 1.8.xyz,....., 4.5.xyz ARE PRESENT
!     RAYMOND ATTA-FYNN
!     EMSL, PACIFIC NORTHWEST NATION LAB
!     ORIGINALLY WRITTEN ON OCTOBER 23, 2009
!     MODIFIED ON JUNE 1, 2011 WITH COMMENTS
!
!
subroutine wham(nbins,n_all,xi_min,xi_max,bin_size,xi_wins,p_biased,pmf)
use general
use evb_mod
implicit none
integer::NW, I, J, K, NSAMP,nbins,n_all
!     NMAX1/n_samplings IS THE MAX WINDOW DIMENSION
!     NMAX2/nbins IS THE MAX NUMBER OF HISTOGRAM BINS
real(kind=8)::DXI, XI_I, KB, TEMP, BETA2,bin_size
real(kind=8)::xi_min,xi_max  ! the borders of the sampling coordinate
real(kind=8)::tmp   ! sum of all histogram values
real(kind=8)::L,DT,TCUT,KS
real(kind=8)::NI(n_all) ! SAMPLED SNAPSHOTS IN A GIVEN WINDOW
real(kind=8)::W(n_all,nbins) ! BIASING POTENTIAL
real(kind=8)::P_BIASED(n_all,nbins) ! BIASED DISTRIBUTION
real(kind=8)::P_UNBIASED(nbins) ! OPTIMAL ESTIMATE OF UNBIASED DISTRIBUTION
real(kind=8)::FI(n_all) ! FI = EXP(-fi*BETA)
real(kind=8)::pmf(nbins)  ! the resulting free energy surface!
real(kind=8)::HIST(nbins) ! DISTRIBUTION (HISTOGRAM) OF A GIVEN SIMULATION
real(kind=8)::EPS,TOL,FI_OLD,PMIN,TMP1,TMP2
real(kind=8)::xi_wins(n_all)  ! values of all umbrella window coordinates
real(kind=8)::V_res(nbins),FAC,NUMERATOR,DENOMINATOR,FI_NEW
real(kind=8)::eps_old  ! for smaller output of convergence parameter: only each order of magnitude
integer::stepinfo  ! number of steps needed for convergence (only for information..)
CHARACTER(len=3)::F*3 ! THIS WAS DEFINED FOR FILE READING PURPOSES
!*****************************EDIT THIS SECTION********************************
xi_min=xi_min                 ! MIN REACTION COORDINATE IN ANGSTROM
xi_max=xi_max                 ! MAX REACTION COORDINATE IN ANGSTROM
DXI=bin_size                   ! BIN WIDTH
nbins=nbins                   ! TOTAL NUMBER OF BINS
L=11.4D0                      ! SIMULATION BOX LENGTH IN ANGSTROM (only for PBC..)
DT=30.D0                      ! TIME INTERVAL FOR SAVING CONFIGURATIONS IN A.U.
DT=DT/(41.3413733D+3)         ! CONVERT TIME STEP FROM AU TO PICOSECONDS
TCUT =1.D0                    ! EQUILIBRATION TIME IN PICOSECONDS
NSAMP =2                      ! SAMPLE EVERY OTHER POINT.
                              ! TO SAMPLE AT A LOWER RATE INCREASE NSAMP.
KB=0.316679D-5                ! BOLTZMANN’S CONSTANT IN A.U./KELVIN
TEMP =kelvin                  ! SIMULATION TEMPERATURE IN KELVIN
BETA2 =1.D0/KB/TEMP            ! INVERSE Kb*T
FAC=627.509D0                 ! CONVERSION FACTOR FROM A.U. TO KCAL/MOL
NW=n_all                ! TOTAL NUMBER OF WINDOWS
TOL=1.D-5                     ! SCF TOLERANCE
eps_old=1.0d0                 ! for writing out of EPS parameter
stepinfo=1                    ! number of steps needed for convergence
!*****************************END EDITING***************************************
!*********************COMPUTATION BEGINS HERE**********************************
!write(*,*) "Begin weighting of histograms with the WHAM method to"
!write(*,*) " generate the PMF surface for the reaction!"
DO I = 1, NW ! LOOP OVER WINDOWS
   XI_I = xi_wins(i) ! COMPUTE CENTRAL COORDINATE OF WINDOW I
!   WRITE(F,'(F3.1)')XI_I ! CHARACTER F TAKES VALUES 1.7, 1.8, 1.9,......,4.5
!     XYZ FILES NAMED 1.7.xyz, 1.8.xyz,......, 4.5.xyz WILL BE OPENED SUCCESSIVELY
!     INSIDE THE SUBROUTINE PROB BELOW
!   CALL PROB(XI_MIN,DXI,N,XI_I,L,DT,TCUT,KS,NSAMP,F,K,HIST,V) ! CALL SUBROUTINE
!   NI(I)=DBLE(K)
!
!     Compute normalization factor for each single distribution
!        
   hist=p_biased(i,:)
   tmp=0
   DO J = 1, nbins
      tmp = tmp + HIST(j)
   END DO
   DO J = 1, nbins ! LOOP OVER HISTORGRAM POINTS IN WINDOW I and calculate restrain for each bin!
      TMP2 = xi_min+bin_size*j    ! REACTION COORDINATE IN ANGSTROM in this bin
      HIST(j) = HIST(j)/TMP             ! NORMALIZED PROBALITY AT TMP2
      v_res(j) = k_force(I)*((TMP2-XI_I))**2          ! BIASING WINDOW POTENTIAL
!      v_res(j) = v_res(j)*FAC**2                ! CONVERT POTENTIAL TO A.U.
      P_BIASED(I,J)=HIST(J) ! BIASED DISTRIBUTION OF WINDOW I AT COORDINATE XI_J
      
      W(I,J)=DEXP(-BETA2*V_res(J)) ! RESTRAINING POTENTIAL OF WINDOW I AT COORDINATE XI_J
!      write(79,*) w(i,j),(tmp2-xi_i)**2,v_res(j)
   END DO
   ni(i)=tmp
END DO
!INITIAL GUESS OF UNITY FOR FI = DEXP(-fi*BETA)
DO I = 1, NW
   FI(I) = 1.D0
END DO
100 CONTINUE

DO J = 1, nbins
   NUMERATOR= 0.D0
   DENOMINATOR = 0.D0
!
!     COMPUTING THE DENOMINATOR OF EQUATION (1)
!
   DO I = 1, NW
      DENOMINATOR = DENOMINATOR + NI(I)*W(I,J)/FI(I)
   END DO
!
!     COMPUTING THE NUMERATOR OF EQUATION (1)
   DO I = 1, NW
      NUMERATOR = NUMERATOR + NI(I)*P_BIASED(I,J)
   END DO

   TMP1=NUMERATOR/DENOMINATOR
 !  write(*,*) TMP1,numerator,denominator
 !  stop
!
!     CONVERT ZERO PROBABILITY TO VERY SMALL POSITIVE NUMBER
!
   IF (TMP1.EQ.0.D0) TMP1=1.D-15
   P_UNBIASED(J)=TMP1
END DO
!
!     COMPUTE NEW FI BASED AND THE OLD USE EPS=SUM{(1-FI_NEW/FI_OLD)**2}
!     AS THE CONVERGENCE PARAMETER
!

EPS = 0.D0
DO I = 1, NW
   FI_OLD = FI(I)
   FI_NEW = 0.D0
   DO J = 1, nbins
      FI_NEW = FI_NEW + W(I,J)*P_UNBIASED(J)
   END DO
   FI(I)=FI_NEW
   EPS = EPS + (1.D0 - FI_NEW/FI_OLD)**2
END DO
!
!     write new residual only if the error has reached a new order of 
!     magnitude! (smaller)
!
if (int(log10(eps)) .ne. int(log10(eps_old))) then
   WRITE(15,*) "Steps:",stepinfo,"remaining error:", EPS
   WRITE(*,*) "Steps:",stepinfo,"remaining error:", EPS
end if
stepinfo=stepinfo+1
eps_old=EPS
IF (EPS.GT.TOL) GOTO 100

WRITE(15,*) "WHAM has converged!!"
WRITE(*,*) "WHAM has converged!!"
!
!     FIND FREE ENERGY SHIFT
!
PMIN = -1.D+8
DO J = 1, nbins
   PMIN=MAX(PMIN,(1.D0/BETA2)*DLOG(P_UNBIASED(J)))
END DO
!     PRINT PMF
OPEN(10,FILE='pmf_wham.dat',STATUS="UNKNOWN")
DO J = 1, nbins
   TMP1= XI_MIN + DXI*(DBLE(J)-0.5D0) ! REACTION COORDINATE IN ANGSTROM
   TMP2=(-(1.D0/BETA2)*DLOG(P_UNBIASED(J))+PMIN)*FAC*bohr ! PMF IN KCAL/MOL (scaled to bohr!)
   pmf(j)=TMP2  ! store the surface for later using   
   WRITE(10,1000)TMP1,TMP2
END DO
1000 FORMAT(2F12.6)
CLOSE(10)

return
end subroutine wham

