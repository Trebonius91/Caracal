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
!     ##################################################################
!     ##                                                              ##
!     ##  module evb_mod  --  parameters for evb_qmdff-calculations   ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "evb_mod" contains all important global informations that are
!     used in the evb-qmdff-parts: qmdffgen,evbopt,dynamic etc.
!

module evb_mod
implicit none
save

integer::nqmdff  !the number of used QMDFF´s
integer n_one,natoms
real(kind=8),allocatable :: xyz (:,:)
real(kind=8),allocatable :: q   (:)
real(kind=8),allocatable :: g_one (:,:)
real(kind=8),allocatable :: c6xy(:,:)
real(kind=8),allocatable :: cn  (:)
integer,allocatable :: at (:)
integer,allocatable :: at2(:)
integer,allocatable :: imass(:)
!integer,allocatable :: molnum(:)  ! for solvent box simulations
!real(kind=8)::dens  ! density (i.e. number of molecules in box)

integer n_two
real(kind=8) ,allocatable :: xyz_two (:,:)
real(kind=8) ,allocatable :: q_two(:)
real(kind=8) ,allocatable :: g_two(:,:)
real(kind=8) ,allocatable :: c6xy_two(:,:)
real(kind=8) ,allocatable :: cn_two(:)
integer,allocatable :: at_two(:)
integer,allocatable :: imass_two(:)

integer n_three
real(kind=8) ,allocatable :: xyz_three (:,:)
real(kind=8) ,allocatable :: q_three(:)
real(kind=8) ,allocatable :: g_three(:,:)
real(kind=8) ,allocatable :: c6xy_three(:,:)
real(kind=8) ,allocatable :: cn_three(:)
integer,allocatable :: at_three(:)
integer,allocatable :: imass_three(:)

real(kind=8)::offa,offb,offc,offd,offe,offf,offg,offh,offi,offj,&
        &offk,offl,offm,offn
! The shifted single-point energy of the QMDFFs
real(kind=8)::E_zero1,E_zero2,E_zero3  

real(kind=8),dimension(94,94)::r0ab,zab,sr42,r094_mod
real(kind=8),dimension(12) :: start_params
real(kind=8),dimension(14) :: evb3_params
!integer::mode,off_choice
character(len=40)::mode1,off_basis,calc
character(len=40)::orca_inp
character(len=:),allocatable::orca_egrad
! the simple dE coupling term
logical::evb_de
! decides, if evb_qmdff or orca_calculation
logical::orca_calculation
! shall plumed be used or not!
logical::use_plumed
! if you want to control that the energy is conserved
logical::energycon
! separate printout of covalent and noncovalent energies for a single QMDFF
logical::energysplit
real(kind=8)::e_cov_split,e_noncov_split
! maximal iterations for levmarq
integer::optsteps
! start-value for levmarq parameter between newton and CG and the alternation
real(kind=8)::lm_par
real(kind=8)::lm_par_change
! threshold-parameter for convergence in levmarq
real(kind=8)::lm_threshold
! numercial differentation-step in levmarq for 3x3-EVB
real(kind=8)::diff_step
! for the evb_dq-gradient
real(kind=8),dimension(:),allocatable::ts_coordinates,ts_coordinates2
logical::use_dq,evb_dq
! for mopac reference-calculations: number of commandlines after coordsection
integer::mopac_maxline
! numerical gradient for 3x3-evb for comparison, sets 1-3 elements not zero 
logical num_grad, full
real(kind=8)::num_grad_step   ! stepsize for numerical gradients
!  parameters for DG-EVB with more than one reference point
logical::dg_evb
! Plumed Inputfile
character(len=80)::plumed_input
integer::dg_evb_points ! number of reference points for DG-EVB calculation
integer::nat6 ! NUMBER OF INTERNAL COORDINATES
integer::num_coord ! number of internal coordinates defined by user
real(kind=8),dimension(:),allocatable::all_ens,allf1_ens,allf2_ens,f_vec,b_vec
real(kind=8),dimension(:,:),allocatable::all_xyz,d_mat,all_int,point_int
real(kind=8),dimension(:,:),allocatable::all_grad,allf1_grad,allf2_grad
real(kind=8),dimension(:,:,:),allocatable::all_hess,allf1_hess,allf2_hess
logical::dg_evb_opt_multi,diff_ana ! if the alpha values should be optimized separately or at once
integer::maxstart !maximal number of local starts for optimization 
integer::dg_mode       ! for DG-EVB (1,2 or 3)
integer::ref_input_unit  ! unit for the ref.input file in DG-EVB calculations
real(kind=8),dimension(:),allocatable::alph_opt,b_opt
! start interval for random number
real(kind=8)::lower_bond,upper_bond
! if the ordering of internal coordinates should be read in manually
logical::read_coord
! if the internal coordinates shall be defined as distance matrix!
logical::dist_matrix
! global array of internal coordinate definitions (D:1; Z:4)
integer,dimension(:,:),allocatable::coord_nr,coord_def
integer,dimension(:,:),allocatable::bonds,angles,ooplanes,diheds
! global integers for numbers of single coordinates (bonds,angles,ooplane,diheds)
integer::nbonds,nangles,noops,ndiheds
! for normal mode analysis
logical::print_nm
! for geometry optimization
integer::geomax
! global parameter for additional info output during programs
logical::more_info
! if only a fragment of the structure shall be considered for calculation
logical::calc_frag
integer,dimension(:),allocatable::indi
integer,dimension(:),allocatable::indi3
integer::nats
! for IRC optimizations
character(len=50):: irc_method
real(kind=8)::irc_steplen
real(kind=8)::irc_eulerlen
integer::irc_maxstep
! test case: mueller brown surface
! name of the DG-EVB reference file 
character(len=80)::dg_ref_file
! if the DG-EVB coefficient matrix D shall be calculated numerically
logical::dg_mat_num
! minimum value below that distant Gaussians are neglected
real(kind=8)::g_thres
! if for DG-EVB mode 3 two alphas per reference point shall be used 
logical::double_alpha
integer::add_alph
! for the new RP-EVB coupling method
logical::rp_evb
integer::rp_evb_points  ! number of reference points with energies, gradients and hessians
integer::rp_irc_points  ! number of structures with energies along the IRC
integer::rp_spl2_dim   ! total number of reference informations per gradient/hessian point
real(kind=8)::lambda  ! the inverse average disctance between two reference frames
real(kind=8)::pareta  ! shape parameter for the damping of V12
real(kind=8)::interp_tol  ! how small the interval for s-determin. shall become
real(kind=8)::rp_ana_step  ! for analytical gradient: step of internal coordinate elongation
real(kind=8),allocatable::rp_spl_ref(:,:) ! reference information for spline interpolation
real(kind=8),allocatable::rp_spl_d2(:,:) ! second derivatives for spline interpolation
real(kind=8),allocatable::rp_spl_s(:)  ! the parameter reference values for spline interp.
real(kind=8),allocatable::rp_spl_ref2(:,:) ! same as above, but for gradients/hessians
real(kind=8),allocatable::rp_spl_d22(:,:) ! same as above, but for gradients/hessians
real(kind=8),allocatable::rp_spl_s2(:) ! same as above, but for gradients/hessians
real(kind=8),allocatable::rp_point_s(:) ! s values of the points with gradients and hessians
real(kind=8),allocatable::rp_point_v12(:)  ! coupling strengths on these points
real(kind=8),dimension(:,:),allocatable::dv12 ! the coupling gradients at the rp_evb points
real(kind=8),dimension(:,:,:),allocatable::d2v12  ! the coupling hessians at the rp_evb points
integer::qmdff_order    ! if QMDFF1=educt, QMDFF2=product or the other way
real(kind=8)::s_transfer  ! the actual value of the reaction coordinate s
real(kind=8)::s_ts    ! position of the TS along the path (in terms of the s-parameter)
real(kind=8)::s_tot    ! the total length of the IRC for rescaling
real(kind=8)::par_epsi  ! slope of the correction potential parabola
real(kind=8)::pre_exp   ! exponential prefactor of the gaussian damping for RP-EVB
real(kind=8)::corr_max   ! maximal value of the parabola prefactor for the QMDFF corrections
real(kind=8)::deltp   !  half size of the whole interval for which QMDFF shall be corrected
real(kind=8)::ddeltp  !  shift of the correction interval borders for the QMDFF hessians
logical::s_bord_man  ! if the borders along the IRC shall be read in manually
real(kind=8)::rp_mid_tot,rp_mid_trans   ! NEW 04.09.2018: only energy region near the TS
real(kind=8)::trans_r_lo,trans_r_hi  ! transition borders right to the ts (s values)
real(kind=8)::trans_l_lo,trans_l_hi   ! transition borders left to the ts (s values)
real(kind=8),dimension(:,:),allocatable::xyz_r_hi,xyz_l_lo  ! xyz structures of RP-EVB borders
real(kind=8),dimension(:),allocatable::int_r_hi,int_l_lo  ! int proj structures of RP-EVB borders
logical::shift_man   !  turn off automatical optimization of QMDFF minimum shift energies
logical::int_grad_plot  ! option to activate printing of internal gradient components to file..
logical::int_coord_plot  ! option to activate printing of internal coordinate components
logical::num_wilson  ! option, to numerically calculate the wilson b matrix for coord. transf.
logical::use_internals ! option to use internal coordinates for geometry optimizations
logical::spline_all  ! option to activate full IRC scan for Xi value for RP-EVB points
real(kind=8)::path_dist_limit  ! if the z-value is larger than this limit, QMDFF1/2 will be assumed
logical::no_evb  ! if the RP-EVB part shall be deactivated and only direct interpolation/QMDFF be used
real(kind=8)::iq_l_lo,iq_l_hi,iq_r_lo,iq_r_hi  ! borders between direct interpolation and QMDFFs
! for Mechanochemistry
logical::add_force
integer::force1_at,force2_at
real(kind=8)::force1_k,force2_k
real(kind=8)::force1_v(3),force2_v(3)
logical::afm_run
integer::afm_fix_at,afm_move_at
real(kind=8)::afm_move_dist
real(kind=8)::afm_move_v(3)
real(kind=8)::afm_move_first(3)
real(kind=8)::afm_k
integer::afm_steps
! for the new RPMD program
integer::nbeads   ! the number of ring polymer beads to be sampled
integer::npaths  ! number of equivalent reaction paths
real(kind=8),dimension(:,:,:),allocatable::p_i,q_i  ! positions and momenta of the system
logical::use_rpmd  ! if the rpmd.x program is used, allow also 1 internal coordinate..
integer,allocatable::rpmd_atoms(:) ! if not all atoms shall be treated with RPMD, define the active ones
real(kind=8),dimension(:,:),allocatable::ts_ref,ed_ref
real(kind=8),dimension(:),allocatable::int_ed,int_ts,xi_scale ! scale int. coordinates to 0..1
real(kind=8),dimension(:),allocatable::form_ref,break_ref  ! reference values for the TS bonds
real(kind=8),dimension(:),allocatable::form_ed,break_ed ! reference values of educts bonds (unimol.)
real(kind=8)::k_force ! umbrella force constant
real(kind=8),allocatable::struc_equi(:,:,:)
real(kind=8)::umbr_lo,umbr_hi ! bonds for umbrella samplings
real(kind=8)::umbr_dist ! distance between single umbrella samplings
character(len=20)::umbr_type  ! type of coordinate system to be used for umbrella sampling
integer::form_num,break_num  ! number of forming and breaking bonds for the reac.type
integer,dimension(:,:),allocatable::bond_form,bond_break  ! forming and breaking bonds
real(kind=8),allocatable::r_refs(:,:)  ! the needed educt-educt distances for the pre-equilibration
integer::shift_atom,shift_coord  ! for the ATOM_SHIFT coordinate 
real(kind=8)::shift_lo,shift_hi  ! for the ATOM_SHIFT coordinate 
logical::fix_atoms  ! if some atoms shall be fixed in dynamics calculations
integer,allocatable::fix_list(:) ! file with list of fixed atoms
integer::fix_num  ! number of fixed atoms
integer::umbr_traj ! number of umbrella trajectories per Xi value
real(kind=8)::t_avg ! TEST for temperature equilibrium
integer::n_samplings,n_over  ! number of sampling windows
logical::writestat ! decides if MD statistics shall be written out
character(len=20)::pmf_method  ! if WHAM or umbrella integration shall be used!
real(kind=8),allocatable::average(:),variance(:)   ! for umbrella integration!
integer::sum_eds   ! total number of educt molcules (1 up to 4)
integer,allocatable::at_ed(:,:)  ! atom numbers of all educts (up to four for ADDITION4)
real(kind=8)::R_inf   ! reference distance of both/all educts
integer,allocatable::n_ed(:) ! number of atoms in the educts (up to four for ADDITION4)
real(kind=8),allocatable::mass_ed(:)  ! masses of the educts (up to four for ADDITION4)
real(kind=8)::beta  ! inverse temperature
real(kind=8)::andersen_time ! application interval for Andersen thermostat
character(len=40)::thermo ! which thermostat shall be used (character)?
!integer::thermostat   ! switch for thermostat: 0: Anderson, 1: GLE
integer::andersen_step  ! after how many steps the Andersen thermostat shall be applied
integer::recr_equi,child_tot,child_interv,child_point,child_evol ! recrossing calc. parameters 
integer::child_times
logical::recross_calc  ! signalize of a recrossing calculation is done currently
real(kind=8),allocatable::gle_mat_A(:,:),gle_mat_C(:,:)  ! the matrices used for the GLE thermostat
integer::gle_ns  ! number of dimensions for GLE thermostat
real(kind=8),allocatable::gle_S(:,:),gle_T(:,:),gle_p(:,:,:,:)  ! GLE thermostat arrays
real(kind=8),allocatable::gle_np(:,:,:,:) !  GLE thermostat arrays 
logical::gen_test  ! perform only initial start structure generations
integer::gen_pr_frac,gen_pr_act   ! only each N'th structure, energy etc shall be written out
real(kind=8)::gen_energies  ! total sum of energies during structure generation
! For box simulations without periodic boundary conditions 
logical::box_nonp
real(kind=8)::box_center(3)
real(kind=8)::box_size
integer::box_fixed(1000)
integer::num_fixed
! For constant pressure npt simulations with a barostat
logical::npt
real(kind=8)::pressure

!   If the starting velocity shall be read in
logical::read_vel
real(kind=8),allocatable::vel_start(:,:)

!   For QMDFF parameter optimization
real(kind=8),allocatable::xyz_ref_lm(:,:,:)

!   the averaged kinetic energy 
real(kind=8)::ekin_avg 
real(kind=8)::ekin2_avg
logical::calc_ekin
integer::ekin_num
integer::ekin_atoms(100000)

!    change mass of one element/isotope
logical::change_mass
real(kind=8)::newmass
character(len=2)::elem_mass
! for orca calculations (ab initio MD)
logical::orca
character(len=100)::orca_com
! for the water SPC model
logical::water_spc
real(kind=8),allocatable::water_pars(:) ! all needed parameters
integer::nwater  ! number of water molecules 
integer,allocatable::water_act(:)  ! to which water molecule the atoms belong
!  for boxes with hard walls
logical::box_walls
real(kind=8)::walldim_x,walldim_y,walldim_z
!   For the Nose-Hoover/Berendsen barostat
integer::barostat   ! switch: 1: Nose-Hoover, 2: Berendsen
real(kind=8)::vel_baro
real(kind=8)::temp_test  ! for andersen testing..
real(kind=8)::press_act   ! the actual pressure to be calculated
real(kind=8)::mass_tot   ! for density calculations (mass of the system)
! for acceleration of calculations: store previous infos
real(kind=8),allocatable::inter_old(:,:),change_tol
integer,allocatable::i_best_old(:)
integer::irc_local  ! readin: how many structures left/right shall be investigated
! for geometry optimization (minima, TS)
logical::newton_raphson  ! if the simpler Newton-Raphson method shall be used instead of P-RFO
integer::maxiter  ! maximal number of geoopt steps
real(kind=8)::stepmax ! maximal step size
real(kind=8)::ethr,gthr,dthr,gmaxthr,dmaxthr  ! convergence criteria (geoopt)
real(kind=8)::irc_ethr,irc_gthr  ! convergence criteria (IRC)
integer::rp_evb_mode  ! if only energies and gradients shall be used for RP-EVB coupling
! for error handling
logical::act_check    ! if error handling is activated at all..
logical::recross_check  ! if checking is activated for recrossing
logical::loose_check   ! samplings are not canceled after single error
integer::err_max,err_act_max
real(kind=8)::energy_tol  ! tolerance for energy in relation to TS energy per bead
!    TEST TEST TEST
real(kind=8),allocatable::error_vec_ens(:)  ! vector with energy parts 
real(kind=8),allocatable::error_vec_v12(:)  ! vector with coupling parts
real(kind=8)::xi_test  ! for manual shutdown
real(kind=8)::no_evb_xi  ! maximum xi value for pure QMDFF for no-evb
! for time measuring 
real(kind=8)::time1,time2
real(kind=8)::time(10),duration
end module evb_mod
