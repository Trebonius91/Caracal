# Makefile for compiling the EVB-QMDFF program
# Edited by Julien Steffen, 2021

#  Fortran compiler for serial version:
#FC = gfortran  # GNU Fortran
# FC = ifort   # Intel Fortran compiler 
#  Fortran compiler for MPI version
 FC = mpif90  # GNU Fortran
# FC = mpifort # Intel Fortran compiler

#  Location of the compiled and linked executables
BINDIR = ~/bin

LIBRARY = libqmdff.a
CURRENT = $(shell pwd)

SRCDIR = src
FFLAGS =  -fno-align-commons -Wno-argument-mismatch -O1  #-ffree-form #-Wall # normal version
#FFLAGS =  -fno-align-commons -g -ffpe-trap=zero,invalid,overflow,underflow  #-ffree-form #-Wall # debug version!
LINKFLAGS = -static-libgcc -fopenmp -llapack -lblas -lfftw3 -fno-align-commons # normal version 
#LINKFLAGS = -static-libgcc -fopenmp -llapack -lblas -lfftw3 -fno-align-commons -g -ffpe-trap=zero,invalid,overflow,underflow # debug version!


# Targets by default
.PHONY: all
.PHONY: clean

# All objects, listed:
# - modules
# - QMDFF subroutines 
# - EVB subroutines
# - Other analytic PES
# - program files
OBJS = general.o evb_mod.o qmdff.o lm_module.o debug.o \
       \
       moment.o ff_e.o eabh0.o thermo.o axis.o thermocal.o heat.o \
       bonde.o gethirsh.o docm5.o cm5mod.o copyc6.o limit.o rsp.o \
       basis0.o basis.o setsto3.o setsto4.o dex2.o ehtfull.o eht.o \
       stints.o fermismear.o gab.o compmodel.o ff_anal.o pqn.o   \
       ff_bond.o pairs.o lina.o ff_clean.o ff_mod.o rdfix.o getf.o \
       hbpara.o checktrip.o doeht.o abdamp.o eabhag.o eabxag.o   \
       eabx.o egrestrain.o ff_nonb.o ff_hb.o ff_eg.o ff_nonbe.o \
       ff_egh.o ff_hb_two.o ff_eg_two.o ff_nonb_two.o ff_hb_three.o \
       ff_eg_three.o ff_nonb_three.o bangle.o valijkl.o valijk.o \
       crossprod.o vecnorm.o omega.o impsc.o dphidr.o domegadr.o \
       ff_neighbor.o sort.o isort.o nneighbor.o outgeo.o ff_set.o \
       modelgeo.o covbond.o atomthere.o avdamp.o getrot.o \
       checkfour.o setnonb.o valel.o setr0.o ncoord.o getc6.o \
       gaurd0.o gaurd.o rdghess.o hess.o g98fake.o trproj.o \
       gtrprojm.o blckmgs.o dsyprj.o htosq.o htpack1.o htpack2.o \
       preig2.o preig4.o rdhess.o rdohess.o rdchess.o hfit.o \
       getpar0.o getparini.o putpar.o pderiv.o getpar.o ffhess.o \
       lmnpre.o pola.o rhftce.o prod.o main_gen.o procload.o \
       pfit.o getxy.o rd0.o rdo0.o rdc0.o rd.o rdo.o rdc.o \
       elem.o asym.o upper.o readl.o readaa.o cutline.o checktype.o \
       rdsolvff0.o rdsolvff.o rdsolvff_two.o rdsolvff_three.o \
       rdwbo.o getring.o samering.o chk.o setvalel.o setzetaandip.o \
       setmetal.o splitmol.o split2.o rotfrag.o calcrotation.o \
       fragmentation.o vadd.o vsub.o vscal.o vlen.o crprod.o warn.o \
       warn2.o lm_good.o rdgwbo.o qmdff_corr.o \
       \
       atommass.o build_dmat.o getkey.o freeunit.o gettext.o \
       upcase.o nextarg.o getword.o basefile.o trimtext.o suffix.o \
       version.o lowcase.o optimize_de.o optimize_dq.o lm_de_func.o \
       lm_dq_func.o deltaq.o dg_evb_init_ser.o solving_lgs.o \
       lm_function.o mdstat.o mdrest.o invert.o verlet.o maxwell.o \
       ranvec.o erfinv.o erf.o erfcore.o temper.o  \
       normal.o calendar.o egqmdff.o eqmdff3.o eqmdff.o \
       fatal.o getxyz.o gmres.o gradient.o gradnum.o hessevb.o \
       hessqmdff.o initial.o init_int.o xyz_2int.o grad2int.o \
       int2grad.o hess2int.o diagonalize.o dist.o ang.o oop.o \
       dihed.o optimize_3evb.o energy_1g.o mat_diag.o dmatinv.o \
       mdinit.o next_geo.o prepare.o promo.o promo_log.o random.o \
       read2int.o read_geo.o read_grad.o read_hess.o read_evb.o \
       geoopt.o calc_freq.o sum_v12.o umbr_coord.o verlet_bias.o \
       umbrella.o wham.o prob.o supdot.o calc_xi.o calc_com.o\
       andersen.o mdinit_bias.o recross.o constrain_q.o constrain_p.o\
       calc_k_t.o bonds_ref.o pre_sample.o get_centroid.o rfft.o  \
       irfft.o rp_evb_init.o build_spline.o interp_spline.o \
       spline_dist.o project_hess.o project_hess_xyz.o egqmdff_corr.o\
       calc_wilson.o calc_dwilson.o opt_ts.o geoopt_int.o pseudoinv.o  \
       int2step.o sign_func.o krondel.o opt_irc.o sum_dv12.o \
       rpmd_check.o recross_serial.o dg_evb_init_par.o help.o\
       int_diff.o fact.o grad_test.o atomname.o lin_reg.o \
       matrix_exp.o cholesky.o cholesky.o rpmd_read.o \
       rpmd_info.o egrad_dg_evb.o egrad_treq.o calc_d2vec.o \
       transrot.o nhc.o orca_grad.o set_periodic.o box_image.o \
       ewald_recip.o bsplgen.o setchunk.o ewald_adjust.o ewald_fft.o \
       bspline.o dftmod.o exp_switch.o nhc_npt.o opt_qmdff_ser.o \
       lm_qmdff.o \
       \
       egrad_ch4oh.o util_ch4oh.o egrad_h3.o \
       numeral.o getnumb.o getstring.o torphase.o egrad_ch4h.o \
       util_ch4h.o egrad_brh2.o egrad_oh3.o util_oh3.o egrad_geh4oh.o\
       util_geh4oh.o egrad_c2h7.o egrad_clnh3.o util_clnh3.o \
       egrad_mueller.o egrad_o3.o egrad_ch4cn.o util_ch4cn.o  \
       egrad_nh3oh.o util_nh3oh.o \
       \
       qmdffgen.o evbopt.o egrad.o dynamic.o evb_qmdff.o rpmd.o \
       irc.o evb_kt_driver.o mult_qmdff.o qmdffopt.o
 
EXES = qmdffgen.x dynamic.x egrad.x evbopt.x evb_qmdff.x rpmd.x \
       irc.x evb_kt_driver.x mult_qmdff.x qmdffopt.x

all: $(OBJS)  $(EXES)

# compile all .f or .f90 files to .o files

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

#Create the library-file
$(LIBRARY): $(OBJS)
	ar -crusv $(LIBRARY) $(OBJS)


#Make finally the exetuables and copy them into the bin directory
#Generate the directory on the fly
# $@ is the actual element on the list
%.x: %.o libqmdff.a  
	${FC} ${LINKFLAGS} -o  $@ $^ $(LINKLIBS) ; strip $@
	$(shell mkdir -p ../bin)
	cp $@ ../bin/
	mv $@ $(BINDIR) 


#remove all object and executable data
clean:
	rm -f *.o $(PROG)
	rm -f *.x $(PROG) 
	rm -f *.mod $(PROG)
	rm -f libqmdff.a
	rm ../bin/*.x

