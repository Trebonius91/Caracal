# Makefile for compiling Caracal
# Edited by Julien Steffen, 2023

#   Fortran compiler for serial version:
#   GNU Fortran compiler
#FC = gfortran  
#   Intel Fortran compiler
#FC = ifort  
#   GNU/Intel Fortran compiler for MPI version
FC = mpifort
#   GNU/Intel C compiler for the ASE wrapper
CC = mpicc

#   Location of the compiled and linked executables
BINDIR = ~/bin

#   The Caracal library file 
LIBRARY = libcaracal.a
CURRENT = $(shell pwd)

#   The Python.h header files for the ASE wrapper
PYINCLUDES = $(shell python3-config --includes)
#   The absolute path of the directory with the Python files
PYTHON_SRC_DIR = $(shell realpath ./C_API)
#   For the final linking: the explicit location of the Python library
PYTHON_LINK = $(shell python3-config --ldflags)
#   The name of the Python library, depending on the current version
PYTHON_LIB = -lpython$(shell python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")


#   Location of the GULP package, if compiled with Caracal
#   Choose either Linux (for Linux) or Darwin (for Mac) 
#OSTYPE = Linux
#   The absolute or relative path to the GULP main directory
#GULPDIR =  ../../gulp-6.2
#   Comment in for the compiler flags 
#GULPFF = -DGULP=def -I${GULPDIR}/Src/${OSTYPE}
#   Comment in for the linking flags
#GULPLINK = ${GULPDIR}/Src/${OSTYPE}/libgulp.a  ${GULPDIR}/Utils/pGFNFF/Src/libpGFNFF.a 

#   The source files directory
SRCDIR = src
#   Fortran Compiler flags, should usually work
FFLAGS = -cpp -fopenmp ${GULPFF} -fno-align-commons -DPARALLEL -fallow-argument-mismatch -O1 
#FFLAGS = -cpp -fopenmp ${GULPFF} -fno-align-commons -fallow-argument-mismatch -O1
#   Fortran Compiler flags, debug version
#FFLAGS =  -fno-align-commons -g -ffpe-trap=zero,invalid,overflow,underflow  #-ffree-form #-Wall # debug version!
#   C Compiler flags for ASE wrapper
CFLAGS = ${PYINCLUDES} -O2 -DPYTHON_SRC_DIR=\"$(PYTHON_SRC_DIR)\"

#   Link against ownstanding Lapack, BLAS and FFTW libraries (probably deprecated)
#LINKFLAGS = -static-libgcc -fopenmp -llapack -lblas -lfftw3 -fno-align-commons # normal version 
#   link against Intel MKL, has been tested with GNU Fortran MPI compiler
LINKFLAGS = -fopenmp -static-libgcc ${PYINCLUDES} \
     -DPYTHON_SRC_DIR=\"$(PYTHON_SRC_DIR)\" -L${MKLROOT}/lib/intel64 \
     -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
     -lgomp -lpthread ${PYTHON_LINK} ${PYTHON_LIB} -lm -ldl 

LINKLIBS = -lgfortran

#   Targets by default
.PHONY: all
.PHONY: clean

#   All objects, listed:
#    - modules
#    - QMDFF subroutines 
#    - EVB subroutines
#    - MCTC subroutines
#    - DFTD3 subroutines
#    - Multicharge subroutines
#    - DFTD4 subroutines
#    - TOMLF subroutines
#    - tblite subroutines
#    - Other analytic PES
#    - program files
OBJS = general.o evb_mod.o qmdff.o lm_module.o debug.o h2co_mod.o \
       pbc_mod.o fftw_mod.o xtb_mod.o inter_mace.o \
       \
       moment.o ff_e.o eabh0.o thermo.o axis.o thermocal.o heat.o \
       bonde.o gethirsh.o docm5.o cm5mod.o copyc6.o limit.o rsp.o \
       basis0_eht.o basis_eht.o setsto3.o setsto4.o dex2.o \
       ehtfull.o eht.o \
       stints.o fermismear_eht.o gab.o compmodel.o ff_anal.o pqn.o   \
       ff_bond.o pairs.o lina.o ff_clean.o ff_mod.o rdfix.o getf.o \
       hbpara.o checktrip.o doeht.o abdamp.o eabhag.o eabxag.o   \
       eabx.o egrestrain.o ff_nonb.o ff_hb.o ff_eg.o ff_nonbe.o \
       ff_egh.o ff_hb_two.o ff_eg_two.o ff_nonb_two.o ff_hb_three.o \
       ff_eg_three.o ff_nonb_three.o bangle.o valijkl.o valijk.o \
       crossprod.o vecnorm.o omega.o impsc.o dphidr.o domegadr.o \
       ff_neighbor.o rsort.o isort.o nneighbor.o outgeo.o ff_set.o \
       modelgeo.o covbond.o atomthere.o avdamp.o getrot.o \
       checkfour.o setnonb.o valel.o setr0.o ncoord_qmdff.o getc6.o \
       gaurd0.o gaurd.o rdghess.o hess.o g98fake.o trproj.o \
       gtrprojm.o blckmgs.o dsyprj.o htosq.o htpack1.o htpack2.o \
       preig2.o preig4.o rdhess.o rdohess.o rdchess.o hfit.o \
       getpar0.o getparini.o putpar.o pderiv.o getpar.o ffhess.o \
       lmnpre.o pola.o rhftce.o prod.o main_gen.o procload.o \
       pfit.o getxy.o rd0.o rdo0.o rdc0.o rd.o rdo.o rdc.o rdv.o \
       rdv0.o rdvhess.o atommass_loc.o \
       elem.o asym.o upper.o readl.o readaa.o cutline.o checktype.o \
       rdsolvff0.o rdsolvff.o rdsolvff_two.o rdsolvff_three.o \
       rdwbo.o getring.o samering.o chk.o setvalel.o setzetaandip.o \
       setmetal.o splitmol.o split2.o rotfrag.o calcrotation.o \
       fragmentation.o vadd.o vsub.o vscal.o vlen.o crprod.o warn.o \
       warn2.o lm_good.o rdgwbo.o qmdff_corr.o \
       \
       atommass.o build_dmat.o getkey.o freeunit.o gettext.o \
       upcase.o nextarg.o getword.o basefile.o trimtext.o suffix.o \
       version0.o lowcase.o optimize_de.o optimize_dq.o lm_de_func.o \
       lm_dq_func.o deltaq.o dg_evb_init_ser.o solving_lgs.o \
       lm_function.o mdstat.o mdrest.o invert.o verlet.o maxwell.o \
       ranvec.o erfinv.o erf.o erfcore.o temper.o matinv3.o \
       normal.o calendar.o egqmdff.o eqmdff3.o eqmdff.o \
       fatal.o getxyz.o gmres.o gradient.o gradnum.o hessevb.o \
       hessqmdff.o prog_initial.o init_int.o xyz_2int.o grad2int.o \
       int2grad.o hess2int.o diagonalize.o dist.o ang.o oop.o \
       dihed.o optimize_3evb.o energy_1g.o mat_diag.o dmatinv.o \
       mdinit.o next_geo.o prepare.o promo.o promo_log.o random.o \
       read2int.o read_geo.o read_grad.o read_hess.o read_pes.o \
       geoopt.o calc_freq.o sum_v12.o umbr_coord.o \
       umbrella.o wham.o prob.o supdot.o calc_xi.o calc_com.o\
       andersen.o recross.o constrain_q.o constrain_p.o\
       calc_k_t.o bonds_ref.o pre_sample.o get_centroid.o rfft.o  \
       irfft.o rp_evb_init.o build_spline.o interp_spline.o \
       spline_dist.o project_hess.o project_hess_xyz.o egqmdff_corr.o\
       calc_wilson.o calc_dwilson.o opt_ts.o geoopt_int.o pseudoinv.o  \
       int2step.o sign_func.o krondel.o opt_irc.o sum_dv12.o \
       rpmd_check.o recross_serial.o dg_evb_init_par.o help.o\
       int_diff.o fact.o grad_test.o atomname.o lin_reg.o \
       matrix_exp.o cholesky.o cholesky.o calc_rate_read.o \
       calc_rate_info.o egrad_dg_evb.o egrad_treq.o calc_d2vec.o \
       transrot.o nhc.o orca_grad.o set_periodic.o box_image.o \
       ewald_recip.o bsplgen.o setchunk.o ewald_adjust.o ewald_fft.o \
       bspline.o dftmod.o exp_switch.o nhc_npt.o opt_qmdff_ser.o \
       lm_qmdff.o water_init.o egrad_water.o ev_coord_init.o \
       external_grad.o custom_init.o custom_grad.o pgfn_init.o \
       egrad_pgfn.o egrad_mace.o mace_init.o  \
       \
       mctc/env/accuracy.o mctc/env/error.o mctc/env/system.o \
       mctc/env/testing.o mctc/env.o mctc/io/codata2018.o mctc/io/constants.o \
       mctc/io/convert.o mctc/io/filetype.o mctc/io/math.o \
       mctc/io/resize.o mctc/io/symbols.o mctc/io/structure/info.o \
       mctc/io/structure.o mctc/io/utils.o \
       mctc/io/read/aims.o  mctc/io/read/ctfile.o \
       mctc/io/read/gaussian.o mctc/io/read/genformat.o mctc/io/read/pdb.o \
       mctc/io/read/qchem.o mctc/io/read/turbomole.o \
       mctc/io/read/vasp.o mctc/io/read/xyz.o mctc/io/read.o \
       mctc/io/resize.o mctc/io/structure.o mctc/io/symbols.o \
       mctc/io/utils.o mctc/io/write/aims.o mctc/io/write/ctfile.o \
       mctc/io/write/gaussian.o mctc/io/write/genformat.o \
       mctc/io/write/pdb.o mctc/io/write/qchem.o mctc/io/write/turbomole.o \
       mctc/io/write/vasp.o mctc/io/write/xyz.o \
       mctc/io/write.o mctc/env.o mctc/io.o mctc/version.o \
       \
       dftd3/dftd3/cutoff.o dftd3/dftd3/blas.o dftd3/dftd3/damping.o \
       dftd3/dftd3/data/covrad.o dftd3/dftd3/data/r4r2.o \
       dftd3/dftd3/data/vdwrad.o  dftd3/dftd3/data.o \
       dftd3/dftd3/reference.o dftd3/dftd3/model.o dftd3/dftd3/ncoord.o \
       dftd3/dftd3/disp.o dftd3/dftd3/damping/atm.o dftd3/dftd3/param.o \
       dftd3/dftd3/damping/mzero.o dftd3/dftd3/damping/optimizedpower.o \
       dftd3/dftd3/damping/rational.o dftd3/dftd3/damping/zero.o  \
       dftd3/dftd3/version.o \
       dftd3/dftd3.o  \
       \
       multicharge/cutoff.o multicharge/data/covrad.o multicharge/data.o \
       multicharge/blas.o multicharge/ewald.o multicharge/lapack.o \
       multicharge/wignerseitz.o multicharge/model.o multicharge/ncoord.o \
       multicharge/output.o multicharge/param/eeq2019.o multicharge/param.o \
       multicharge/version.o multicharge.o \
       \
       dftd4/dftd4/cutoff.o dftd4/dftd4/blas.o dftd4/dftd4/charge.o \
       dftd4/dftd4/damping.o dftd4/dftd4/data/covrad.o dftd4/dftd4/data/en.o \
       dftd4/dftd4/data/hardness.o dftd4/dftd4/data/r4r2.o \
       dftd4/dftd4/data/zeff.o  dftd4/dftd4/data.o dftd4/dftd4/reference.o \
       dftd4/dftd4/model.o dftd4/dftd4/ncoord.o \
       dftd4/dftd4/disp.o dftd4/dftd4/numdiff.o dftd4/dftd4/damping/atm.o \
       dftd4/dftd4/damping/rational.o dftd4/dftd4/param.o \
       dftd4/dftd4/version.o  dftd4/dftd4.o \
       \
       tomlf/tomlf/constants.o tomlf/tomlf/datetime.o tomlf/tomlf/error.o \
       tomlf/tomlf/utils/io.o tomlf/tomlf/utils.o tomlf/tomlf/type/value.o \
       tomlf/tomlf/structure/list.o tomlf/tomlf/structure/map.o \
       tomlf/tomlf/structure/node.o tomlf/tomlf/structure/array_list.o \
       tomlf/tomlf/structure/ordered_map.o tomlf/tomlf/structure.o \
       tomlf/tomlf/type/array.o tomlf/tomlf/type/keyval.o \
       tomlf/tomlf/type/table.o tomlf/tomlf/type.o \
       tomlf/tomlf/build/keyval.o \
       tomlf/tomlf/build/array.o tomlf/tomlf/build/merge.o \
       tomlf/tomlf/build/table.o tomlf/tomlf/build/path.o  tomlf/tomlf/build.o \
       tomlf/tomlf/de/token.o tomlf/tomlf/terminal.o tomlf/tomlf/diagnostic.o \
       tomlf/tomlf/de/context.o tomlf/tomlf/de/abc.o tomlf/tomlf/de/lexer.o \
       tomlf/tomlf/de/parser.o tomlf/tomlf/de.o tomlf/tomlf/ser.o \
       tomlf/tomlf/utils/sort.o tomlf/tomlf/version.o \
       tomlf/tomlf.o \
       \
       tblite/adjlist.o tblite/basis/type.o tblite/basis/ortho.o \
       tblite/basis/slater.o tblite/basis.o tblite/blas/level1.o \
       tblite/blas/level2.o tblite/blas/level3.o tblite/blas.o \
       tblite/container/cache.o tblite/scf/info.o tblite/integral/type.o \
       tblite/wavefunction/spin.o tblite/scf/potential.o \
       tblite/wavefunction/type.o tblite/container/type.o \
       tblite/output/format.o \
       tblite/container/list.o tblite/container.o \
       tblite/disp/cache.o tblite/cutoff.o tblite/disp/type.o \
       tblite/disp/d3.o tblite/disp/d4.o tblite/disp.o \
       tblite/wavefunction/fermi.o tblite/wavefunction/mulliken.o \
       tblite/coulomb/ewald.o tblite/wignerseitz.o tblite/coulomb/cache.o \
       tblite/coulomb/type.o tblite/coulomb/charge/type.o \
       tblite/coulomb/charge/effective.o tblite/coulomb/charge/gamma.o \
       tblite/coulomb/charge.o tblite/data/covrad.o tblite/ncoord/type.o \
       tblite/ncoord/gfn.o tblite/coulomb/multipole.o \
       tblite/coulomb/thirdorder.o  tblite/coulomb.o \
       tblite/xtb/coulomb.o tblite/lapack/getrf.o tblite/lapack/getri.o \
       tblite/lapack/getrs.o tblite/lapack.o tblite/scf/mixer/type.o \
       tblite/scf/mixer/broyden.o tblite/scf/mixer.o \
       tblite/scf/solver.o \
       tblite/scf/iterator.o  tblite/container/type.o \
       tblite/container/list.o tblite/container.o \
       tblite/ncoord/exp.o tblite/data/covrad_ceh.o \
       tblite/ncoord/ceh_std.o tblite/ncoord/ceh_en.o tblite/ncoord.o \
       tblite/toml.o tblite/param/serde.o \
       tblite/param/charge.o tblite/param/dispersion.o \
       tblite/data/paulingen.o \
       tblite/param/element.o tblite/param/halogen.o \
       tblite/param/hamiltonian.o tblite/param/multipole.o \
       tblite/param/repulsion.o tblite/param/thirdorder.o  \
       tblite/param/mask.o tblite/param.o tblite/ceh/h0.o \
       tblite/ceh/calculator.o tblite/context/logger.o \
       tblite/context/terminal.o tblite/context/solver.o \
       tblite/lapack/sygvd.o tblite/lapack/sygst.o tblite/lapack/potrf.o \
       tblite/lapack/sygvr.o tblite/lapack/solver.o \
       tblite/context/type.o tblite/context.o tblite/integral/trafo.o \
       tblite/integral/diat_trafo.o tblite/integral/dipole.o \
       tblite/integral/multipole.o tblite/data/atomicrad.o tblite/xtb/spec.o \
       tblite/xtb/h0.o tblite/repulsion/type.o tblite/classical/halogen.o \
       tblite/repulsion/effective.o tblite/repulsion.o tblite/xtb/calculator.o \
       tblite/wavefunction/guess.o tblite/wavefunction.o tblite/scf.o \
       tblite/external/field.o tblite/timer.o tblite/ceh/ceh.o \
       tblite/version.o tblite/output/ascii.o tblite/output/property.o \
       tblite/results.o tblite/xtb/singlepoint.o tblite/xtb/gfn1.o \
       tblite/xtb/gfn2.o tblite/xtb/ipea1.o \
       tblite/api/version.o tblite/api/utils.o tblite/api/error.o \
       tblite/api/table.o tblite/api/result.o tblite/api/structure.o \
       tblite/api/context.o tblite/api/param.o \
       tblite/xtb.o tblite/mesh/lebedev.o tblite/solvation/born.o \
       tblite/solvation/data.o tblite/solvation/type.o tblite/solvation/alpb.o \
       tblite/solvation/cpcm_dd.o tblite/solvation/cpcm.o  \
       tblite/solvation/input.o  tblite/solvation.o tblite/data/spin.o \
       tblite/spin.o tblite/api/calculator.o tblite/api/container.o \
       wibsort.o egrad_tblite.o \
       aenet/ext/arglib.o aenet/ext/chebyshev.o aenet/ext/feedforward.o \
       aenet/ext/io.o aenet/ext/sortlib.o aenet/ext/sfbasis.o \
       aenet/ext/lclist.o aenet/ext/symmfunc.o aenet/ext/timing.o \
       aenet/ext/unittest.o aenet/ext/xsflib.o aenet/constants.o \
       aenet/aeio.o aenet/geometry.o aenet/input.o aenet/sfsetup.o \
       aenet/trainset.o aenet/potential.o aenet/random.o \
       aenet/parallel.o aenet/optimize.o aenet/aenet.o \
       aenet/initialize.o aenet/energy.o egrad_aenet.o \
       \
       egrad_ch4oh.o util_ch4oh.o egrad_h3.o \
       numeral.o getnumb.o getstring.o torphase.o egrad_ch4h.o \
       util_ch4h.o egrad_brh2.o egrad_oh3.o util_oh3.o egrad_geh4oh.o\
       util_geh4oh.o egrad_c2h7.o egrad_clnh3.o util_clnh3.o \
       egrad_mueller.o egrad_o3.o egrad_ch4cn.o util_ch4cn.o  \
       egrad_nh3oh.o util_nh3oh.o main_h2co.o util_h2co.o \
       egrad_c2h6cn.o \
       \
       C_API/python_calls.o \
       \
       qmdffgen.o evbopt.o dynamic.o explore.o calc_rate.o \
       black_box.o poly_qmdff.o stick_coeff.o
 
EXES = qmdffgen.x dynamic.x evbopt.x explore.x calc_rate.x \
       black_box.x poly_qmdff.x stick_coeff.x

all: $(OBJS)  $(EXES)

#    compile all .f or .f90 files to .o files

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

#    Create the library-file
$(LIBRARY): $(OBJS)
	ar -crusv $(LIBRARY) $(OBJS)


#    Make finally the exetuables and copy them into the bin directory
#    Generate the directory on the fly
# $@ is the actual element on the list
%.x: %.o libcaracal.a ${GULPLINK} 
	${FC} ${LINKFLAGS} $^ $(LINKLIBS) -o $@ ; strip $@
	$(shell mkdir -p ../bin)
	cp $@ ../bin/
	mv $@ $(BINDIR) 


#    remove all object and executable data
clean:
	rm -f *.o $(PROG)
	rm -f */*.o $(PROG)
	rm -f */*/*.o $(PROG)
	rm -f */*/*/*.o $(PROG)
	rm -f *.x $(PROG) 
	rm -f *.mod $(PROG)
	rm -f libqmdff.a
	rm ../bin/*.x

