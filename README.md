

<p align="center">
<img src="https://github.com/Trebonius91/Caracal/blob/main/manual/figures/logo.png" alt="drawing" width="300"/>
</p>

<p align="center">
<sub><sup>(Logo reprinted with permission from https://doi.org/10.1021/acs.jctc.3c00568. Copyright 2023 American Chemical Society.)</sup></sub>
</p>

# Caracal
Ring polymer molecular dynamics and rate constant calculations on black-box and machine-learning potential energy surfaces

## About Caracal

Caracal is a free open-source software package that enables a wide variety of molecular dynamics applications.
Unbiased and biased molecular dynamics trajectories can be sampled, unperiodic as well as periodic (NVE, NVT, NpT) setups are possible.

Classical (Newtonian) as well as ring polymer molecular dynamics (RPMD) can be used.
A central feature is the automized setup of potential energy surfaces for gas phase reactions using the EVB-QMDFF methodology.
Further, numerous different machine-learning interatomic potentials (MLIPs) can be utilized for all tasks.

Quantum mechanical derived force fields (QMDFFs) can be generated in black box fashion for arbitrary chemical systems, two QMDFFs can be coupled by
different EVB (Empirical Valence Bond) coupling methods.

A growing number of MLIPs is included or plugged into Caracal: Behler-Parrinello neural networks via the Aenet program, message-passing atomic cluster expansion (MACE) and the universal model for atoms (UMA). Caracal supports the self-consistent fine tuning of MACE MLIPs for biased and unbiased samplings.

The well-known class of semiempirical GFN-xTB methods is implemented now as well, enabling the black-box simulation of arbitrary chemical systems.

Further, a steadily growing number of analytical PES representations of gas phase reaction systems from the literature is integrated into Caracal, mainly for the calculation of rate constants. Alternatively, you can couple your own
PES routine to the program (either by defining a syscall or directly integrating the subroutine into the Caracal code) and perform biased or unbiased RPMD with it.

If you use Caracal for your research, I would be very glad if you could cite the ***original Caracal paper*** ([DOI: 10.1021/acs.jctc.3c00568](https://pubs.acs.org/doi/10.1021/acs.jctc.3c00568)). 
Within it and its supporting information, more details concerning the underlying theory and algorithms can be found as well.

Included (and sometimes slightly modified) software:

mctc_lib: [https://github.com/grimme-lab/mctc-lib](https://github.com/grimme-lab/mctc-lib)

dftd3_lib: [https://github.com/dftbplus/dftd3-lib](https://github.com/dftbplus/dftd3-lib)

dftd4: [https://github.com/dftd4/dftd4](https://github.com/dftd4/dftd4)

multicharge: [https://github.com/grimme-lab/multicharge](https://github.com/grimme-lab/multicharge)

toml-f: [https://github.com/toml-f/toml-f](https://github.com/toml-f/toml-f)

tblite: [https://github.com/tblite/tblite](https://github.com/tblite/tblite)

aenet: [http://ann.atomistic.net/](http://ann.atomistic.net/)

## License

Caracal is distributed unter the terms of the [MIT license](https://opensource.org/licenses/mit-license):

```
   CARACAL - Ring polymer molecular dynamics and rate constant calculations
             on black-box generated potential energy surfaces

   Copyright (c) 2023 by Julien Steffen (mail@j-steffen.org)
                         Stefan Grimme (grimme@thch.uni-bonn.de) (QMDFF code)

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
   DEALINGS IN THE SOFTWARE.
```

## Obtaining Caracal

The program package can be obtained by getting a copy of this repository using git:
```
git clone git://github.com/Trebonius91/CARACAL.git
```
## Installation:

Caracal is written in Fortran 90, with some minor parts in Fortran 77.
For compilation, the following dependencies are required:

- Standard Fortran 90 compiler (gfortran, ifort, etc.)
- Lapack and BLAS
- FFTW
- MPI

If MACE and UMA shall be used, the respective Python libraries as well the 
atomic simulation environment (ASE) libraries need to be installed.
Missing li6raries have no impact on the general functionality of Caracal.

## Compiling from Source

The Makefile, which is located in the main directory, should be modified to meet your system requirements, choosing a suitable set of compiler and required libraries
(a separate configure file will be added in the future). 
Since Caracal now includes an API to call machine-learned interatomic potential (MLIP) energies+gradients from ASE, up to date Fortran and C compilers need to be present, as well as a recent Python version. 
If the latter differs too much from the first, broken dependencies can occur during linking.

After, this, copy it to the src directory and run
```
$ make
```
Currently, parallel make (-jN) is not recommended since problems with nested Fortran modules might arise.
If the build is successful, the different programs of the package are located in the bin directory.

If you want to use pGFN-FF as a PES description, the GULP program must be included as a patch. See [pGFN-FF with GULP](https://github.com/Trebonius91/CARACAL/wiki/dynamic#example-4-hh2-with-external-program) for how GULP can be obtained and the Caracal patch can be done.

## Executing the programs

In order to understand the handling of the different programs in Caracal, type the help option, e.g.
```
$ dynamic.x -h
```
Now, general instructions how to use the program as well as a list of available keywords is shown.

## Tutorials and Examples

A detailed set of tutorials describing the different programs, including a set of different examples covering relevant application cases can be found in the [Caracal Wiki](https://github.com/Trebonius91/CARACAL/wiki)

## Further Questions, Suggestions or Bug Reports

If you have further questions, suggestions what shall be added to Caracal, or detected a bug / unexpected behavior within Caracal, please do not hesitate to contact me! (either here on github, or via mail (mail@j-steffen.org)).        

## Included methods (citations)

Here, literature referring to the different methods included into Caracal are given. Please cite the one that you used during your research with Caracal.

- QMDFF: S. Grimme; J. Chem. Theory Comput. 2014, 10, 4497-4514 (DOI: [10.1021/ct500573f](https://doi.org/10.1021/ct500573f))
- dE-EVB: B. Hartke; S. Grimme; Phys. Chem. Chem. Phys. 2015, 17, 16715-16718 (DOI: [10.1039/C5CP02580J](https://doi.org/10.1039/C5CP02580J))
- DG-EVB/dQ-EVB: J. Steffen; B. Hartke; J. Chem. Phys. 2017, 147, 161701 (DOI: [10.1063/1.4979712](https://doi.org/10.1063/1.4979712))
- TREQ: J. Steffen; J. Chem. Phys. 2019, 150, 154105 (DOI: [10.1063/1.5092589](https://doi.org/10.1063/1.5092589))
- pGFN-FF (with GULP): J. Gale et al.; J. Chem. Theory Comput. 2021, 17, 12, 7827-7849 (DOI: [10.1021/acs.jctc.1c00832](https://doi.org/10.1021/acs.jctc.1c00832))
- GFN1-xTB: S. Grimme et al.; J. Chem. Theory Comput. 2017, 13, 1989-2009 (DOI: [10.1021/acs.jctc.7b00118](https://doi.org/10.1021/acs.jctc.7b00118)) 
- GFN2-xTB: C. Bannwarth et al.; J. Chem. Theory Comput. 2019, 15, 1652-1671 (DOI: [10.1021/acs.jctc.8b01176](https://doi.org/10.1021/acs.jctc.8b01176))
- ALPB solvation: S. Ehlert et al.; J. Chem. Theory Comput. 2021, 17, 4250-4261 (DOI: [10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471))
- Behler-Parrinello NNs via aenet: N. Artrith et al.: Comput. Mater. Sci. 2016, 114, 135-150 (DOI: [10.1016/j.commatsci.2015.11.047](https://doi.org/10.1016/j.commatsci.2015.11.047))
- MACE via ASE: I. Batatia et al.: arXiv:2206.07697, 2022 (DOI: [10.48550/arXiv.2206.07697](https://doi.org/10.48550/arXiv.2206.07697))
- UMA via ASE: B. M. Wood et al.: arXiv:2506.23971, 2025 (DOI: [10.48550/arXiv.2506.23971](https://doi.org/10.48550/arXiv.2506.23971)) 

## Future improvements

- Inclusion of additional MLIPs for their easy applications in unbiased and biased RPMD simulations
- Development of more targeted fine-tuning methodologies for an effective black-box fine-tuning of MLIPs 
- A full rate theory of RPMD rate constant calculations, for arbitrary unimolecular and bimolcular reactions
- Development of multiscale approaches including MLIPs, such as ML/MM approaches, with MLIP as inner and, e.g., GFN-FF as outer method.

## Latest changes

- 10/09/2023: The FFT routines were updated to the more recent Fortran2003/C interface as recommended by the FFTW developers in order to fix a bug when using the Intel MKL libraries
- 10/09/2023: The keyword RANDOM\_SEED [value] was added for the programs dynamic and calc\_rate in order to enable exact reproduction of calculations for debugging or benchmarks
- 10/30/2023: The program file 'poly_qmdff.f90' which was missing so far, was now added.
- 11/03/2023: Release version 1.1, addition of GFN-xTB methods (and related keywords)
- 12/18/2023: The input format for QMDFF force field has been slightly changed. A new QMDFF input section was introduced (see Caracal wiki and examples)
- 12/20/2023: An optional call to the GULP program as been added, pGFN-FF calculations can now be done!
- 01/15/2024: Internal coordinate analyses and Wilson matrix calculations can be done with explore.x, using the 'pes topol' option
- 02/26/2024: The position along the reaction path where the recrossing is calculated in calc_rate.x can now be chosen manually.
- 04/29/2024: The explore program has been improved. The conjugate gradient algirithm for geometry optimization added, better output for geometry optimization. IR intensities for pGFN-FF, GFN-xTB and QMDFF in frequency calculation.
- 06/21/2024: Neural network potentials as implemented in the aenet program are now available as an integrated PES function. Further, periodic structures can be given efficiently by the VASP POSCAR geometry format.
- 06/24/2024: Umbrella samplings and recrossing calculations were optimized, both can now be restarted more efficiently, the parallelized recrossing is now more straightforward.
- 07/20/2024: The dynamic.x program can now be parallelized with MPI to speed up the potential energy/gradient calculations. So far, only ANN potentials can be calculated with MPI.
- 07/20/2024: Recrossing calculations in calc_rate.x are now calculated with MPI automatically if more than one core is used, the keyword MPI in the RECROSS section is not needed anymore.
- 08/10/2024: Behler-Parrinello neural networks (2nd generation) are now included via the aenet program source code (accessible via the aenet_ann command).
- 06/04/2025: A new program stick_coeff.x enabling the automated calculation of molecular sticking coefficients on surfaces has been added.
- 07/09/2025: The inclusion of several modern machine-learned interatomic potential (MLIP) distributions like MACE has been made possible by a C-wrapper to the ASE program.
- 02/10/2026: The targeted printout of training set information for MLIP traininings or fine-tunings during calc_rate.x calculations as been added.
