

<p align="center">
<img src="https://github.com/Trebonius91/Caracal/blob/main/manual/figures/logo.png" alt="drawing" width="300"/>
</p>

<p align="center">
<sub><sup>(Logo reprinted with permission from https://doi.org/10.1021/acs.jctc.3c00568. Copyright 2023 American Chemical Society.)</sup></sub>
</p>

# Caracal
Ring polymer molecular dynamics and rate constant calculations on black-box potential energy surfaces

## About Caracal

Caracal is a free open-source software package that enables a wide variety of molecular dynamics applications.
Unbiased and biased molecular dynamics trajectories can be sampled, unperiodic as well as periodic (NVE, NVT, NpT) setups are possible.

Classical (Newtonian) as well as ring polymer molecular dynamics (RPMD) can be used.
A central feature is the automized setup of potential energy surfaces for gas phase reactions using the EVB-QMDFF methodology.

Quantum mechanical derived force fields (QMDFFs) can be generated in black box fashion for arbitrary chemical systems, two QMDFFs can be coupled by
different EVB (Empirical Valence Bond) coupling methods.

Simple energy difference coupling methods as well as more sophisticated methods like distributed gaussian (DG)-EVB or transition region corrected reaction path EVB-QMDFF (TREQ) are availabe for that purpose.

Especially the TREQ method allows for black-box generation of high quality PES descriptions, the whole process of PES setup and rate constant calculation with RPMD is realized in the black-box program within Caracal.

The well-known class of semiempirical GFN-xTB methods is implemented now as well, enabling the black-box simulation of arbitrary chemical systems.

QMDFFs of single molecules can be polymerized to describe complex multicomponent mixtures or solutions, containing molecules like transition metal complexes for which it is hard to get force field parametrizations but are easy to setup with QMDFF.

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

## Compiling from Source

The Makefile, which is located in the main directory, should be modified to meet your system requirements, choosing a suitable set of compiler and required libraries
(a separate configure file will be added in the future). After, this, copy it to the src directory and run
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

## Future improvements

- Development of local versions of TREQ, enabling the calculation of larger systems without complicated and error-prone sets of internal coordinates
- Exact unimolecular rate constants: implementing theoretically stringent prefactors for the different unimolecular mechanisms
- Reactions in solvents: Implementing a "QM/MM" scheme, with the reactive part described by TREQ and sampled with RPMD and the solvent described by GFN-FF and classical dynamics

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
