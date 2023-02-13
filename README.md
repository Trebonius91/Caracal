# EVB-QMDFF
RPMD and rate constant calculations on black-box potential energy surfaces

## About EVB-QMDFF

EVB-QMDFF is a free open-source software package that enables a wide variety of molecular dynamics applications.
Quantum mechanical derived force fields (QMDFFs) can be generated in black box fashion for arbitrary chemical systems, two QMDFFs can be coupled by 
different EVB (Empirical Valence Bond) coupling methods as well as reaction-path based methods such as TREQ.
Based on the EVB-QMDFF potential energy surfaces, Ring Polymer Molecular Dynamics (RPMD) can be run for several MD applications such as molecular force
simulations or high-quality chemical reaction rate constant calculations.

## License

EVB-QMDFF is distributed unter the terms of the [MIT license](https://opensource.org/licenses/mit-license):

```
   EVB-QMDFF - RPMD molecular dynamics and rate constant calculations on
               black-box generated potential energy surfaces

   Copyright (c) 2021 by Julien Steffen (steffen@pctc.uni-kiel.de)
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

## Obtaining EVB-QMDFF

The program package can be obtained by getting a copy of this repository using git:
```
git clone git://github.com/Trebonius91/EVB-QMDFF.git
```
## Installation:

EVB-QMDFF is written in Fortran 90, with some minor parts in Fortran 77.
For compilation, the following dependencies are required:

- Standard Fortran 90 compiler (gfortran, ifort, etc.)
- Lapack and BLAS  
- FFTW 
- MPI

## Compiling from Source

The Makefile, which is located in the main directory, should be modified to meet your sytem requirements 
(a separate configure file will be added in the future). After, this, copy ot tp the src directory and run
```
$ make 
```
If the build is successful, the different programs of the package are located in the bin directory. 

## Executing the programs

In order to understand the handling of the different programs in EVB-QMDFF, type the help option, e.g. 
```
$ dynamic.x -h
```
Now, general instructions how to use the program as well as a list of available keywords is shown.

## Tutorials and Examples

A detailed set of tutorials describing the different programs, including a set of different examples covering relevant application cases can be found in the [EVB-QMDFF Wiki](https://github.com/Trebonius91/EVB-QMDFF/wiki)

## Future improvements

- A configure file for automatic detection of system requirements will be written.
- A manual (PDF and HTML) with detailed handling and background information will be made available.
- A number of test jobs for ensuring correct functionality of the program will be provided.
- A number of example calculations serving as templates for your applications on the several topics that is EVB-QMDFF able to handle will be added. 

