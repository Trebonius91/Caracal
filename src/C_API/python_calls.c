//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//   CARACAL - Ring polymer molecular dynamics and rate constant calculations
//             on black-box generated potential energy surfaces
//
//   Copyright (c) 2023 by Julien Steffen (mail@j-steffen.org)
//                         Stefan Grimme (grimme@thch.uni-bonn.de) (QMDFF code)
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the "Software"),
//   to deal in the Software without restriction, including without limitation
//   the rights to use, copy, modify, merge, publish, distribute, sublicense,
//   and/or sell copies of the Software, and to permit persons to whom the
//   Software is furnished to do so, subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be included in
//   all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//   DEALINGS IN THE SOFTWARE.
//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//
//     C-file python_calls: Wrapper routines needed for a direct 
//     communication between Fortran and Python for communication
//     with ASE
//

#include <Python.h>
#include <stdio.h>
#include <stdbool.h>

static PyObject *pModule = NULL;

//     
//     The initialization routine: Read in the force field and define the 
//     structure of the system.
//
void init_mace(char *mlip_file, char *coord_file, bool *set_disp) {
   PyObject *pName, *pFunc, *pArgs;

//
//     Initiaize Python interpreter
//
   Py_Initialize();
   if (!Py_IsInitialized()) {
      fprintf(stderr, "Python interpreter not found!");
      return;
   }
//
//     Source directory of Python file added to sys.path
//
   PyRun_SimpleString("import sys\nsys.path.append(\"" PYTHON_SRC_DIR "\")\n");

//
//     Now call the actual Python code, first define the function!
//
   pName = PyUnicode_FromString("call_mace");
   pModule = PyImport_Import(pName);
   Py_DECREF(pName);

    if (pModule == NULL) {
        pModule = PyImport_ImportModule("call_mace");
        if (!pModule) {
            PyErr_Print();
            return;
        }
    }


   if (pModule != NULL) {
//
//     Get function from module
//
      pFunc = PyObject_GetAttrString(pModule, "init_mace");
      if (pFunc && PyCallable_Check(pFunc)) {

//
//     The transferred object with all input information
//         
         PyObject *pArgs = PyTuple_New(3);
//
//     Pack MLIP filename into Python string
//
         PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(mlip_file));
//
//     Pack structure filename into Python string
//
         PyTuple_SetItem(pArgs, 1, PyUnicode_FromString(coord_file));
//
//     Pack dispersion switch into Python string
//
         PyTuple_SetItem(pArgs, 2, PyBool_FromLong(*set_disp ? 1 : 0));
//
//     Call Python function
//
         PyObject *pResult = PyObject_CallObject(pFunc, pArgs);
         Py_DECREF(pArgs);
         if (pResult == NULL) {
            PyErr_Print();
            fprintf(stderr, "Error in MACE initialization python routine! \n");
            exit(EXIT_FAILURE);
         } else {
            Py_DECREF(pResult);
         }
       //  PyObject_CallObject(pFunc, pArgs);
       //  Py_DECREF(pArgs);
      } else {
         PyErr_Print();
      }
//      Py_DECREF(pFunc);
//      Py_DECREF(pModule);
   } else {
      PyErr_Print();
   }
}


//
//     The energy and gradient calculation routine: give the structure and the 
//     unit cell and the number of atoms, return the energy and gradient
//
void ase_mace(double *coords, double *unitcell, double *energy, double *gradient, int *natoms) {
   PyObject *pName, *pFunc;
   PyObject *pArgs, *pCoordList, *pUnitcellList;
   PyObject *pResult, *pEnergy, *pGradList;
   int i, j, n = *natoms;

//
//     Now call the actual Python code, first define the function!
//   
   if (pModule != NULL) {
//
//     Obtain function definition from Python
//
      pFunc = PyObject_GetAttrString(pModule, "ase_mace");
      if (pFunc && PyCallable_Check(pFunc)) {
//
//     Build list of coordinates for input
// 
         pCoordList = PyList_New(n);
         for (i=0; i<n; i++) {
            PyObject *row = PyList_New(3); 
 
            for (j = 0; j < 3; j++) {
               PyList_SetItem(row, j, PyFloat_FromDouble(coords[3*i + j]));          
            } 
            PyList_SetItem(pCoordList,i,row);
         }
//
//     Build the list with the unit cell shape (3,3)
//
         pUnitcellList = PyList_New(3);
         for (i = 0; i < 3; i++) {
             PyObject *row = PyList_New(3);
             for (j = 0; j < 3; j++) {
                 PyList_SetItem(row, j, PyFloat_FromDouble(unitcell[3*i + j]));
             }
             PyList_SetItem(pUnitcellList, i, row);
         }
//
//     Put together the arguments for the Python routine call
//
         pArgs = PyTuple_New(3);
         PyTuple_SetItem(pArgs, 0, pCoordList);
         PyTuple_SetItem(pArgs, 1, pUnitcellList);
         PyTuple_SetItem(pArgs, 2, PyLong_FromLong(n));
//
//     Call the Python energy+gradient calculation ASE function
//
         pResult = PyObject_CallObject(pFunc, pArgs);

//
//     Check if the results tuple is reasonable
//
         if (pResult == NULL) {
            PyErr_Print();
            fprintf(stderr, "Error in MACE energy+gradient routine!\n");
            exit(EXIT_FAILURE);
            return;
         }

         if (!PyTuple_Check(pResult) || PyTuple_Size(pResult) != 2) {
            fprintf(stderr, "ase_mace did not return a 2-item tuple.\n");
            Py_DECREF(pResult);
            return;
         }

         Py_DECREF(pArgs);
         if (pResult != NULL) {
//
//     Extract the results: expecting (energy, gradient)
//
            pEnergy = PyTuple_GetItem(pResult, 0);
            *energy = PyFloat_AsDouble(pEnergy);

            pGradList = PyTuple_GetItem(pResult, 1);
            for (i = 0; i < n; i++) {
               PyObject *row = PyList_GetItem(pGradList, i);
               for (j = 0; j < 3; j++) {
                  gradient[3*i + j] = PyFloat_AsDouble(PyList_GetItem(row, j));
               }
            }
            Py_DECREF(pResult);
         } else {
            PyErr_Print();
         }
         Py_DECREF(pFunc);
      } else {
         PyErr_Print();
      }

   } else {
      fprintf(stderr, "Python interpreter not found!\n");
      PyErr_Print();
   }
}
