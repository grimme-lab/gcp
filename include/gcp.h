/* This file is part of mctc-gcp.
* SPDX-Identifier: LGPL-3.0-or-later
*
* mctc-gcp is free software: you can redistribute it and/or modify it under
* the terms of the Lesser GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* mctc-gcp is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* Lesser GNU General Public License for more details.
*
* You should have received a copy of the Lesser GNU General Public License
* along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.
**/

#ifndef GCP_H
#define GCP_H

#ifdef __cplusplus
#define GCP_API_ENTRY extern "C"
#else
#define GCP_API_ENTRY extern
#endif

GCP_API_ENTRY void c_gcp_call(int* n,             // No. of atoms
                              double* xyz,        // coordinates (3,n) in Bohr
                              double* lat,        // lattice matrix
                              int* iz,            // element numbers
                              double* gcp_e,      // gcp energy
                              double* gcp_g,      // gcp gradient
                              double* gcp_glat,   // gcp lattice gradient
                              bool* dograd,       // flag: evaluate gradient
                              bool* dohess,       // flag: evaluate hessian (file: gcp_hessian)
                              bool* pbc,          // flag: periodic boundary conditions
                              const char* method, // method name
                              bool* echo,         // flag: verbose output
                              bool* parfile       // flag: print extended parameter file (gcp.param)
                              );

GCP_API_ENTRY void setr0ab(int*    max_elem,      // max element number
                           double* autoang,       // a.u. to Angstrom (0.52917726)
                           double* r              // vector of size max_elem x max_elem
                          );

#endif

