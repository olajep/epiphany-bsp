/*
File: host_spmd.h

This file is part of the Epiphany BSP library.

Copyright (C) 2014 Buurlage Wits
Support e-mail: <info@buurlagewits.nl>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
and the GNU Lesser General Public License along with this program,
see the files COPYING and COPYING.LESSER. If not, see
<http://www.gnu.org/licenses/>.
*/

#include <host_bsp.h>
#include <host_bsp_inspector.h>

#include <stdlib.h>
#include <stdio.h>

#include "common.h"

#define DEBUG

// information on matrix and procs
int N = -1;
int M = -1;

// FIXME: Taken from our Zee library
// MATRIX TYPES
typedef struct
{
    int i;
    int j;
    float val;
} z_sp_entry_f;

// IMPORTANT: The task assumes the nonzeros be stored in column-major
// order. This makes it much easier for the fan-in and fan-out stages.
typedef struct
{
    int n;
    int m;
    int nz;
    z_sp_entry_f* entries;
} z_sp_mat_f;

typedef struct
{
    int n;
    int m;
    float* entries;
} z_mat_f;

// MATRIX FUNCTIONS
z_sp_mat_f eye(int n)
{
    z_sp_mat_f A;

    A.n = n;
    A.m = n;
    A.nz = n;

    A.entries = malloc(sizeof(z_sp_entry_f) * A.nz);
    for (int k = 0; k < A.nz; ++k) {
        A.entries[k].i = k;
        A.entries[k].j = k;
        A.entries[k].val = k;
    }

    return A;
}

z_mat_f rand(int n, int m)
{
    z_mat_f A;

    A.n = n;
    A.m = m;

    A.entries = malloc(sizeof(float) * n * m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A.entries[k] = rand() / RAND_MAX;

    return A;
}


z_mat_f zeros(int n, int m)
{
    z_mat_f A;

    A.n = n;
    A.m = m;

    A.entries = malloc(sizeof(float) * n * m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            A.entries[k] = 0.0f;

    return A;
}

void z_mat_destroy(z_mat_f A)
{
    free(A.entries);
}


void z_sp_mat_destroy(z_sp_mat_f A)
{
    free(A.entries);
}

//----------------------------------------------------------------------
// SpMV on Host
//----------------------------------------------------------------------

int main(int argc, char **argv)
{
    srand(12345);

    // COMPUTE u = Av
    // where u, v are vectors of dimension n
    // and A is a sparse (n x m) matrix
    
    // dimensions
    int n = 50;
    int m = 50;

    // initialize matrices and vectors
    z_sp_mat_f Id = eye(n);
    z_mat_f v = rand(n, 1);

    // initializes components (indices) of u
    z_mat_f u = zeros(n, 1);
    for(int i = 0; i < n; ++i)
        u.entries[i] = i;

    // initialize the BSP system
    bsp_init("bin/e_spmd.srec", argc, argv);
    bsp_begin(bsp_nprocs());

    // distribute the matrix
    switch (bsp_nprocs()) {
        case 16:
            N = 4;
            M = 4;
            break;

        case 64:
            N = 8;
            M = 8;
            break;

        default:
            fprintf(stderr, "Unsupported processor count, please add values\
for N and M in the host program.");
            return -1;
    }

    int* offsets = malloc(sizeof(int) * bsp_nprocs());
    int* offsets_cyc = malloc(sizeof(int) * bsp_nprocs());

    for (int i = 0; i < bsp_nprocs(); ++i)
    {
        offsets[i] = 0; 
        offsets_cyc[i] = 0;
    }

    // write matrix triples to random processors
    // to ensure that thi SpMV methods works with any distribution
    for (int i = 0; i < Id.nz; ++i) {
        int s = rand() % bsp_nprocs();
        ebsp_write(s,
                &Id.entries[i].i,
                (void*)(LOC_MAT_I + offsets[s] * sizeof(int)),
                sizeof(int));
        ebsp_write(s,
                &Id.entries[i].j,
                (void*)(LOC_MAT_J + offsets[s] * sizeof(int)),
                sizeof(int));
        ebsp_write(s,
                &Id.entries[i].val,
                (void*)(LOC_MAT_ENTRIES + offsets[s] * sizeof(float)),
                sizeof(float));
        offsets[s]++;
    }

    // write offsets as counts
    // .. FIXME
 
    for (int i = 0; i < bsp_nprocs(); ++i)
        offsets[i] = 0; 

    // write vector components to random processors
    // and write owner cyclically
    for (int i = 0; i < n; ++i) {
        int s = rand() % bsp_nprocs();
        ebsp_write(s,
                &i,
                (void*)(LOC_V_IDXS + offsets[s] * sizeof(int)),
                sizeof(int));
        ebsp_write(s,
                &v.entries[i],
                (void*)(LOC_V_ENTRIES + offsets[s] * sizeof(float)),
                sizeof(float));

        int t = i % bsp_nprocs;
        ebsp_write(t,
                &s,
                (void*)(LOC_V_OWNERS + offsets_cyc[t] * sizeof(int)),
                sizeof(int));

        offsets[s]++;
        offsets_cyc[t]++;
    }

    for (int i = 0; i < bsp_nprocs(); ++i)
    {
        offsets[s]++;
        offsets_cyc[i] = 0;
    }

    // write u owner to processors
    for (int i = 0; i < n; ++i) {
        int s = rand() % bsp_nprocs();
        int t = i % bsp_nprocs;
        ebsp_write(s,
                &i,
                (void*)(LOC_U_IDXS + offsets[t] * sizeof(int)),
                sizeof(int));
        ebsp_write(t,
                &s,
                (void*)(LOC_U_OWNERS + offsets_cyc[t] * sizeof(int)),
                sizeof(int));
        offsets[t]++;
        offsets_cyc[t]++;
    }

    free(offsets);
    free(offsets_cyc);

#ifdef DEBUG
    //ebsp_inspector_enable();
#endif

    ebsp_spmd();

    // read result
 
    bsp_end();

    return 0;
}
