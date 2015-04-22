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
#include <assert.h>

#include "zee_matrix.h"
#include "common.h"

#define DEBUG

// information on matrix and procs
int N = -1;
int M = -1;

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
    int n = 16;
    int m = 16;

    // initialize matrices and vectors
    z_sp_mat_f A = eye(n);

    z_mat_f v = rand_mat(m, 1);

    // initializes components (indices) of u
    z_mat_f u = zeros(n, 1);
    for(int i = 0; i < n; ++i)
        u.entries[i] = i;


    // initialize the BSP system
    bsp_init("bin/e_spmv.srec", argc, argv);

    int nprocs = bsp_nprocs();
    bsp_begin(nprocs);

    // distribute the matrix
    switch (nprocs)
    {
        case 16:
            N = 4;
            M = 4;
            break;

        case 64:
            N = 8;
            M = 8;
            break;

        default:
            fprintf(stderr, "Unsupported processor count, please add values"
                            "for N and M in the host program.");
            return -1;
    }

    // partition in e.g. equal blocks and write to epiphany
    int chunk = A.nz / nprocs;
    int chunk_v = A.n / nprocs;
    int chunk_u = A.m / nprocs;

    int* row_idx = malloc(A.m * sizeof(int));
    int* col_idx = malloc(A.n * sizeof(int));
    int* v_idx = malloc(chunk_v * sizeof(int));
    int* u_idx = malloc(chunk_u * sizeof(int));

    SPMV_DOWN_TAG tag;
    int tagsize = sizeof(SPMV_DOWN_TAG);
    ebsp_set_tagsize(&tagsize);

    int offset = 0;
    int offset_v = 0;
    int offset_u = 0;
    for (int pid = 0; pid < nprocs; pid++)
    {
        // FIXME: fix size of last chunk
        if (pid == nprocs - 1) {
            chunk = A.nz % (A.nz / nprocs);
            if (chunk == 0)
                chunk = A.nz / nprocs;

            chunk_v = A.n % (A.n / nprocs);
            if (chunk_v == 0)
                chunk_v = A.n / nprocs;

            chunk_u = A.m % (A.m / nprocs);
            if (chunk_u == 0)
                chunk_u = A.m / nprocs;
        }
        
        tag = TAG_ROWS;
        ebsp_send_down(pid, &tag, &A.n, sizeof(int));

        tag = TAG_COLS;
        ebsp_send_down(pid, &tag, &A.m, sizeof(int));

        tag = TAG_MAT;
        ebsp_send_down(pid, &tag, &A.entries[offset], chunk * sizeof(float));

        // construct row_idx and col_idx
        for(int i = 0; i < A.m; ++i)
            row_idx[i] = 0;

        for(int j = 0; j < A.n; ++j)
            col_idx[j] = 0;

        for (int i = 0; i < chunk; ++i) {
            if (offset + i > A.nz)
                break;

            row_idx[A.entries[offset + i].i] = 1;
            col_idx[A.entries[offset + i].j] = 1;
        }

        int r = 0;
        for (int i = 0; i < A.n; ++i)
            if (row_idx[i] != 0)
                row_idx[r++] = i;

        int c = 0;
        for (int j = 0; j < A.m; ++j)
            if (col_idx[j] != 0)
                col_idx[c++] = j;

        tag = TAG_ROW_IDX;
        ebsp_send_down(pid, &tag, &row_idx[0], r * sizeof(int));

        tag = TAG_COL_IDX;
        ebsp_send_down(pid, &tag, &col_idx[0], c * sizeof(int));

        tag = TAG_MAT_INC;
        // FIXME: INC should be local indices
        ebsp_send_down(pid, &tag, &A.inc[offset], chunk * sizeof(float));

        for (int i = 0; i < chunk_v; ++i)
            v_idx[i] = offset_v + i;

        for (int j = 0; j < chunk_u; ++j)
            u_idx[j] = offset_u + j;

        tag = TAG_V_IDX;
        ebsp_send_down(pid, &tag, &v_idx[0], chunk_v * sizeof(int));

        tag = TAG_U_IDX;
        ebsp_send_down(pid, &tag, &u_idx[0], chunk_u * sizeof(int));

        tag = TAG_V_VALUES;
        ebsp_send_down(pid, &tag, &v.entries[offset_v], chunk_v * sizeof(float));

        offset += chunk;
        offset_v += chunk_v;
        offset_u += chunk_u;
    }

    free(u_idx);
    free(v_idx);
    free(col_idx);
    free(row_idx);

    ebsp_spmd();

    // read result
    tagsize = sizeof(int);
    ebsp_set_tagsize(&tagsize);

    int nmsgs = -1;
    int nbytes = -1;

    ebsp_qsize(&nmsgs, &nbytes);

    int status = -1;
    int up_tag = -1;
    for (int i = 0; i < nmsgs; i++)
    {
        ebsp_get_tag(&status, &up_tag);
        printf("up_tag: %i\n", up_tag);
        ebsp_move(&u.entries[up_tag], sizeof(float));
    }

    bsp_end();

    printf("v:\n");
    z_mat_pretty_print(v);

    printf("u:\n");
    z_mat_pretty_print(u);

    z_sp_mat_destroy(&A);
    z_mat_destroy(&v);
    z_mat_destroy(&u);

    return 0;
}
