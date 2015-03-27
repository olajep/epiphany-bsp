/*
File: e_spmd.h

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

#include <e_bsp.h>
#include "e-lib.h"

#include <math.h>
#include "common.h"

// Function for switch from col major to row major
int row_major_cmp(const void* lhs, const void* rhs)
{
    int diff_i = (z_triplet*)lhs->i - (z_triplet*)rhs->i;
    if (diff_i == 0)
        return (z_triplet*)lhs->j - (z_triplet*)rhs->j;
    else
        return diff_i
}

// SPARSE MATRIX-VECTOR MULTIPLICATION
// -----------------------------------
// Algorithm:
// ==========
//
// The algorithm consists of 4 steps. 
// (1) Obtain non-local vector components v_j
//      (a) number of non-local js
//      (b) indices of non-local js
//      (c) owners of non-local js
//      (d) values of non-local js
// (2) Compute contributions (u_i)_s
// (3) Send contributions (u_i)_s to owner of u_i
// (4) Accumulate contributions and obtain vector u = Av
//
// Memory:
// =======
//
// Note that since we assume a completely arbitrary distribution, and we 
// will allocate the memory statically we have the following worst-
// cases:
// - entries: O(nz) * 3
// -       v: O(m) * 2
// -       u: O(n) * 2
// -   v_own: O(m / p)
// -   u_own: O(n / p)
// - (u_i)_s: O(np)
//
// Especially the last one is troubling, but can be solved by using
// a message passing system, dynamic allocation or limiting the load-
// imbalance epsilon. For now we just make the arrays as large as they
// need to be for the worst case. Later we will provide a real spmv 
// that is optimized and can be used in real applications.

typedef struct
{
    int i;
    float value;
} z_pair;

typedef struct
{
    int i;
    int j; 
    float value;
} z_triplet;

int main()
{
    bsp_begin();

    int nprocs = bsp_nprocs(); 
    int p = bsp_pid();

    // information about distribution
    // note: these are defined for the local (sub)objects
    int nz = 0;
    int nv = 0;
    int nu = 0;

    // information on matrix size
    int rows = 0;
    int cols = 0;

    // FIXME: OBTAIN FROM BSP MESSAGE FROM ARM
    // nz = bsp_move();
    // nz = bsp_move();
    // nz = bsp_move();
    // nz = bsp_move();
    // nz = bsp_move();
    // nz = bsp_move();
    // switch ase
    // ...

    // FIXME: think about rounding (cols / nprocs)
    // store lhs vector
    z_pair* v_pairs = malloc(nv * sizeof(z_pair));
    int* v_owners = malloc(cols / nprocs);

    // store matrix
    z_triplet* a_triplets = malloc(nz * sizeof(z_triplet));

    // store rhs result vector
    z_pair* u_pairs = malloc(nu * sizeof(z_pair));
    int* u_owners = malloc(rows / nprocs);

    //-------------------------------------------------------------------------
    // (1) OBTAIN v_j IF A_ij =/= 0 for some i
    // note that we can loop over the vector simultaneously to check since 
    // triplets are stored in *column major order* we can thus assume it is
    // sorted by j
 
    int v_get_count = 0;
    int vk = 0;
    int* v_get = 0;

    // (a) we see how many vector components we have to obtain from remote loc
    for (int i = 0; i < nz; ++i)
    {
        if (i != (nz - 1) && a_triplets[i].j == a_triplets[i + 1].j)
            continue;

        // we know that this i represents the last one in a column
        // we check if aij and vk are identical
        while (v_pairs[vk].i < a_triplets[i].j)
            ++vk;

        // if we dont have the vector component,
        // we append ais[i] to local v_idxs en obtain it later
        if (v_pairs[vk].i != a_triplets[i].j)
            ++v_get_count;
    }

    v_get = malloc(v_get_count * sizeof(z_pair));
    vk = 0;
    int k = 0;

    // (b) we set the indices of the vector components that we have to obtain
    // FIXME: put in single loop with (a)?
    for (int i = 0; i < nz; ++i)
    {
        if (i != (nz - 1) && a_triplets[i].j == a_triplets[i + 1].j)
            continue;

        while (v_pairs[vk].i < a_triplets[i].j)
            ++vk;

        if (v_pairs[vk].i != a_triplets[i].j)
            v_get[k++].i = ja_triplets[i].j;
    }

    int* v_get_owners = malloc(v_get_count * sizeof(int));
    bsp_push_reg((void*)v_get_owners, v_get_count * sizeof(int));
    bsp_sync();

    // (c) obtain owners of the indices which we need to obtain
    for (int i = 0; i < v_get_count; ++i)
        bsp_get(v_get[i].i % nprocs,
                v_get_owners,
                i * sizeof(int),
                sizeof(int));

    // registering here saves one l cost
    bsp_push_reg((void*)v_get_owners, v_get_count * sizeof(z_pair));
    bsp_sync();

    // (d) finally obtain the actual vector components

    // obtain owner of v_{v_idxs[i]} 
    for (int i = 0; i < v_get_count; ++i)
        bsp_get(v_get[i].i % nprocs,
                v_get_owners,
                i * sizeof(z_pair) + offsetof(z_pair, value),
                sizeof(int));
    bsp_sync();
 
    //-------------------------------------------------------------------------
    // (2) COMPUTE (u_i)_s 

    // we now have the u_i we need.. how do we know the number of local rows
    // first we switch to row major ordering (taking O(nz log(nz)) running time
    // but we focus here on reducing memory footprint so we allow this
    qsort(a_triplets, nz, sizeof(z_triplet), row_major_cmp);

    // repeat the same trick as in (a), compute which u_is we have

    //-------------------------------------------------------------------------
    // (3) SEND u_i TO OWNER. HOW DO YOU WE KNOW WHO IS THE OWNER OF u_i
    // DISTRIBUTE OWNERS CYCLICALLY! (VIA ARM!)
    
    // ...
    
    //-------------------------------------------------------------------------
    // (4) COMPUTE u_i

    // ...

    bsp_end();

    return 0;
}
