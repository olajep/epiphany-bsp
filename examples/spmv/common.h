/*
File: common.h

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

#define ROWS 50 
#define COLS 50

//-----------------------------------------------------------------------------
// FIXME: EVERYTHING IN THIS SECTION SHOULD BE OBSOLETE 
//-----------------------------------------------------------------------------
#define LOC_MAT_I 0x4000 // O(nz)

// here we assume nz is defined locally (both on host and task)
// FIXME: very ugly implementation, will be fixed by dynamic allocation,
// mailbox system, message passing, or load imbalance constraints.
// We can only store 2048 numbers in one bank..
#define LOC_MAT_J (LOC_MAT_I + nz * sizeof(int)) // O(nz)
#define LOC_MAT_ENTRIES (LOC_MAT_J + nz * sizeof(int)) // O(nz)
#define LOC_NZ (LOC_MAT_ENTRIES + nz * sizeof(float)) // O(1)

#define LOC_V_OWNERS (LOC_NZ + sizeof(int)) //O(m / p)
#define LOC_V_IDXS (LOC_V_OWNERS + ((m / nprocs) + 1) * sizeof(int)) // O(m)
#define LOC_V_ENTRIES (LOC_V_IDXS + m * sizeof(int)) // O(m) 
#define LOC_V_GET_FROM (LOC_V_ENTRIES + m * sizeof(float)) // O(m)
#define LOC_V_COUNT (LOC_V_GET_FROM + m * sizeof(float)) // O(1)

#define LOC_U_IDXS 0
#define LOC_U_ENTRIES 0
#define LOC_U_OWNERS 0
#define LOC_U_COUNT 0
//-----------------------------------------------------------------------------
