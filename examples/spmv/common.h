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

typedef enum
{
    // Phase (1)
    TAG_ROWS = 1,       // O(1)
    TAG_COLS = 2,       // O(1)
    TAG_MAT = 3,        // O(nz)
    TAG_MAT_INC = 4,    // O(nz)
    TAG_ROW_IDX = 5,   // O(nrows)
    TAG_COL_IDX = 6,   // O(ncols)
    TAG_U_IDX = 7,     // O(nu)
    TAG_V_IDX = 8,     // O(nv)
    TAG_V_VALUES = 9  // O(nv)
} SPMV_DOWN_TAG;
