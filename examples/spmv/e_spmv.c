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

int main()
{
    bsp_begin();

    int n = bsp_nprocs(); 
    int p = bsp_pid();


    // (1) OBTAIN v_j IF A_ij =/= 0 for some i
    // note that we can loop over the vector simultaneously to check

    // ...
 
    // (2) COMPUTE (u_i)_s 

    // ...

    // (3) SEND u_i TO OWNER. HOW DO YOU WE KNOW WHO IS THE OWNER OF u_i
    // DISTRIBUTE OWNERS CYCLICALLY! (ON ARM?)
    
    // ...
    
    // (4) COMPUTE u_i

    // ...

    bsp_end();

    return 0;
}
