/*
File: e_bsp.h

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

#pragma once

#include "common.h"
#include "_ansi.h"

//const char registermap_buffer_shm_name[] = REGISTERMAP_BUFFER_SHM_NAME;

/** Starts the BSP program.
 */
void bsp_begin();

/** Finalizes and cleans up the BSP program.
 */
void bsp_end();

/** Returns the number of available processors.
 *
 *  @return nprocs: An integer indicating the number of available processors.
 */
int bsp_nprocs();

/** Returns the processor ID of the local core
 *
 *  pid: An integer indicating the id of the local core
 */
int bsp_pid();

/** Time in seconds that has passed since bsp_begin() was called
 *
 *  t: A floating point value indicating the passed time in seconds.
 */
float bsp_time();

float bsp_remote_time();

/** Terminates a superstep, and starts all communication. The computation is 
 *  halted until all communication has been performed.
 */
void bsp_sync();

/** Registers a variable. Takes effect at the next sync.
 */
void bsp_push_reg(const void* variable, const int nbytes);

/** Put a variable to another processor
 * Buffered version: the data in src is saved at the moment of the call
 * and it is transferred to the other core at the next sync
 */
void bsp_put(int pid, const void *src, void *dst, int offset, int nbytes);

/** Unbuffered version of bsp_put. The data is transferred immediately.
 */
void bsp_hpput(int pid, const void *src, void *dst, int offset, int nbytes);

/** Gets a variable from a processor.
 * Buffered version: the communication occurs during the next sync
 * So after calling bsp_get the caller does not get the data untill
 * the next bsp_sync has finished.
 */
void bsp_get(int pid, const void *src, int offset, void *dst, int nbytes);

/** Unbuffered version of bsp_get
 * The data is tranferred immediately and there is no guarantee
 * about the state of the other processor at that moment.
 */
void bsp_hpget(int pid, const void *src, int offset, void *dst, int nbytes);

/** ebsp_message outputs a debug message by sending it to shared memory
 * So that the host processor can output it to the terminal
 * The attributes in this definition make sure that the compiler checks the
 * arguments for errors
 */
void _EXFUN(ebsp_message, (const char *, ...)
               _ATTRIBUTE ((__format__ (__printf__, 1, 2))));

