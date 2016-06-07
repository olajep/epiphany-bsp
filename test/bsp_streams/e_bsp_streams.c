/*
This file is part of the Epiphany BSP library.

Copyright (C) 2014-2015 Buurlage Wits
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
#include <limits.h>
#include "../common.h"

int main() {
    bsp_begin();

    int s = bsp_pid();

    // Old streaming API

    int* upstream = 0;
    int* upstreamDouble = 0;
    int chunk_size = ebsp_open_up_stream((void**)&upstream, 0);
    ebsp_open_up_stream((void**)&upstreamDouble, 1);
    if (s == 0)
        ebsp_open_up_stream((void**)&upstreamDouble, 1);
    // expect: ($00: BSP ERROR: tried creating opened stream)
    int chunks = 4;

    int* downchunk;
    int* downchunkB;
    int* downchunkDouble;

    ebsp_open_down_stream((void**)&downchunk, 2);
    ebsp_open_down_stream((void**)&downchunkB, 3);
    ebsp_open_down_stream((void**)&downchunkDouble, 4);

    if (s == 0)
        ebsp_open_down_stream((void**)0, 7);
    // expect: ($00: BSP ERROR: stream does not exist)
    if (s == 0)
        ebsp_open_up_stream((void**)0, 7);
    // expect: ($00: BSP ERROR: stream does not exist)
    if (s == 0)
        ebsp_open_down_stream((void**)0, 1);
    // expect: ($00: BSP ERROR: mixed up and down streams)

    for (int i = 0; i < chunks; ++i) {
        ebsp_move_chunk_down((void**)&downchunk, 2, 0);
        ebsp_move_chunk_down((void**)&downchunkB, 3, 0);
        ebsp_move_chunk_down((void**)&downchunkDouble, 4, 1);
        for (int j = 0; j < chunk_size / sizeof(int); ++j) {
            upstream[j] = (i % 2 == 1) ? downchunk[j] : downchunkB[j];
            upstreamDouble[j] = 2 * downchunk[j];
        }
        ebsp_move_chunk_up((void**)&upstreamDouble, 1, 1);
        ebsp_move_chunk_up((void**)&upstream, 0, 0);
    }

    EBSP_MSG_ORDERED("%i", downchunk[0]);
    // expect_for_pid: (3)

    EBSP_MSG_ORDERED("%i", downchunkB[0]);
    // expect_for_pid: (12)

    ebsp_close_up_stream(0);
    ebsp_close_up_stream(1);
    ebsp_close_down_stream(2);
    ebsp_close_down_stream(3);
    ebsp_close_down_stream(4);

    if (s == 0)
        ebsp_close_up_stream(0);
    // expect: ($00: BSP ERROR: tried to close closed stream)
    if (s == 0)
        ebsp_close_down_stream(2);
    // expect: ($00: BSP ERROR: tried to close closed stream)
    if (s == 0)
        ebsp_close_down_stream(7);
    // expect: ($00: BSP ERROR: stream does not exist)
    if (s == 0)
        ebsp_close_up_stream(7);
    // expect: ($00: BSP ERROR: stream does not exist)

    // New streaming API
    int tokensize  = ebsp_stream_open(5);
    int tokensize2 = ebsp_stream_open(6);

    if (tokensize != tokensize2)
        ebsp_message("Invalid token size at ebsp_stream_open");

    // Double buffered upstream
    int* up1 = ebsp_malloc(tokensize);
    int* up2 = ebsp_malloc(tokensize);

    // First stream down from 6 and copy it into 5
    for (;;) {
        int* buffer;
        int size = ebsp_stream_move_down(6, (void**)&buffer, 1);
        if (size == 0)
            break;

        for (int j = 0; j < tokensize / sizeof(int); ++j)
            up1[j] = buffer[j];

        ebsp_stream_move_up(5, up1, size, 0);
        // swap buffers
        int* tmp = up1;
        up1 = up2;
        up2 = tmp;
    }

    // Now stream down from 5, double the values, and copy it into 6
    ebsp_stream_seek(5, INT_MIN); // go back to start
    ebsp_stream_seek(6, INT_MIN); // go back to start
    for (;;) {
        int* buffer;
        int size = ebsp_stream_move_down(5, (void**)&buffer, 1);
        if (size == 0)
            break;

        for (int j = 0; j < tokensize / sizeof(int); ++j)
            up1[j] = 2 * buffer[j];

        ebsp_stream_move_up(6, up1, size, 0);
        // swap buffers
        int* tmp = up1;
        up1 = up2;
        up2 = tmp;
    }

    ebsp_stream_close(5);
    ebsp_stream_close(6);

    ebsp_free(up1);
    ebsp_free(up2);

    bsp_end();

    return 0;
}

