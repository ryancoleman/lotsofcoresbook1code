/*
  Helper files
  (C) Jason Sewall : Intel
*/
/*

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <getopt.h>
#include "timeseries.hpp"
#include <unistd.h>
static const char usage_str[] = "USAGE:\t%s [-s start] [-e end] [-l] [-h] idxfile1 idxfile2\n";

static void usage(const char *name)
{
    die(usage_str, basename(name));
}

static void help(const char *name)
{
    fprintf(stderr, usage_str, name);
    fprintf(stderr, "DESCRIPTION\n"
            "\t Compare timseries results with vortex particles\n");
    fprintf(stderr, "OPTIONS\n"
            "\t-s,--start <frame>\n\t    Start at frame <frame> (default 1)\n"
            "\t-e,--end <frame>\n\t    Run up to <frame> (not inclusive, defaults to end of shorter)>\n"
            "\t-l,--last-only\n\t    Only test last common frame in inputs (ignores -s and -e options)\n"
            "\t-h,--help\n\t    print this help message\n"
            );
}

bool compare(int nx,
             int ny,
             double *l2,
             double *linf,
             int    *linfarg,
             int numstream,
             int frameno,
             const void *frame1, const void *frame2)
{
    const int stride = nx*ny;

    for(int s = 0; s < numstream; ++s)
    {
        l2[s]               = 0.0;
        linf[s]             = 0.0;
        linfarg[s]          = -1;
        const double *base1 = ((const double*)frame1) + s*stride;
        const double *base2 = ((const double*)frame2) + s*stride;

        for(int i = 0; i < stride; ++i)
        {
            l2[s]   += (base1[i] - base2[i])*(base1[i] - base2[i]);
            if(linf[s] < std::abs(base1[i] - base2[i]))
            {
                linf[s]    = std::abs(base1[i] - base2[i]);
                linfarg[s] = i;
            }
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    int  start_frame = 0;
    int  end_frame   = -1;
    bool last_only   = false;
    option opts[]    =
    {
        {"start",              required_argument, 0, 's'},
        {"end",                required_argument, 0, 'e'},
        {"last-only",          required_argument, 0, 'l'},
        {"help",               false,             0, 'h'},
        {0,                    0,                 0,   0},
    };

    int opt;
    while((opt = getopt_long(argc, argv, "s:e:lh", opts, 0)) != -1)
    {
        switch(opt)
        {
        case 0:
            break;
        case 's':
            start_frame = atoi(optarg);
            if(start_frame < 0)
                die("--[s]start is %d, must be >= 0\n", start_frame);
            break;
        case 'e':
            end_frame = atoi(optarg);
            if(end_frame < 0)
                die("--[e]nd is %d, must be >= 0\n", end_frame);
            break;
        case 'l':
            last_only = true;
            break;
        case 'h':
            help(argv[0]);
            exit(0);
        default:
            usage(argv[0]);
        }
    }

    if(optind >= argc + 1)
        die("Expected 2 arguments (index files) after options\n");

    timeseries_reader idx1;

    if(!idx1.load(argv[optind]))
        die("Can't load (first) index file %s\n", argv[optind]);

    timeseries_reader idx2;

    if(!idx2.load(argv[optind + 1]))
        die("Can't load (second) index file %s\n", argv[optind + 1]);

    int last_frame = std::min(idx1.frames_n, idx2.frames_n);
    if(end_frame != -1)
        last_frame = std::min(last_frame, end_frame);

    if(last_only)
        start_frame = last_frame-1;
    for(int current_frame = start_frame; current_frame < last_frame; ++current_frame)
    {
        double time1;
        size_t size1;
        const void *fr1 = idx1.get_frame(current_frame, &time1, &size1);
        if(!fr1)
            die("Woah, couldn't get frame %d from idx1", current_frame);

        double time2;
        size_t size2;
        const void *fr2 = idx2.get_frame(current_frame, &time2, &size2);
        if(!fr2)
            die("Woah, couldn't get frame %d from idx2", current_frame);

        if(time1 != time2)
            die("Frame %d: times differ! (First = %le, second = %le, first-second = %le)\n", current_frame, time1, time2, time1-time2);

        if(size1 != size2)
            die("Frame %d: sizes differ! (First = %zu, second = %zu)\n", current_frame, size1, size2);

        double l2[4];
        double linf[4];
        int    linfarg[4];

        size_t size;
        const int nx1 = ((const int*)(idx1.get_static("nx", &size)))[0];
        const int ny1 = ((const int*)(idx1.get_static("ny", &size)))[0];

        const int nx2 = ((const int*)(idx2.get_static("nx", &size)))[0];
        const int ny2 = ((const int*)(idx2.get_static("ny", &size)))[0];

        if(nx1 != nx2)
            die("Differing grid x dimensions! (First = %d, second = %d)\n", nx1, nx2);
        if(ny1 != ny2)
            die("Differing grid y dimensions! (First = %d, second = %d)\n", ny1, ny2);

        compare(nx1, ny1, l2, linf, linfarg, 4, current_frame, fr1, fr2);
        fprintf(stderr, "Frame %d\n", current_frame);
        for(int s = 0; s < 4; ++s)
            fprintf(stderr, "          %d err:l2 = %le linf = %le (inf @ %d)\n", s, std::sqrt(l2[s]), linf[s], linfarg[s]);

    }
    if(idx1.frames_n - last_frame > 0)
        fprintf(stderr, "[warning] First has %d more frames unchecked\n", idx1.frames_n - last_frame);

    if(idx2.frames_n - last_frame > 0)
        fprintf(stderr, "[warning] Second has %d more frames unchecked\n", idx2.frames_n - last_frame);

    return EXIT_SUCCESS;
}
