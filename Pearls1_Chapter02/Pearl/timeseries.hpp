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
#ifndef __TIMESERIES_HPP__
#define __TIMESERIES_HPP__

#include "array-macros.hpp"

struct timeseries_writer
{
    bool   initialize(const char *in_prefix, size_t max_pack_size);
    void   finish();
    size_t append(const void *data, size_t data_size);
    bool   comment(const char *str);
    bool   new_frame(double t, size_t size_hint);
    bool   new_static(const char *name, size_t size_hint);

    char   *format_string;
    size_t  max_pack_size;

    FILE *index_file;

    FILE   *current_pack;
    char   *current_pack_name;
    size_t  current_pack_size;
    int     pack_no;

    char current_entry_str[1024];
};

struct file_entry
{
    char   *name;
    FILE   *fp;
    size_t  map_bytes;
    void   *map_base;
    size_t  lowest_offset;
    size_t  highest_offset;
};

struct frame_entry
{
    double  t;
    size_t start_offset;
    size_t end_offset;
    int    file_no;
};

struct static_entry
{
    char   *name;
    size_t  start_offset;
    size_t  end_offset;
    int     file_no;
};

struct timeseries_reader
{
    bool load(const char *index_file);
    int refresh();

    const void *get_frame(int frameno, double *t, size_t *size) const;
    const void *get_static(int staticno, const char **name, size_t *size) const;
    const void *get_static(const char *name, size_t *size) const;

    char         *prefix;
    FILE         *index_file;
    frame_entry  *frames;
    int           frames_n;
    int           frames_n_allocd;

    static_entry *statics;
    int           statics_n;
    int           statics_n_allocd;
    file_entry   *files;
    int           files_n;
    int           files_n_allocd;
};
#endif
