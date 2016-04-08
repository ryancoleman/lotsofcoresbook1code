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
#include "timeseries.hpp"
#include <cstring>
#include <cstdarg>
#include <cstdlib>
#include <sys/mman.h>
#include <algorithm>
#include <sys/stat.h>
#include <cerrno>
#include <cassert>

#undef EXTEND_ARRAY

inline void xdie(const char *fmt, ...)
{
    va_list val;
    va_start(val, fmt);
    vfprintf(stderr, fmt, val);
    va_end(val);
    exit(EXIT_FAILURE);
}

static void make_path(const char *str)
{
    char buff[1024];
    char *current = buff;
    while(*str)
    {
        *current = *str;
        ++current;
        if(current - buff >= 1023)
            xdie("Prefix path is too long (allow 1023, got %d\n", current-buff);
        if(*str == '/')
        {
            *current = 0;
            int dirres = mkdir(buff,  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if(dirres == -1 && errno != EEXIST)
            {
                perror("mkdir");
                xdie("mkdir failed!");
            }
        }
        ++str;
    }
}

bool timeseries_writer::initialize(const char *in_prefix, size_t mps)
{
    max_pack_size = mps;

    make_path(in_prefix);

    char buff[1024];
    snprintf(buff, 1023, "%s.idx", in_prefix);

    index_file = fopen(buff, "w");
    if(!index_file)
        return false;

    snprintf(buff, 1023, "%s%%05d.pak", in_prefix);
    format_string = strdup(buff);

    current_pack      = 0;
    current_pack_name = 0;
    pack_no           = -1;

    current_entry_str[0] = 0;

    return true;
}

void timeseries_writer::finish()
{
    if(current_pack)
    {
        fprintf(index_file, "%s %zu\n", current_entry_str, current_pack_size);
        fclose(current_pack);
    }
    fclose(index_file);

    free(format_string);
    free(current_pack_name);
}

static bool flush_frame(timeseries_writer *tw)
{
    if(tw->current_entry_str[0])
    {
        assert(tw->current_pack);
        fprintf(tw->index_file, "%s %zu\n", tw->current_entry_str, tw->current_pack_size);
        fflush(tw->index_file);
        tw->current_entry_str[0] = 0;
        return true;
    }
    return false;
}

bool timeseries_writer::comment(const char *str)
{
    flush_frame(this);

    if(*str)
        fputs("# ", index_file);

    for(; *str; ++str)
    {
        fputc(*str, index_file);
        if(*str == '\n')
        {
            if(*(str + 1))
                fputs("# ", index_file);
            else
                return true;
        }
    }
    fputc('\n', index_file);
    fflush(index_file);
    return true;
}

static bool check_file(timeseries_writer *tw, size_t size_hint)
{
    flush_frame(tw);

    if(!tw->current_pack || (size_hint < tw->max_pack_size && tw->current_pack_size + size_hint > tw->max_pack_size))
    {
        if(tw->current_pack)
        {
            fclose(tw->current_pack);
            free(tw->current_pack_name);
        }

        char buff[1024];
        snprintf(buff, 1023, tw->format_string, ++tw->pack_no);
        tw->current_pack_name = strdup(basename(buff));
        tw->current_pack      = fopen(buff, "wb");
        if(!tw->current_pack)
            return false;
        tw->current_pack_size = 0;
    }
    return true;
}

bool timeseries_writer::new_frame(double t, size_t size_hint)
{
    if(!check_file(this, size_hint))
        return false;

    assert(current_entry_str[0] == 0);
    snprintf(current_entry_str, 1023, "f %20.14lf %s %zu", t, current_pack_name, current_pack_size);
    return true;
}

bool timeseries_writer::new_static(const char *name, size_t size_hint)
{
    if(!check_file(this, size_hint))
        return false;

    assert(current_entry_str[0] == 0);
    snprintf(current_entry_str, 1023, "s %s %s %zu", name, current_pack_name, current_pack_size);
    return true;
}

size_t timeseries_writer::append(const void *data, size_t data_size)
{
    if(current_entry_str[0] == 0)
        return 0;
    size_t wrote = fwrite(data, data_size, 1, current_pack);
    current_pack_size += wrote*data_size;
    return wrote*data_size;
}

bool timeseries_reader::load(const char *index_filename)
{
    index_file = fopen(index_filename, "r");
    if(!index_file)
        return false;
    struct stat st;
    int stat_res = fstat(fileno(index_file), &st);
    if(stat_res != 0)
        return false;
    if(!(S_ISREG(st.st_mode) || S_ISLNK(st.st_mode)))
        return false;

    const char *back = strrchr(index_filename, '/');
    prefix           = back ? strndup(index_filename, back-index_filename+1) : strdup("");

    frames          = 0;
    frames_n        = 0;
    frames_n_allocd = 0;

    statics          = 0;
    statics_n        = 0;
    statics_n_allocd = 0;

    files          = 0;
    files_n        = 0;
    files_n_allocd = 0;

    refresh();

    return true;
}

#define EXTEND_ARRAY(name, num, n_allocd)       \
    if(name##_n + num >= n_allocd)              \
    {                                           \
        n_allocd = (name##_n + num)*2;          \
        void *m  = realloc(name, sizeof(name[0])*n_allocd); \
        name     = (typeof(name)) m;            \
    }

static bool read_frame(timeseries_reader *tsr, char *file, size_t *low_offset, size_t *high_offset)
{
    EXTEND_ARRAY(tsr->frames, 1, tsr->frames_n_allocd);
    int num_read = fscanf(tsr->index_file, "%lf %1023s %zu %zu", &tsr->frames[tsr->frames_n].t, file, &tsr->frames[tsr->frames_n].start_offset, &tsr->frames[tsr->frames_n].end_offset);
    if(num_read == 4)
    {
        *low_offset   = tsr->frames[tsr->frames_n].start_offset;
        *high_offset  = tsr->frames[tsr->frames_n].end_offset;
        return true;
    }

    return false;
}

static bool read_static(timeseries_reader *tsr, char *file, size_t *low_offset, size_t *high_offset)
{
    EXTEND_ARRAY(tsr->statics, 1, tsr->statics_n_allocd);
    char buff[1024];
    int num_read = fscanf(tsr->index_file, "%1023s %1023s %zu %zu", buff, file, &tsr->statics[tsr->statics_n].start_offset, &tsr->statics[tsr->statics_n].end_offset);
    if(num_read == 4)
    {
        tsr->statics[tsr->statics_n].name = strdup(buff);
        *low_offset                       = tsr->statics[tsr->statics_n].start_offset;
        *high_offset                      = tsr->statics[tsr->statics_n].end_offset;
        return true;
    }

    return false;
}

int timeseries_reader::refresh()
{
    char buff[1024];
    int  nread     = 0;
    int  back_file = files_n;
    while(!feof(index_file))
    {
        off_t last_offs = ftello64(index_file);

        size_t  low_offset;
        size_t  high_offset;
        int    *file_no;
        char    current = fgetc(index_file);
        bool    reset = false;
        switch(current)
        {
        case '#':
            while(!feof(index_file) && current != '\n')
                current = fgetc(index_file);
            if(current == '\n')
                continue;
            reset = true;
            break;
        case 'f':
            {
                bool okay = read_frame(this, buff, &low_offset, &high_offset);
                if(okay)
                {
                    file_no = &(frames[frames_n].file_no);
                    ++frames_n;
                    ++nread;
                }
                reset = !okay;
            }
            break;
        case 's':
            {
                bool okay = read_static(this, buff, &low_offset, &high_offset);
                if(okay)
                {
                    file_no = &(statics[statics_n].file_no);
                    ++statics_n;
                    ++nread;
                }
                reset = !okay;
            }
            break;
        }

        if(reset)
        {
            fseeko64(index_file, last_offs, SEEK_SET);
            break;
        }

        if(!files_n || strcmp(files[files_n-1].name, buff) != 0)
        {
            EXTEND_ARRAY(files, 1, files_n_allocd);
            files[files_n].name           = strdup(buff);
            files[files_n].fp             = 0;
            files[files_n].map_bytes      = 0;
            files[files_n].map_base       = 0;
            files[files_n].lowest_offset  = low_offset;
            files[files_n].highest_offset = high_offset;
            ++files_n;
        }
        else
        {
            files[files_n-1].lowest_offset  = std::min(files[files_n-1].lowest_offset,  low_offset);
            files[files_n-1].highest_offset = std::max(files[files_n-1].highest_offset, high_offset);
        }
        // currently, we assume that files appear in strictly increasing order in the index file
        *file_no = files_n-1;
    }

    if(back_file)
    {
        if(files[back_file-1].map_bytes != files[back_file-1].highest_offset)
        {
            munmap(files[back_file-1].map_base, files[back_file-1].map_bytes);
            void *new_map = mmap(files[back_file-1].map_base, files[back_file-1].highest_offset,  PROT_READ, MAP_PRIVATE, fileno(files[back_file-1].fp), 0);
            if(new_map == (void*)-1)
            {
                perror("mmap");
                xdie("Couldn't mmap file: %s\n", files[back_file-1].name);
            }
            files[back_file-1].map_bytes = files[back_file-1].highest_offset;
        }
    }

    for(int fi = back_file; fi < files_n; ++fi)
    {
        char buff[1024];
        snprintf(buff, 1023, "%s%s", prefix, files[fi].name);
        files[fi].fp        = fopen(buff, "r");
        if(!files[fi].fp)
        {
            perror("fopen");
            xdie("Couldn't open file: %s\n", buff);
        }

        files[fi].map_bytes = files[fi].highest_offset;
        files[fi].map_base  = mmap(0, files[fi].map_bytes, PROT_READ, MAP_PRIVATE, fileno(files[fi].fp), 0);
        if(files[fi].map_base == (void*)-1)
        {
            perror("mmap");
            xdie("Couldn't mmap file: %s\n", buff);
        }
    }

    return nread;
}

const void *timeseries_reader::get_frame(int frameno, double *t, size_t *size) const
{
    const frame_entry *fr     = frames + frameno;
    const file_entry  *fi     = files + fr->file_no;
    *t                        = fr->t;
    *size                     = fr->end_offset - fr->start_offset;
    return (const char*)fi->map_base + fr->start_offset;
}

const void *timeseries_reader::get_static(int staticno, const char **name, size_t *size) const
{
    const static_entry *st = statics + staticno;
    const file_entry   *fi = files + st->file_no;
    *name                  = st->name;
    *size                  = st->end_offset - st->start_offset;
    return (const char*)fi->map_base + st->start_offset;
}

const void *timeseries_reader::get_static(const char *name, size_t *size) const
{
    const char *outname;
    for(int i = 0; i < statics_n; ++i)
        if(strcmp(statics[i].name, name) == 0)
            return get_static(i, &outname, size);

    *size = 0;
    return 0;
}
