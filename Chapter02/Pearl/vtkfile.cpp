/*
  A 2D hydro code
  (C) Romain Teyssier : CEA/IRFU           -- original F90 code
  (C) Pierre-Francois Lavallee : IDRIS      -- original F90 code
  (C) Guillaume Colin de Verdiere : CEA/DAM -- for the C version
  (C) Jason Sewall : Intel -- for the 'pcl-hydro' C++ version
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
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/stat.h>
#include "arch.hpp"

typedef unsigned char byte;

static const char s_CharPlusSign = '+';
static const char s_CharSlash    = '/';

static char SixBitToChar(byte b);
static char *ToBase64(unsigned char *data, int length);

static char SixBitToChar(byte b) {
    char c;
    if (b < 26) {
        c = (char) ((int) b + (int) 'A');
    } else if (b < 52) {
        c = (char) ((int) b - 26 + (int) 'a');
    } else if (b < 62) {
        c = (char) ((int) b - 52 + (int) '0');
    } else if (b == 62) {
        c = s_CharPlusSign;
    } else {
        c = s_CharSlash;
    }
    return c;
}

static char *ToBase64(unsigned char *data, int length) {
    int padding = length % 3;
    int blocks = (length - 1) / 3 + 1;
    size_t lalloc;
    char *s;
    int i;

    if (length == 0)
        return NULL;

    if (padding > 0)
        padding = 3 - padding;

    // lalloc = (blocks * 4 + 1 + 16);
    lalloc = blocks;
    lalloc *= 4;
    lalloc += 17;

    s = (char*)malloc(lalloc);
    if (s == NULL) {
        fprintf(stderr, "Length=%d, blocks=%d lalloc=%ld\n", length, blocks, lalloc);
    }
    assert(s != NULL);

    for (i = 0; i < blocks; i++) {
        bool finalBlock = i == blocks - 1;
        bool pad2 = false;
        bool pad1 = false;
        if (finalBlock) {
            pad2 = padding == 2;
            pad1 = padding > 0;
        }

        int index = i * 3;
        byte b1 = data[index];
        byte b2 = pad2 ? (byte) 0 : data[index + 1];
        byte b3 = pad1 ? (byte) 0 : data[index + 2];

        byte temp1 = (byte) ((b1 & 0xFC) >> 2);

        byte temp = (byte) ((b1 & 0x03) << 4);
        byte temp2 = (byte) ((b2 & 0xF0) >> 4);
        temp2 += temp;

        temp = (byte) ((b2 & 0x0F) << 2);
        byte temp3 = (byte) ((b3 & 0xC0) >> 6);
        temp3 += temp;

        byte temp4 = (byte) (b3 & 0x3F);

        index = i * 4;
        s[index] = SixBitToChar(temp1);
        s[index + 1] = SixBitToChar(temp2);
        s[index + 2] = pad2 ? '=' : SixBitToChar(temp3);
        s[index + 3] = pad1 ? '=' : SixBitToChar(temp4);
    }
    s[blocks * 4] = (byte) 0;
    return s;
}

#define BINARY 1
#undef MPI
static void vtkwpvd(int nout, char *r) {
    char n[1024];
    char vfname[1024];
    int i;
    FILE *vf = NULL;
    char tmp[10];

    vf = fopen("Hydro.pvd", "w");
    fprintf(vf, "<?xml version=\"1.0\"?>\n");
    fprintf(vf, " <VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(vf, "  <Collection>\n");

    for (i = 1; i <= nout; i++) {
        sprintf(tmp, "%06d", i);
        sprintf(n, "Dep/%c%c%c%c", tmp[0], tmp[1], tmp[2], tmp[3]);
        sprintf(n, "%s/%c%c", n, tmp[4], tmp[5]);
        sprintf(vfname, "%s/Hydro_%04d.pvtr", n, i);
        fprintf(vf, "  <DataSet timestep=\"%d\" part=\"0\" file=\"%s\"  name=\"Asmb:FRAME\"/>\n", i, vfname);
    }

    fprintf(vf, " </Collection>\n");
    fprintf(vf, "</VTKFile>\n");
    fclose(vf);
}

static void vtknm(char *n, int me, int nout) {
    char tmp[10];

    sprintf(tmp, "%06d", nout);
    sprintf(n, "Dep");
    if (me == 0) {
        mkdir(n, 0777);
    }
    sprintf(n, "%s/%c%c%c%c", n, tmp[0], tmp[1], tmp[2], tmp[3]);
    if (me == 0) {
        mkdir(n, 0777);
    }
    sprintf(n, "%s/%c%c", n, tmp[4], tmp[5]);

    if (me == 0) {
        mkdir(n, 0777);
    }
}

void vtkfile(int step, const REAL_T *q, const int n[2], const int padding, const int ystride, const int varstride, const double dx) {
    char name[1024];
    char vfrname[1024];
    FILE *fic, *vf;
    int i, j, nv;

    enum {ID = 0, IU = 1, IV = 2, IP = 3};

    // First step : create the directory structure ONLY using PE0
#ifdef MPI
    if (H.nproc > 1) MPI_Barrier(MPI_COMM_WORLD);
#endif
    vtknm(vfrname, 0, step); // create the directory structure
    // if (0 == 0) fprintf(stderr, "%s\n", vfrname);
#ifdef MPI
    if (H.nproc > 1) MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Write a domain per PE
    sprintf(name, "%s/Hydro_%05d_%04d.vtr", vfrname, 0, step);
    fic = fopen(name, "w");
    if (fic == NULL) {
        fprintf(stderr, "Ouverture du fichier %s impossible\n", name);
        exit(1);
    }
    fprintf(fic, "<?xml version=\"1.0\"?>\n");
    fprintf(fic, "<VTKFile type=\"RectilinearGrid\" byte_order=\"LittleEndian\">\n");
    fprintf(fic, " <RectilinearGrid WholeExtent=\" %d %d %d %d %d %d\">\n",
            0, n[0], 0, n[1], 0, 1);
    fprintf(fic, "  <Piece Extent=\" %d %d %d %d %d %d\" GhostLevel=\"0\">\n",
            0, n[0], 0, n[1], 0, 1);
    fprintf(fic, "   <Coordinates>\n");

    fprintf(fic, "    <DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\">\n");
    for (i = 0; i <= n[0]; i++) {
        fprintf(fic, "%f ", i * dx);
    }
    fprintf(fic, "\n");
    fprintf(fic, "    </DataArray>\n");
    fprintf(fic, "    <DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\">\n");
    for (j = 0; j <= n[1]; j++) {
        fprintf(fic, "%f ", j * dx);
    }
    fprintf(fic, "\n");
    fprintf(fic, "    </DataArray>\n");
    fprintf(fic, "    <DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\">\n");
    fprintf(fic, "%f %f\n", 0., 1. * dx);
    fprintf(fic, "    </DataArray>\n");
    fprintf(fic, "   </Coordinates>\n");
    name[0] = 0;
    for (nv = 0; nv <= IP; nv++) {
        if (nv == ID)
            sprintf(name, "%s varID", name);
        if (nv == IU)
            sprintf(name, "%s varIU", name);
        if (nv == IV)
            sprintf(name, "%s varIV", name);
        if (nv == IP)
            sprintf(name, "%s varIP", name);
    }

    // declaration of the variable list
    fprintf(fic, "   <CellData Scalars=\"%s\">\n", name);
    name[0] = 0;
    for (nv = 0; nv <= IP; nv++) {
        if (nv == ID)
            sprintf(name, "varID");
        if (nv == IU)
            sprintf(name, "varIU");
        if (nv == IV)
            sprintf(name, "varIV");
        if (nv == IP)
            sprintf(name, "varIP");

        //Definition of the cell values
#if BINARY == 1
        fprintf(fic,
                "    <DataArray Name=\"%s\" type=\"Float32\" format=\"binary\" encoding=\"base64\" NumberOfComponents=\"1\">\n",
                name);
        {
            // float tuold[h->net_n[0] * h->net_n[1]];
            float *tuold = NULL;
            char *r64;
            size_t p = 0, lst;

            assert((n[0] * n[1]) > 0);
            tuold = (float *) calloc(n[0] * n[1] + 16, sizeof(float));
            assert(tuold != NULL);

            for (j = 0; j < n[1]; j++) {
                for (i = 0; i < n[0]; i++) {
                    tuold[p++] = (float) q[nv * varstride + (j + padding) * ystride + i + padding];
                }
            }
            // Header = size of the following items
            assert(p <= n[0] * n[1]);

            p *= sizeof(float);
            r64 = ToBase64((byte *) & p, sizeof(int));
            lst = strlen(r64);
            fwrite(r64, 1, lst, fic);
            free(r64);
            r64 = ToBase64((byte *) tuold, p);
            lst = strlen(r64);
            fwrite(r64, 1, lst, fic);
            free(r64);
            free(tuold);
        }
#else
        fprintf(fic, "    <DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\" NumberOfComponents=\"1\">\n", name);

        // the image is the interior of the computed domain
            for (j = 0; j < n[1]; j++) {
                for (i = 0; i < n[0]; i++) {
                    fprintf(fic, "%lf ", q[nv * (nt[0]*nt[1]) + (j + padding) * nt[0] + i + padding]);
                }
                fprintf(fic, "\n");
            }
#endif
        fprintf(fic, "    </DataArray>\n");
    }
    fprintf(fic, "   </CellData>\n");
    fprintf(fic, "  </Piece>\n");
    fprintf(fic, " </RectilinearGrid>\n");
    fprintf(fic, "</VTKFile>\n");
    fclose(fic);

    // At this stage we can write VTK containers. Since only one file is
    // necessary even for multiple domains, it has to be written by one
    // PE only.

#ifdef MPI
    if (H.nproc > 1) MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (0 == 0) {
        sprintf(name, "outputvtk_%05d.pvtr", step);
        sprintf(name, "%s/Hydro_%04d.pvtr", vfrname, step);
        vf = fopen(name, "w");
        if (vf == NULL) {
            fprintf(stderr, "Ouverture du fichier %s impossible\n", name);
            exit(1);
        }
        fprintf(vf, "<?xml version=\"1.0\"?>\n");
        fprintf(vf, "<VTKFile type=\"PRectilinearGrid\" byte_order=\"LittleEndian\">\n");
        fprintf(vf, "<PRectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\"  GhostLevel=\"0\" >\n", n[0], n[1], 1);
        fprintf(vf, " <PCellData>\n");
        for (nv = 0; nv <= IP; nv++) {
            name[0] = '\0';
            if (nv == ID)
                sprintf(name, "varID");
            if (nv == IU)
                sprintf(name, "varIU");
            if (nv == IV)
                sprintf(name, "varIV");
            if (nv == IP)
                sprintf(name, "varIP");

#if BINARY == 1
            fprintf(vf,
                    "  <PDataArray Name=\"%s\" type=\"Float32\" format=\"binary\" encoding=\"base64\" NumberOfComponents=\"1\"/>\n",
                    name);
#else
            fprintf(vf, "  <PDataArray Name=\"%s\" type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>\n", name);
#endif
        }
        fprintf(vf, " </PCellData>\n");
        fprintf(vf, " <PCoordinates>\n");
        fprintf(vf, "  <PDataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>\n");
        fprintf(vf, "  <PDataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>\n");
        fprintf(vf, "  <PDataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\"/>\n");
        fprintf(vf, " </PCoordinates>\n");
        for (i = 0; i < 1; i++) {
            // int box[8];
            // memset(box, 0, 8 * sizeof(int));
            //            CalcSubSurface(0, H.n[0], 0, H.n[1], 0, H.nproc - 1, 0, box, i, 0);
            sprintf(name, "Hydro_%05d_%04d.vtr", i, step);
            // fprintf(vf, " <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>\n", box[XMIN_BOX],
            //         box[XMAX_BOX], box[YMIN_BOX], box[YMAX_BOX], 0, 1, name);
            fprintf(vf, " <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>\n", 0, n[0], 0, n[1], 0, 1, name);

        }
        fprintf(vf, "</PRectilinearGrid>\n");
        fprintf(vf, "</VTKFile>\n");
        fclose(vf);

        // We make the time step available only now to ensure consistency
        vtkwpvd(step, "Dep");
    }
}
