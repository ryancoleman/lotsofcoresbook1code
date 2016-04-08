/*
Copyright (c) The University of Tennessee.  All rights reserved.


$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer listed in this license in the documentation and/or other materials provided with the distribution.

- Neither the name of the copyright holders nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

The copyright holders provide no reassurances that the source code provided does not infringe any patent, copyright, or any other intellectual property rights of third parties.  The copyright holders disclaim any liability to any recipient for claims brought against recipient by any third party for infringement of that parties intellectual property rights.
*/
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#define ROOT 0

int pdlltinfo(int* n, int* nb, int* nprow, int* npcol, int* memsize){
    int rank;
    int error;
    int temp[5];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank == ROOT){
        //// #1
        temp[0] = atoi(getenv("PDLLT_N"));
        if(temp[0] <= 0){
            error = 1; goto ERROR;
        }
        //// #2
        temp[1] = atoi(getenv("PDLLT_NB"));
        if(temp[1] <= 0){
            error = 2; goto ERROR;
        }
        //// #3
        temp[2] = atoi(getenv("PDLLT_P"));
        if(temp[2] <= 0){
            error = 3; goto ERROR;
        }
        //// #4
        temp[3] = atoi(getenv("PDLLT_Q"));
        if(temp[3] <= 0){
            error = 4; goto ERROR;
        }
        //// #5
        temp[4] = atoi(getenv("PDLLT_MEMSIZE"));
        if(temp[4] <= 0){
            error = 5; goto ERROR;
        }
        error = 0; 
        MPI_Bcast(&error, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        //// broadcast parameters
        MPI_Bcast(temp, 5, MPI_INT, ROOT, MPI_COMM_WORLD);

    ERROR: // because it look nice
        if(error != 0){
            MPI_Bcast(&error, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            return error;
        }
    }else{
        //// get error status
        MPI_Bcast(&error, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        if(error != 0) return error;
        //// get parameters
        MPI_Bcast(temp, 5, MPI_INT, ROOT, MPI_COMM_WORLD);
    }
    *n = temp[0];
    *nb = temp[1];
    *nprow = temp[2];
    *npcol = temp[3];
    *memsize = temp[4]*1024;
    return 0;
}
