#ifndef SCALAPACK_H
#define SCALAPACK_H 1

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>






#ifndef MIN
#define MIN(x,y)  (((x) < (y)) ? (x) : (y) )
#endif

#ifndef MAX
#define MAX(x,y)  (((x) > (y)) ? (x) : (y) )
#endif

#ifndef MOD
#define MOD(x,y)  ((x) % (y))
#endif


#define REAL_PART 0
#define IMAG_PART 1

#define    BLOCK_CYCLIC_2D     1

#define    DTYPE_             0                   /* Descriptor Type */
#define    CTXT_              1                     /* BLACS context */
#define    M_                 2             /* Global Number of Rows */
#define    N_                 3          /* Global Number of Columns */
#define    MB_                4                 /* Row Blocking Size */
#define    NB_                5              /* Column Blocking Size */
#define    RSRC_              6            /* Starting Processor Row */
#define    CSRC_              7         /* Starting Processor Column */
#define    LLD_               8           /* Local Leading Dimension */
#define    DLEN_              9                 /* Descriptor Length */




#include "scalapack_i.h"
#include "scalapack_z.h"
#include "scalapack_c.h"
#include "scalapack_d.h"
#include "scalapack_s.h"






#endif
