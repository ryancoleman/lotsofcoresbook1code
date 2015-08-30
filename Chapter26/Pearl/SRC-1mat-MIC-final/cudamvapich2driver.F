      PROGRAM CUDAMVAPICH2DRIVER
*
*  Internal Parameters
*  ===================
*
*  TOTMEM   INTEGER, default = 2000000
*           TOTMEM is a machine-specific parameter indicating the
*           maximum amount of available memory in bytes.
*           The user should customize TOTMEM to his platform.  Remember
*           to leave room in memory for the operating system, the BLACS
*           buffer, etc.  For example, on a system with 8 MB of memory
*           per process (e.g., one processor on an Intel iPSC/860), the
*           parameters we use are TOTMEM=6200000 (leaving 1.8 MB for OS,
*           code, BLACS buffer, etc).  However, for PVM, we usually set
*           TOTMEM = 2000000.  Some experimenting with the maximum value
*           of TOTMEM may be required.
*
*  INTGSZ   INTEGER, default = 4 bytes.
*  DBLESZ   INTEGER, default = 8 bytes.
*           INTGSZ and DBLESZ indicate the length in bytes on the
*           given platform for an integer and a double precision real.
*  MEM      DOUBLE PRECISION array, dimension ( TOTMEM / DBLESZ )
*
*           All arrays used by SCALAPACK routines are allocated from
*           this array and referenced by pointers.  The integer IPA,
*           for example, is a pointer to the starting element of MEM for
*           the matrix A.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            DBLESZ, INTGSZ, MEMSIZ
      INTEGER            N, NB, NP, IPA
      integer*8 TOTMEM
      parameter(TOTMEM=8*1024*1024*1024)
      PARAMETER          ( DBLESZ = 8, INTGSZ = 4,
     $                     MEMSIZ = TOTMEM / DBLESZ)
      INTEGER            IAM,NPROCS

*     .. Local Arrays ..
      INTEGER            DESCA( DLEN_ ), IERR( 1 )
      DOUBLE PRECISION,allocatable,dimension(:) ::     MEM

*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT, 
     $                   BLACS_GRIDINFO, DESCINIT, BLACS_BARRIER, 
     $                   BLACS_GRIDEXIT, BLACS_EXIT

*      write (*,*), "first line"

      allocate( MEM( TOTMEM/DBLESZ ) )

*      write (*,*), "second line"

      CALL BLACS_PINFO( IAM, NPROCS )
*      write (*,*), "third line"
      call cublasInit(IAM,NPROCS)

*      write (*,*), "after blacs_pinfo cublasinit"

      NPROW = 2
      NPCOL = 2
*
*        Define process grid
*
         CALL BLACS_GET( -1, 0, ICTXT )
         CALL BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
         CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

*      write (*,*), "after blacs_get  blacs_gridinit blacs_gridinfo"

      memsize = 1000
      N = 100
      NB = 2 
      NP = 50 
      IERR( 1 ) = 0
      IPA = 1
 

*
*              Initialize the array descriptor for the matrix A on host men(ipa)
*
               CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, ICTXT,
     $                        NP, IERR( 1 ) )

*      write (*,*), "after descinit"

               CALL BLACS_BARRIER( ICTXT, 'All' )

*      write (*,*), "after descinit blacs_barrier"

               call pdcopymatrix( MEM(ipa),1,1,descA)


      write (*,*), "after descinit blacs_barrier  pdcopymatrix"


*      call cublasShutdown()
      CALL BLACS_GRIDEXIT( ICTXT )
      CALL BLACS_EXIT( 0 )
*      call cublasShutdown()

      write (*,*), "after blacs_gridexit blacs_exit cublasshutdown"

      stop
      end






