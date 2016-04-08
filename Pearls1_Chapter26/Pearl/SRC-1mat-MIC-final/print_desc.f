      subroutine  print_desc( descA )
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      integer descA(DLEN_)

      print*,'descA(DTYPE_) ', descA(DTYPE_)
      print*,'descA(CTXT_) ', descA(CTXT_)
      print*,'descA(M_) ', descA(M_)
      print*,'descA(N_) ', descA(N_)
      print*,'descA(MB_) ', descA(MB_)
      print*,'descA(NB_) ', descA(NB_)
      print*,'descA(RSRC_) ', descA(RSRC_)
      print*,'descA(CSRC_) ', descA(CSRC_)
      print*,'descA(LLD_) ', descA(M_)

      return
      end subroutine
