      double precision function sym_trace_bas_op(r, ang_max, maxf,
     $     n, itype)
C$Id: sym_tr_bs_op.F 19708 2010-10-29 18:04:21Z d3y133 $
      implicit none
c
      integer ang_max, maxf, itype, n
      double precision r(maxf, maxf, 0:ang_max)
c
      integer i
      double precision trace
c
      trace = 0.0d0
      do i = 1, n
         trace = trace + r(i, i, itype)
      enddo
c
      sym_trace_bas_op = trace
c
      end
