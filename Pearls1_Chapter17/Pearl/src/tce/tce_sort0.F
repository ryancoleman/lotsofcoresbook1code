      subroutine tce_sort_0(unsorted,sorted,factor)
c
c $Id: tce_sort0.F 19706 2010-10-29 17:52:31Z d3y133 $
c
      implicit none
      double precision sorted
      double precision unsorted
      double precision factor
      sorted = unsorted * factor
      return
      end
