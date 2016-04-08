c
c     Print header: TDDFT gradient
c     HvD, NG (2012)
c
      subroutine tddft_grad_print_header(luout)
c
      implicit none
c
#include "util.fh"
c
      integer luout
c
      write(luout,*)
      write(luout,*)
      call util_print_centered(luout,'NWChem TDDFT Gradient Module',
     +                           40, .true.)
      write(luout,*)
      write(luout,*)
c
      end 
c $Id: tddft_grad_print_header.F 25750 2014-06-08 07:56:45Z d3y133 $
