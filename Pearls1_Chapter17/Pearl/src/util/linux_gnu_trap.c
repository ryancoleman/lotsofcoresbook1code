#include <fpu_control.h>
/*---------------------------------------------------------*\
$Id: linux_gnu_trap.c 19707 2010-10-29 17:59:36Z d3y133 $

  This routine turns on floating point traps for
  1) Illegal instructions
  2) Divide by zero
  3) and Floating Point Overflow 
\*---------------------------------------------------------*/
static void __attribute__ ((constructor))
     trapfpe ()
{
    __setfpucw (_FPU_DEFAULT &
      ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM));
}
