/*
 $Id: sgi_flush6.c 19707 2010-10-29 17:59:36Z d3y133 $
*/
#include <stdio.h>
void sgi_flush6_(void)
{
  int return_code;
  return_code = fflush(stdout);
  if (return_code == 0) return;
  (void) perror("Error detected in sgi_flush6:");
}
