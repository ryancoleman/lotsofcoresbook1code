#include "zorder2d.h"

inline int dialate_1(int x)
{
//x &= 0x0000ffff;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
x = (x ^ (x <<  8)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
x = (x ^ (x <<  4)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
x = (x ^ (x <<  2)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
x = (x ^ (x <<  1)) & 0x55555555; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
return x;
}

//
// return index of x,y
inline int zorder2d(int x, int y) {
  return (dialate_1(y)<<1)|dialate_1(x);
}
