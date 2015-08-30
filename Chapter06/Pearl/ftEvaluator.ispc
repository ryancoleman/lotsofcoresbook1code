#include "ftTypes.h"

#define SOA_LAYOUT

export uniform int8* uniform ispc_alloc(uniform size_t nbytes)
{
    return uniform new uniform int8[nbytes];
}

export void evaluate(const uniform size_t cPacket,
                     const uniform size_t compCount,
                     const uniform int code[],
                           uniform float *uniform compVals,
                           uniform float *uniform stack,
                           uniform float *uniform results,
                           uniform int8  *uniform status)
{
    uniform int sp;     // stack pointer
    uniform int index;  // system/component index
    uniform int pc = 0; // program counter (index into code array)
    uniform bool stop;  // flag to break out of while loop

    uniform float lvalues[PACKET_SIZE];
    uniform float rvalues[PACKET_SIZE];

    sp = -1;
    pc = 0;
    stop = false;

    foreach (i = 0 ... PACKET_SIZE) {
        status[i] = GoodStatus;
    }

    cwhile (true) {
        switch (code[pc++]) {
        case CompType:
            index = code[pc++];
            ++sp;
            #ifdef SOA_LAYOUT
            foreach (i = 0 ... PACKET_SIZE) {
                stack[(PACKET_SIZE*sp)+i] = compVals[(PACKET_SIZE*index)+i];
            }
            #else
            foreach (i = 0 ... PACKET_SIZE) {
                stack[(PACKET_SIZE*sp)+i] = compVals[(numComps*i)+index];
            }
            #endif
            break;

        case ConstType:
            // push constant value in next code word on stack as float
            ++sp;
            foreach (i = 0 ... PACKET_SIZE) {
                stack[(PACKET_SIZE*sp)+i] = *((float *)&code[pc]);
            }
            pc++;
            break;

        case AndType:
            foreach (i = 0 ... PACKET_SIZE) {
                rvalues[i] = stack[(PACKET_SIZE*sp)+i];
                lvalues[i] = stack[(PACKET_SIZE*(sp-1))+i];
                cif (lvalues[i] != 0.0f) {
                    stack[(PACKET_SIZE*(sp-1))+i] =
                            lvalues[i]*rvalues[i];
                }
            }
            sp--;

            break;

        case OrType:  // assumes statistical independence
            foreach (i = 0 ... PACKET_SIZE) {
                rvalues[i] = stack[(PACKET_SIZE*sp)+i];
                lvalues[i] = stack[(PACKET_SIZE*(sp-1))+i];
                stack[(PACKET_SIZE*(sp-1))+i] = lvalues[i] + rvalues[i] -
                                                 lvalues[i] * rvalues[i];
            }
            sp--;
            break;

        case EndType:
            // end of code, check stack pointer and return
            foreach (i = 0 ... PACKET_SIZE) {
                cif (sp != 0) {
                    status[i] = StackErrStatus;
                } else {
                    results[i] = stack[0+i];
                }
            }
            stop = true; // exit while loop
            break;

        default:
            foreach (i = 0 ... PACKET_SIZE) {
                status[i] = BadCodeStatus;
            }
            stop = true; // exit while loop
            break;
        }

        cif (stop) break;
    }

}
