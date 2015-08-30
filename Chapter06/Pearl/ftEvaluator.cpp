/******************************************************************************
 * Copyright (c) 2014, Jefferson Amstutz
 * Copyright (c) 2014, SURVICE Engineering
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  * The names of contributors cannot be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#include <math.h>
#include <iostream>
using std::cerr;
using std::endl;

#include "string.h"

#include "ftEvaluator.h"

namespace ispc
{
    extern "C" void* ispc_alloc(size_t nbytes);
    extern "C" void evaluate(const size_t cPacket,
                             const size_t compCount,
                             const int code[],
                             float    *compVals,
                             float    *stack,
                             float    *results,
                             ftStatus *status);
}

ftEvaluator::ftEvaluator() :
    m_parser(0),
    m_numEval(1),
    m_system(""),
    m_evalInfo(NULL)
{
    /*no-op*/
}

ftEvaluator::~ftEvaluator()
{
    delete m_evalInfo;
    delete m_parser;
}

bool ftEvaluator::openFile(std::string fileName)
{
    m_parser = new ftParser(fileName.c_str());
    if (!m_parser->isReaderReady()) {
        cerr << "Failed to open file " << fileName << endl;
        delete m_parser;
        m_parser = 0;
        return false;
    }

    if (!m_parser->parse()) {
        cerr << "Parsing error!" << endl;
        delete m_parser;
        m_parser = 0;
        return false;
    }

    return true;
}

unsigned ftEvaluator::getNumberOfEvaluations()
{
    return m_numEval;
}

void ftEvaluator::setNumberOfEvaluations(unsigned n)
{
    m_numEval = n > 0 ? n : 1;// exclude '0' as a valid input
}

unsigned ftEvaluator::getNumberOfComponents()
{
   return m_evalInfo == NULL ? 0 : m_evalInfo->componentCount;
}

void ftEvaluator::setComponentData(float *compVals)
{
    m_compVars = compVals;
}

int ftEvaluator::getCompressedComponentIndex(int componentIndex)
{
   return m_evalInfo == NULL ? -1 : m_evalInfo->index[componentIndex];
}

bool ftEvaluator::setSystemToEvaluate(const std::string &system)
{
    delete m_evalInfo;
    m_evalInfo = m_parser->generateCode(system);
    if (!m_evalInfo) {
        return false;  // error already reported
    }
    m_system = system;

    return true;
}

bool ftEvaluator::evaluate(float *&results, ftStatus *&status)
{
    if (m_numEval == 0 || m_evalInfo == NULL) {
        return false;
    }

    // allocate memory
    float *stack =
        (float*)ispc::ispc_alloc(m_numEval*m_evalInfo->stackSize*sizeof(int));

    int *code = (int*)ispc::ispc_alloc(m_evalInfo->code.size()*sizeof(int));
    memcpy(code, m_evalInfo->code.data(), m_evalInfo->code.size()*sizeof(int));

	for (int j = 0; j < 10; ++j) {
#ifndef USE_ISPC
#  ifdef EVAL_THREADING
    #pragma omp parallel for schedule(runtime)
#  endif
    for (size_t i = 0; i < m_numEval; i++) {
        // offset into linear arrays (results, status)
        const long int ao = i;

        // offset into stack array
        const long int so = m_evalInfo->stackSize * i;

        // offset into component array
        const long int co = m_evalInfo->componentCount * i;

        evaluateCode(code, i, m_compVars + co, stack + so, status + ao);

        // retrieve results
        if (status[ao] == GoodStatus) {
            // results on top of stack
            results[ao] = *(stack + so);
        }
    }
#else
    int nPackets = m_numEval / PACKET_SIZE;// # of packets
    int sSize    = m_evalInfo->componentCount * PACKET_SIZE;// size of section
#  ifdef EVAL_THREADING
    #pragma omp parallel for schedule(runtime)
#  endif
    for (size_t i = 0; i < nPackets; i++) {
        // offset into linear arrays (results, status)
        const long int ao = PACKET_SIZE * i;

        // offset into stack array
        const long int so = m_evalInfo->stackSize * i;

        // offset into component array
        const long int co = sSize * i;

        ispc::evaluate(i, m_evalInfo->componentCount, code, m_compVars + co,
                       stack + so, results + ao, status + ao);
    }
#endif
	}

    delete stack;
    delete code;

    return true;
}

void ftEvaluator::evaluateCode(const int code[], size_t cEval,
                               float *compVals, float *stack, ftStatus *status)
{
    int sp;        // stack pointer
    int index;     // system/component index
    float lvalue;  // cache for left value (stack[sp])
    float rvalue;  // cache for right value (stack[sp + 1])

    *status = GoodStatus;
    sp = -1;
    int pc = 0;
    for (;;) {
        switch (code[pc++]) {
        case CompType:
            index = code[pc++];
            stack[++sp] = compVals[index];
            break;

        case ConstType:
            // push constant value in next code word on stack as float
            stack[++sp] = *((float *)&code[pc++]);
            break;

        case AndType:
            rvalue = stack[sp--];
            lvalue = stack[sp];
            if (lvalue != 0.0f) {
                stack[sp] = lvalue * rvalue;
            }
            break;

        case OrType:  // assumes statistical independence
            rvalue = stack[sp--];
            lvalue = stack[sp];
            stack[sp] = lvalue + rvalue - lvalue * rvalue;
            break;

        case EndType:
            // end of code, check stack pointer and return
            if (sp != 0) {
                *status = StackErrStatus;
            }
            return;  // exit for loop

        default:
            *status = BadCodeStatus;
            return;  // exit for loop
        }
    }
}
