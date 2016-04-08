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
#ifndef FTEVALUATOR_H
#define FTEVALUATOR_H

#include <string>
#include "ftParser.h"

#include "ftTypes.h"


typedef __uint8_t ftStatus;

class ftEvaluator
{
public:

    ftEvaluator();

    ~ftEvaluator();

    // Open a fault-tree file
    //
    // Returns true if the file was successfully opened and parsed.
    //
    // NOTE: This will cause the file to be parsed into memory. It will report
    //       an error on stderr if there was an issue opening the file or with
    //       parsing it.
    bool openFile(std::string fileName);

    // Get/set the number of evaluations (fault-tree instances) that will be
    // evaluated when evaluate() is called;
    unsigned getNumberOfEvaluations();
    void     setNumberOfEvaluations(unsigned n);

    // Set pointer to component values array (the size of this array needs
    // to be the number of components times the number of evaluations)
    unsigned getNumberOfComponents();
    void setComponentData(float *compVals);

    // Translate actual component index to compressed component index
    int getCompressedComponentIndex(int componentIndex);

    // Set the system that will be evaluated, returns if the system was actually
    // set ('true' when no error from parser returning code array)
    bool setSystemToEvaluate(const std::string &system);

    // Do the evaluations for each instance of the tree, return buffers of
    // resulting system values and result status (1 per tree evaluation)
    //
    // returns false if number of evaluations zero or evaluations info not set
    //
    // NOTE: client code is responsible for deleting the buffers later
    bool evaluate(float *&results, ftStatus *&status);

private:

    void evaluateCode(const int code[],
        size_t cEval, float *compVals, float *stack, ftStatus *status);

    // Fault tree parser that this evaluator will ask for the evaluation code
    // that will be executed
    ftParser *m_parser;

    // Number of evaluations
    size_t m_numEval;

    // Component values
    float *m_compVars;

    // Status result values
    float *m_status;

    // System to be evaluated
    std::string m_system;

    // Evaluation info
    ftParser::EvaluationInfo *m_evalInfo;
};

#endif // FTEVALUATOR_H
