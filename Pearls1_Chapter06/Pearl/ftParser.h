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
#ifndef FTPARSER_H
#define FTPARSER_H

#include <string>
#include <list>
#include <vector>

#include "ftNode.h"

class ftFileReader;
class ftNamePool;


class ftParser
{
public:
    struct EvaluationInfo {
        std::vector<int> code;      // expression code
        int stackSize;              // size of stack needed
        int componentCount;         // number of components in expression
        std::vector<int> index;     // component index to compressed index
    };

    // Constructor
    ftParser(const char *fileName);

    // Destructor
    ~ftParser();

    bool isReaderReady(void);
    bool parse(void);
    EvaluationInfo *generateCode(const std::string &system);

    // intermediate functions for testing (should be private)
    bool getSystems(void);
    bool assignIndexes(void);

    // accessor functions for testing (should be removed)
    const ftNamePool &systemNames(void) const;
    const std::vector<ftNode *> &systemNodes(void) const;
    const ftNamePool &CompNames(void) const;
    const std::string compName(int index) const;

private:
    struct OperatorType {
        char op;
        NodeType type;
    };

    bool getSystem(void);
    ftNode *getExpression(void);
    ftNode *getOrOperand(void);
    ftNode *getAndOperand(void);
    ftNode *getOperator(OperatorType *operatorType);
    ftNode *getOperand(void);
    std::string getName(void);
    void reportError(const char *message, int positionOffset = 0);


    // Internal parse tree node type /////////////////////////////////////////
    //TODO: add private data
    ftFileReader *m_reader;
    ftNamePool *m_systemNames;
    std::vector<ftNode *> m_systemNodes;
    ftNamePool *m_CompNames;
    int m_level;  // parentheses level
};

#endif // FTPARSER_H
