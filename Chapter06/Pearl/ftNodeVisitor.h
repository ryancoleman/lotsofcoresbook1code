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
#ifndef FTVISITOR_H
#define FTVISITOR_H

#include <string>
#include <vector>

struct ftNode;


//===========================
//  BASE NODE VISITOR CLASS
//===========================

class ftNodeVisitor
{
public:
    ftNodeVisitor(void);

    virtual void visitNameNode(ftNode *node);
    virtual void visitComponentNode(ftNode *node);
    virtual void visitSystemNode(ftNode *node);
    virtual void visitConstantNode(ftNode *node);
    virtual void visitOperatorNode(ftNode *node);
    virtual void visitOtherNode(ftNode *node);

    bool stopped(void) const;
    bool modified(void) const;
    bool hasError(void) const;
    const std::string &errorMessage(void) const;

    void traverse(ftNode *node);

protected:
    void errorMessageAppend(const std::string &errorMessage);

    bool m_modified;                // node tree has been modified flag
    bool m_stop;                    // set when node wants processing to stop
    std::string m_errorMessage;     // error message
};


//====================================
//  GET ASSIGN INDEXES VISITOR CLASS
//====================================

class ftNamePool;

class ftAssignIndexesNodeVisitor : public ftNodeVisitor
{
public:
    ftAssignIndexesNodeVisitor(const ftNamePool &systemNames,
        ftNamePool &compNames);

    virtual void visitNameNode(ftNode *node);
    virtual void visitOtherNode(ftNode *node);

private:
    const ftNamePool &m_systemNames;    // reference to system names pool
    ftNamePool &m_compNames;            // pointer to component names pool
};


//====================================
//  GENERATE CODE NODE VISITOR CLASS
//====================================

class ftGenerateCodeNodeVisitor : public ftNodeVisitor
{
public:
    ftGenerateCodeNodeVisitor(const ftNamePool &systemNames,
        const ftNamePool &compNames, const std::vector<ftNode *> &systemNodes,
        std::vector<int> &code, int &stackSize, std::vector<int> &index,
        int &componentCount);

    void visitNameNode(ftNode *node);
    void visitComponentNode(ftNode *node);
    void visitSystemNode(ftNode *node);
    void visitConstantNode(ftNode *node);
    void visitOperatorNode(ftNode *node);
    void visitOtherNode(ftNode *node);

private:
    const ftNamePool &m_systemNames;    // reference to system names pool
    const ftNamePool &m_compNames;      // reference to component names pool
    const std::vector<ftNode *> &m_systemNodes;   // system expression trees
    std::vector<int> &m_code;           // reference to code array
    int &m_stackSize;                   // reference to stack size
    std::vector<int> &m_index;          // reference to translation index array
    int &m_componentCount;              // referehce to component count
    int m_stackPointer;                 // local stack pointer
};


#endif // FTVISITOR_H
