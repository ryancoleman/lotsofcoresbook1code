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
#include "ftNode.h"
#include "ftNodeVisitor.h"
#include "ftNamePool.h"


//===========================
//  BASE NODE VISITOR CLASS
//===========================

ftNodeVisitor::ftNodeVisitor(void) :
    m_modified(false),
    m_stop(false)
{
}

void ftNodeVisitor::visitNameNode(ftNode *node)
{
    // no default action for this node type
}

void ftNodeVisitor::visitSystemNode(ftNode *node)
{
    // no default action for this node type
}

void ftNodeVisitor::visitComponentNode(ftNode *node)
{
    // no default action for this node type
}

void ftNodeVisitor::visitConstantNode(ftNode *node)
{
    // no default action for this node type
}

void ftNodeVisitor::visitOperatorNode(ftNode *node)
{
    traverse(node->op.lhs);
    if (!stopped()) {
        traverse(node->op.rhs);
    }
}

void ftNodeVisitor::visitOtherNode(ftNode *node)
{
    // no default action for these node types
}

bool ftNodeVisitor::stopped(void) const
{
    return m_stop;
}

bool ftNodeVisitor::modified(void) const
{
    return m_modified;
}

bool ftNodeVisitor::hasError(void) const
{
    return !m_errorMessage.empty();
}

const std::string &ftNodeVisitor::errorMessage(void) const
{
    return m_errorMessage;
}

void ftNodeVisitor::traverse(ftNode *node)
{
    switch (node->type) {
    case NameType:
        visitNameNode(node);
        break;

    case CompType:
        visitComponentNode(node);
        break;

    case SystemType:
        visitSystemNode(node);
        break;

    case ConstType:
        visitConstantNode(node);
        break;

    case AndType:
    case OrType:
        visitOperatorNode(node);
        break;

    default:
        visitOtherNode(node);
    }
}

void ftNodeVisitor::errorMessageAppend(const std::string &errorMessage)
{
    if (!m_errorMessage.empty()) {
        m_errorMessage.append("\n");
    }
    m_errorMessage.append(errorMessage);
}


//=========================================
//  GET ASSIGN INDEXES NODE VISITOR CLASS
//=========================================

ftAssignIndexesNodeVisitor::ftAssignIndexesNodeVisitor(
        const ftNamePool &systemNames, ftNamePool &compNames) :
    ftNodeVisitor(),
    m_systemNames(systemNames),
    m_compNames(compNames)
{
}

void ftAssignIndexesNodeVisitor::visitNameNode(ftNode *node)
{
    int index = m_systemNames.getIndex(node->text);
    if (index == ftNamePool::NotFound) {
        // not a system, so it is a component,
        // add to component name pool and save its index in node
        node->type = CompType;
        node->name.index = m_compNames.addName(node->text);
    } else {
        // a system, so save its index in node
        node->type = SystemType;
        node->name.index = index;
    }
}

void ftAssignIndexesNodeVisitor::visitOtherNode(ftNode *node)
{
    errorMessageAppend("Invalid node found in tree");
    m_stop = true;
}


//====================================
//  GENERATE CODE NODE VISITOR CLASS
//====================================

ftGenerateCodeNodeVisitor::ftGenerateCodeNodeVisitor(
        const ftNamePool &systemNames, const ftNamePool &compNames,
        const std::vector<ftNode *> &systemNodes, std::vector<int> &code,
        int &stackSize, std::vector<int> &index, int &componentCount) :
    m_systemNames(systemNames),
    m_compNames(compNames),
    m_systemNodes(systemNodes),
    m_code(code),
    m_stackSize(stackSize),
    m_index(index),
    m_componentCount(componentCount),
    m_stackPointer(0)
{
    m_stackSize = 0;
    m_index.assign(m_compNames.size(), -1);
    m_componentCount = 0;
}

void ftGenerateCodeNodeVisitor::visitNameNode(ftNode *node)
{
    m_stop = true;
    errorMessageAppend("Error: Name node type in expression tree");
}

void ftGenerateCodeNodeVisitor::visitComponentNode(ftNode *node)
{
    m_code.push_back(CompType);
    int index = m_index[node->name.index];
    if (index == -1) {
        index = m_componentCount++;
        m_index[node->name.index] = index;
    }
    m_code.push_back(index);
    m_stackPointer++;
    if (m_stackPointer > m_stackSize) {
        m_stackSize = m_stackPointer;
    }
}

void ftGenerateCodeNodeVisitor::visitSystemNode(ftNode *node)
{
    traverse(m_systemNodes[node->name.index]);
}

void ftGenerateCodeNodeVisitor::visitConstantNode(ftNode *node)
{
    union {
        int intValue;
        float floatValue;
    } u;

    m_code.push_back(ConstType);
    u.floatValue = node->constant.value;
    m_code.push_back(u.intValue);
    m_stackPointer++;
    if (m_stackPointer > m_stackSize) {
        m_stackSize = m_stackPointer;
    }
}

void ftGenerateCodeNodeVisitor::visitOperatorNode(ftNode *node)
{
    ftNodeVisitor::visitOperatorNode(node);
    m_code.push_back(node->type);
    m_stackPointer--;
}

void ftGenerateCodeNodeVisitor::visitOtherNode(ftNode *node)
{
    m_stop = true;
    errorMessageAppend("Error: Other node type in expression tree");
}
