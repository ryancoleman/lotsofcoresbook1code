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
#include <iostream>

#include "ftFileReader.h"
#include "ftNode.h"
#include "ftParser.h"
#include "ftNamePool.h"
#include "ftNodeVisitor.h"


ftParser::ftParser(const char *fileName):
    m_reader(new ftFileReader(fileName)),
    m_systemNames(new ftNamePool),
    m_CompNames(new ftNamePool),
    m_level(0)
{

}

ftParser::~ftParser()
{
    delete m_reader;
    delete m_systemNames;
    delete m_CompNames;
    for (unsigned i = 0; i < m_systemNodes.size(); i++) {
        delete m_systemNodes[i];
    }
}


bool ftParser::isReaderReady(void)
{
    return m_reader->isFileOpen();
}


bool ftParser::parse(void)
{
    return getSystems() && assignIndexes();
}


bool ftParser::getSystems(void)
{
    for (;;) {
        m_reader->skipWhiteSpace();
        if (m_reader->getCurrentCharacter() == -1) {
            break;  // end-of-file
        }
        if (!getSystem()) {
            return false;
        }
    }
    return true;
}


bool ftParser::assignIndexes(void)
{
    ftAssignIndexesNodeVisitor assignIndexesNodeVisitor(*m_systemNames,
        *m_CompNames);

    for (int i = 0; i < m_systemNames->size(); i++)
    {
        assignIndexesNodeVisitor.traverse(m_systemNodes[i]);
        if (assignIndexesNodeVisitor.hasError()) {
            std::cerr << "Assign Index Error-System: '"
                << m_systemNames->getName(i) << "' "
                << assignIndexesNodeVisitor.errorMessage();
            return false;
        }
    }
    return true;
}


bool ftParser::getSystem(void)
{
    // need to check for end of input

    std::string systemName = getName();
    if (systemName.empty()) {
        reportError("expected system name", 1);
        return false;
    }
    bool exists;
    int index = m_systemNames->addName(systemName, &exists);
    if (exists) {
        reportError("duplicate system name");
        return false;
    }

    m_reader->skipWhiteSpace();
    char c = m_reader->getCurrentCharacter();
    if (c != '=') {
        reportError("expected '=' after system name", 1);
        return false;
    }
    m_reader->advancePosition();

    ftNode *systemExpression = getExpression();
    if (systemExpression == NULL) {
        return false;  // error already reported
    }

    m_systemNodes.resize(index + 1);
    m_systemNodes[index] = systemExpression;

    // input ready for reading next system
    return true;
}


ftParser::EvaluationInfo *ftParser::generateCode(const std::string &system)
{
    int index = m_systemNames->getIndex(system);
    if (index == ftNamePool::NotFound) {
        std::cerr << "Generate Code Error: system '" << system << "' not found!"
            << std::endl;
        return NULL;
    }

    EvaluationInfo *evaluationInfo = new EvaluationInfo;
    ftGenerateCodeNodeVisitor generateCodeNodeVisitor(*m_systemNames,
        *m_CompNames, m_systemNodes, evaluationInfo->code,
        evaluationInfo->stackSize, evaluationInfo->index,
        evaluationInfo->componentCount);
    generateCodeNodeVisitor.traverse(m_systemNodes[index]);
    if (generateCodeNodeVisitor.hasError()) {
        std::cerr << "Generate Code Error: "
            << generateCodeNodeVisitor.errorMessage() << std::endl;
        delete evaluationInfo;
        return NULL;
    }
    evaluationInfo->code.push_back(EndType);
    return evaluationInfo;
}


ftNode *ftParser::getExpression(void)
{
    static OperatorType operatorType[] = {
        { '|', OrType },
        { ')', ParenType },
        { ';', EndType },
        { '\0' }
    };

    ftNode *lhs = getOrOperand();
    if (lhs == NULL) {
        return NULL;  // error already reported
    }

    for (;;) {
        ftNode *op = getOperator(operatorType);
        if (op == NULL) {
            reportError("expected operator or ';'", 1);
            delete lhs;
            return NULL;  // error already reported
        }
        if (m_level == 0) {
            if (op->type == EndType) {
                delete op;
                return lhs;
            } else if (op->type == ParenType) {
                reportError("no matching opening parentheses");
                delete op;
                delete lhs;
                return NULL;
            }
        } else {  // within parentheses
            if (op->type == EndType) {
                reportError("expected closing parentheses");
                delete op;
                delete lhs;
                return NULL;
            } else if (op->type == ParenType) {
                delete op;
                return lhs;
            }
        }

        ftNode *rhs = getOrOperand();
        if (rhs == NULL) {
            delete op;
            delete lhs;
            return NULL;  // error already reported
        }

        op->op.lhs = lhs;
        op->op.rhs = rhs;

        lhs = op;
    }
}


ftNode *ftParser::getOrOperand(void)
{
    static OperatorType operatorType[] = {
        { '&', AndType },
        { '\0' }
    };

    ftNode *lhs = getAndOperand();
    if (lhs == NULL) {
        return NULL;  // error already reported
    }

    for (;;) {
        ftNode *op = getOperator(operatorType);
        if (op == NULL) {  // not AND operator?
            // let caller decide what to do with next character
            return lhs;  // return what was parsed so far
        }

        ftNode *rhs = getAndOperand();
        if (rhs == NULL) {
            delete op;
            delete lhs;
            return NULL;  // error already reported
        }

        op->op.lhs = lhs;
        op->op.rhs = rhs;

        lhs = op;
    }
}


ftNode *ftParser::getAndOperand(void)
{
    ftNode *node;

    m_reader->skipWhiteSpace();

    char c = m_reader->getCurrentCharacter();
    if (c == '(') {
        m_reader->advancePosition();
        m_level++;
        node = getExpression();
        m_level--;
        if (!node) {
            return NULL;  // error already reported
        }
        return node;
    }

    return getOperand();
}


ftNode *ftParser::getOperator(OperatorType *operatorType)
{
    m_reader->skipWhiteSpace();

    char c = m_reader->getCurrentCharacter();
    while (operatorType->op) {
        if (c == operatorType->op) {
            m_reader->advancePosition();
            return new ftNode(operatorType->type);
        }
        ++operatorType;
    }
    // this is not necessarily an error (caller will determine)
    return NULL;
}


ftNode *ftParser::getOperand(void)
{
    float value;
    bool ok;
    ftNode *node;

    m_reader->skipWhiteSpace();

    if (m_reader->getFloatValue(&value, &ok)) {
        if (!ok) {
            reportError("bad number");
            return NULL;
        }

        node = new ftNode(ConstType);
        node->constant.value = value;
        return node;
    }

    std::string name = getName();
    if (!name.empty()) {
        node = new ftNode(NameType);
        node->text = name;
        node->name.index = -1;
        return node;
    }

    reportError("expected name or number");
    return NULL;  
}


// returns empty string if not valid name
std::string ftParser::getName(void)
{
    std::string name;
    char c;

    if (isalpha(c = m_reader->getCurrentCharacter())) {
        do {
            name.append(1, c);
            m_reader->advancePosition();
            c = m_reader->getCurrentCharacter();
        } while (isalnum(c) || c == '_');
    }
    return name;
}


void ftParser::reportError(const char *message, int positionOffset)
{
    std::cerr << "Input Error: " << message << " at "
        << m_reader->getLineNumber() << ":"
        << m_reader->getPosition() + positionOffset << std::endl;
}


// accessor functions for testing (should be removed)

const ftNamePool &ftParser::systemNames() const
{
    return *m_systemNames;
}

const std::vector<ftNode *> &ftParser::systemNodes() const
{
    return m_systemNodes;
}

const ftNamePool &ftParser::CompNames() const
{
    return *m_CompNames;
}


const std::string ftParser::compName(int index) const
{
    return m_CompNames->getName(index);
}

