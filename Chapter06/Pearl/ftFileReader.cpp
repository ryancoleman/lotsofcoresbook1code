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
#include "ftFileReader.h"
#include <cctype>
#include <cerrno>
#include <cstdlib>


ftFileReader::ftFileReader(const char *filnam):
    m_num(0),
    m_pos(1)
{
    m_file = new std::ifstream(filnam, std::ios_base::in);
}


ftFileReader::~ftFileReader(void)
{
    delete m_file;
}


bool ftFileReader::isFileOpen()
{
    return m_file->is_open();
}


char ftFileReader::getCurrentCharacter(void)
{
    while (m_pos >= (int)m_line.length()) {
        if (m_pos == (int)m_line.length()) {
            return '\n';
        }
        if (m_pos > (int)m_line.length()) {
            // read new line
            if (!std::getline(*m_file, m_line, '\n')) {
                return -1;  // end-of-file
            }
            m_pos = 0;
            m_num++;
        }
    }
    return m_line[m_pos];
}


void ftFileReader::advancePosition(void)
{
    m_pos++;
}


char ftFileReader::getCurrentCharacterAndAdvance(void)
{
    char c = getCurrentCharacter();
    advancePosition();
    return c;
}


void ftFileReader::skipWhiteSpace(void)
{
    while (isspace(getCurrentCharacter())) {
        getCurrentCharacter();
        advancePosition();
    }
}


// returns: true and ok == true if good number
// returns: true and ok == false if bad number
// returns: false if no number
bool ftFileReader::getFloatValue(float *value, bool *ok)
{
    char c = getCurrentCharacter();
    if (c != '-' && !isdigit(c)) {
        return false;
    }

    const char *beginPtr = m_line.c_str() + m_pos;
    char *endPtr;
    errno = 0;
    *value = strtof(beginPtr, &endPtr);
    *ok = errno != ERANGE;

    m_pos += endPtr - beginPtr;  // advance past number
    return true;
}


int ftFileReader::getLineNumber(void)
{
    return m_num;
}


int ftFileReader::getPosition(void)
{
    return m_pos;
}

