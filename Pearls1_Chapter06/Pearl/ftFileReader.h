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
#ifndef FTFILEREADER_H
#define FTFILEREADER_H

#include <string>
#include <fstream>

class ftFileReader
{
    std::string m_line;
    std::ifstream *m_file;
    int m_num;
    int m_pos;

    ftFileReader(void) {}
    ftFileReader(const ftFileReader &) {}
    const ftFileReader &operator =(const ftFileReader &) { return *this; }

public:
    explicit ftFileReader(const char *filnam);
    ~ftFileReader(void); 

    bool isFileOpen(void);
    char getCurrentCharacter(void);
    void advancePosition(void);
    char getCurrentCharacterAndAdvance(void);
    void skipWhiteSpace(void);
    bool getFloatValue(float *value, bool *ok);

    int getLineNumber(void);
    int getPosition(void);
};

#endif // FTFILEREADER_H
