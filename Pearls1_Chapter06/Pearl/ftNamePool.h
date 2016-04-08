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
#ifndef NAMEPOOL_H
#define NAMEPOOL_H

#include <string>
#include <vector>
#include <map>


typedef std::map<std::string,int> nameMap;
typedef nameMap::iterator nameItr;
typedef nameMap::const_iterator constNameItr;

class ftNamePool {
    nameMap namemap;
    std::vector<nameItr> entries;
    std::vector<int> deleted;

    ftNamePool(const ftNamePool &) {}
    const ftNamePool &operator =(const ftNamePool &) { return *this; }

public:
    ftNamePool(void);

    enum {
        NotFound = -1
    };

    int addName(const std::string &name);
    int addName(const std::string &name, bool *exists);
    int deleteName(const std::string &name);
    const char *getName(int index) const;
    int getIndex(const std::string &name) const;
    int count(void) const;
    int size(void) const;
};


#endif // NAMEPOOL_H
