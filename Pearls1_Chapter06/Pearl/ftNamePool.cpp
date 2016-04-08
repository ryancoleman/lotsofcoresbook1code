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
#include "ftNamePool.h"

ftNamePool::ftNamePool(void) 
{

}

int ftNamePool::addName(const std::string &name)
{
    bool exists;
    return addName(name, &exists);
}


// adds a name to the pool
int ftNamePool::addName(const std::string &name, bool *exists)
{
    std::pair<nameItr, bool> ret;
    int index;
    bool used;

    if (deleted.size() > 0)
    {
        index = deleted.back();  // use last deleted index
        used = true;
    }
    else
    {
        index = namemap.size();  // new index at end
        used = false;
    }
    ret = namemap.insert(std::pair<std::string, int>(name, index));
    if ((*exists = ret.second == false))
    {
        // item already exists, return index
        return (*ret.first).second;
    }

    if (used)
    {
        deleted.pop_back();  // remove from deleted list
    }
    else
    {
        entries.resize(index + 1);  // increase size
    }
    entries[index] = ret.first;
    return index;
}


// deletes a name from the pool
int ftNamePool::deleteName(const std::string &name)
{
    nameItr it = namemap.find(name);
    if (it == namemap.end())
    {
        return NotFound;
    }
    int index = (*it).second;
    namemap.erase(it);
    entries[index] = namemap.end();
    deleted.push_back(index);
    return index;
}


// returns the name for an index
const char *ftNamePool::getName(int index) const
{
    if (index < 0 || index >= (int)entries.size())
    {
        return NULL;
    }

    nameItr entry = entries[index];
    if (entry == namemap.end())
    {
        return NULL;
    }
    return (*entries[index]).first.c_str();
}


// returns the index of a name
int ftNamePool::getIndex(const std::string &name) const
{
    constNameItr it = namemap.find(name);
    if (it == namemap.end())
    {
        return NotFound;
    }
    return (*it).second;
}


// returns the number of names in the name pool
int ftNamePool::count(void) const
{
    return namemap.size();
}


// returns the size of the name pool (including deleted entries)
int ftNamePool::size(void) const
{
    return entries.size();
}

