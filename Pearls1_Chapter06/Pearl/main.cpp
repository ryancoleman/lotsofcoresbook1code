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
using std::cerr;
using std::cout;
using std::endl;

#include <random>
using std::random_device;

#include <omp.h>

#include "ftEvaluator.h"

namespace ispc
{
    extern "C" void* ispc_alloc(size_t nbytes);
}

void printUsage()
{
    cerr << "Usage: fte <input file> <top-level system> <num trees>" << endl;
}

int main(int argc, const char* argv[])
{
    if (argc < 4 || argc > 4) {
        printUsage();
        return 0;
    }

    // Parse inputs and fault-tree file ///////////////////////////////////////

    // Create the fault tree evaluator
    ftEvaluator evaluator;

    // Open the fault-tree file from the command line
    if (!evaluator.openFile(argv[1])) {
        return 1;
    }

    // Get the number of evaluations from the command line
    size_t nEvals = atoi(argv[3]);

    // Set the system that will be evaluated
    if (!evaluator.setSystemToEvaluate(argv[2])) {
        cerr << "Failed to use system '" << argv[2] << "' to evaluate!" << endl;
        return 1;
    }

    // Get the number of components that were parsed out of the fault tree file
    size_t nComps = evaluator.getNumberOfComponents();

	cout << "    Top-level system: '" << argv[2] << "'" << endl;
	cout << "Number of components: " << nComps << endl;
	cout << "In-flight memory cap: " << CHUNK_SIZE << "MB" << endl;
	cout << endl;
	cout.flush();

    // Calcualte number of chunks to be executed //////////////////////////////

    // Chunk size (in bytes)
    size_t cSize   = CHUNK_SIZE * 1024ul * 1024ul;

    // Number of components per chunk
    size_t cChunk  = cSize / sizeof(float);

    // Total size of component buffer (in bytes)
    size_t tBytes  = nComps * nEvals * sizeof(float);

    // Total number of chunks
    size_t nChunks = tBytes / cSize;

    // Number of evaluations per chunk
    size_t cEvals  = cChunk / nComps;

    // Buffer of component states that will be used for evaluations
    float *compVals = nullptr;

    // Create a results buffer and status buffer
    float    *results = (float*)ispc::ispc_alloc(nEvals*sizeof(float));
    ftStatus *status  = (ftStatus*)ispc::ispc_alloc(nEvals*sizeof(ftStatus));

    // Check to see if we only need one chunk
    if (nChunks == 0) {
        // Allocate what we need and set nChunks to 1
        nChunks  = 1;
        compVals = (float*)ispc::ispc_alloc(tBytes);
        results  = (float*)ispc::ispc_alloc(nEvals*sizeof(float));
        status   = (ftStatus*)ispc::ispc_alloc(nEvals*sizeof(ftStatus));

        // Set number of evaluations to the number of evaluations we have
        evaluator.setNumberOfEvaluations(nEvals);
    } else {
        // Allocate one chunk's worth of bytes
        compVals = (float*)ispc::ispc_alloc(cSize);
        results  = (float*)ispc::ispc_alloc(cEvals*sizeof(float));
        status   = (ftStatus*)ispc::ispc_alloc(cEvals*sizeof(ftStatus));

        // Set number of evaluations to the number of evaluations per chunk
        evaluator.setNumberOfEvaluations(cEvals);
    }

    double start, end;

    // Track the time it takes to evaluate the trees
    start = omp_get_wtime();

    // Do the evaluations
    cout << "Doing " << nEvals << " evaluations with system '" << argv[2]
         << "'..." << endl;
    cout.flush();

    for (size_t i = 0; i < nChunks; ++i) {

#ifdef INIT_COMPONENTS
        cout << "Initializing random component values for " << nComps << "x"
             << nEvals << " components..." << endl;

        start = omp_get_wtime();

        // Create random number generator for randomizing component states
        random_device rd;
        std::uniform_real_distribution<float> rng(0.f, 1.f);

        // Initialize the buffer with random component states
        #pragma omp parallel for
        for (size_t i = 0; i < nComps*nEvals; ++i)
            compVals[i] = rng(rd);

        // Calculate the runtime cost of generating random numbers
        end = omp_get_wtime();
        cout << "...finished in " << end - start << "s" << endl;
        cout.flush();
#endif

        // Set the buffer in the evaluator for reading component states
        evaluator.setComponentData(compVals);

        if (!evaluator.evaluate(results, status)) {
            cerr << "Evaluation failed!" << endl;

            // Cleanup memory
            delete compVals;
            delete results;
            delete status;

            return 1;
        }
    }

    // Calculate the runtime cost of evaluating the trees
    end = omp_get_wtime();
    cout << "...finished everything in " << end - start << "s" << endl;
    cout.flush();

#ifdef INIT_COMPONENTS
    // Average all of the results together and output to the console
    float final = 0;

    for (size_t i = 0; i < nEvals; ++i) {
        if (status[i] != GoodStatus) {
            cerr << "Bad evaluation at: " << i << endl;
            continue;
        }
        final += results[i];
    }

    final /= static_cast<float>(nEvals);

    cout << "Average pk: " << final << endl;
#endif

    // Cleanup memory
    delete compVals;
    delete results;
    delete status;

    return 0;
}
