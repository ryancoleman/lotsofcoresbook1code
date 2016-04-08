
#include <stdlib.h>
#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <intel-coi/sink/COIPipeline_sink.h>
#include <intel-coi/sink/COIProcess_sink.h>
#include <intel-coi/common/COIMacros_common.h>
#include <intel-coi/common/COISysInfo_common.h>
#include <intel-coi/common/COIEvent_common.h>
#include "common.h"

int main(int , char**)
{
    TSLOG_DEBUG("Hello from the sink!");

    UNUSED_ATTR COIRESULT result;

    result = COIPipelineStartExecutingRunFunctions();

    assert(result == COI_SUCCESS);

    COIProcessWaitForShutdown();
    return 0;
}

COINATIVELIBEXPORT
void receiveData(uint32_t         in_BufferCount, // The number of buffers passed to the run function.
                 void**           in_ppBufferPointers, // An array that is in_BufferCount in length that contains the sink side virtual addresses for each buffer passed in to the run function.
                 uint64_t*        in_pBufferLengths, // An array that is in_BufferCount in length of uint32_t integers describing the length of each passed in buffer in bytes.
                 void*            in_pMiscData, // Pointer to the MiscData passed in when the run function was enqueued on the source.
                 uint16_t         in_MiscDataLength, // Length in bytes of the MiscData passed in when the run function was enqueued on the source.
                 void*            in_pReturnValue, // Pointer to the location where the return value from this run function will be stored.
                 uint16_t         in_ReturnValueLength) // Length in bytes of the user-allocated ReturnValue pointer.
{
    UNREFERENCED_PARAM(in_BufferCount);
    UNREFERENCED_PARAM(in_ppBufferPointers);
    UNREFERENCED_PARAM(in_pBufferLengths);
    UNREFERENCED_PARAM(in_pMiscData);
    UNREFERENCED_PARAM(in_MiscDataLength);
    for (size_t i = 0; i < in_BufferCount; ++i) {
        TSLOG_DEBUG("sink: buffer received size " << in_pBufferLengths[i]/1024 << " KB");
    }
}
