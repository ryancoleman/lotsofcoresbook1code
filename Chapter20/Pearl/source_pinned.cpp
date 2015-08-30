#include <stdlib.h>
#include <stdio.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include <string.h>
#include <iostream>
#include <vector>
#include <map>
#include <intel-coi/source/COIProcess_source.h>
#include <intel-coi/source/COIEngine_source.h>
#include <intel-coi/source/COIPipeline_source.h>
#include <intel-coi/source/COIBuffer_source.h>
#include <intel-coi/source/COIEvent_source.h>
#include "common.h"

int main()
{
    const char* SINK_NAME = "sink_mic";
    const size_t KB = 1024;
    const size_t MB = 1024 * 1024;
    size_t LOOPS = 3;
    COIRESULT result = COI_ERROR;
    uint32_t engineCount = 0;
    result = COIEngineGetCount(COI_ISA_MIC, &engineCount);
    if (result != COI_SUCCESS) {
        TSLOG_ERROR(COIResultGetName(result));
        return -1;
    }
    TSLOG_INFO("Engine count: " << engineCount);

    if (engineCount < 1) {
        TSLOG_ERROR("No engine found.");
        return -1;
    }

    COIENGINE engine = NULL;
    result = COIEngineGetHandle(COI_ISA_MIC, 0, &engine);
    if (result != COI_SUCCESS) {
        TSLOG_ERROR(COIResultGetName(result));
        return -1;
    }
    TSLOG_DEBUG("Got engine handle");

    COI_ENGINE_INFO info;
    result = COIEngineGetInfo(engine, sizeof(info), &info);
    if (result != COI_SUCCESS) {
        TSLOG_ERROR(COIResultGetName(result));
        return -1;
    }
    TSLOG_INFO("Xeon Phi info: " << info.NumCores << " cores, " << info.PhysicalMemory/MB << "MB memory");

    COIPROCESS proc = NULL;
    result = COIProcessCreateFromFile(
                                      engine,
                                      SINK_NAME,
                                      0,
                                      NULL,
                                      false,
                                      NULL,
                                      true,
                                      NULL,
                                      0,
                                      NULL,
                                      &proc
                                      );
    if (result != COI_SUCCESS)
        {
        TSLOG_ERROR("COIProcessCreateFromFile result " << COIResultGetName(result));
        return -1;
    }

    double tstart = 0;
    double ttime = 0;
    const char* funcName = "receiveData";
    COIFUNCTION func[1];
    COIPIPELINE pipeline;
    CHECK_COI_RESULT(COIPipelineCreate(proc, NULL, 0, &pipeline));
    TSLOG_INFO("Created pipeline");
    CHECK_COI_RESULT(COIProcessGetFunctionHandles(proc,
                                                  1,
                                                  &funcName,
                                                  func
                                                  ));
    TSLOG_INFO("Got handle to sink function: " << funcName);

    const char misc_data[] = "Hello COI";
    int strlength = strlen(misc_data) + 1;
    char* return_value = (char*) malloc(strlength);
    if (return_value == NULL) {
        TSLOG_ERROR("failed to allocate return value");
        return -1;
    }

    std::vector<int> bufferSize;
    for (size_t i = 1; i <= 1024; i = i * 2) {
        bufferSize.push_back(i);
    }

    {
        tstart = dtime();
        COIEVENT completion_event;

        CHECK_COI_RESULT(COIPipelineRunFunction(
                         pipeline, func[0],
                         0, NULL, NULL,
                         0, NULL,
                         misc_data, strlength,
                         return_value, strlength,
                         &completion_event
                         ));
        CHECK_COI_RESULT(COIEventWait(
                         1,
                         &completion_event,
                         -1,
                         true,
                         NULL, NULL
                         ));

        ttime = dtime() - tstart;
        TSLOG_INFO("Run function time: " << ttime*1000 << " ms");
    }
    typedef std::map<std::string, size_t> UnitMap;
    UnitMap unitMap;
    unitMap["KB"] = 1024;
    unitMap["MB"] = 1024 * 1024;
    const char cookie[] = "MAGIC";
    const char cookieEnd[] = "END";

    {
        for (size_t iter = 0; iter < LOOPS; ++iter) {
            TSLOG_INFO("\n>>> Buffer Type: COI_BUFFER_PINNED" << ": " << iter);
            for (UnitMap::const_iterator iter = unitMap.begin(); iter != unitMap.end(); ++iter) {
                float toMB = iter->second / (1024.0 * 1024.0);
                for (size_t i = 0; i < bufferSize.size(); ++i) {
                    size_t bytes = bufferSize[i] * iter->second;
                    COIBUFFER buffer;
                    char* data;
                    posix_memalign((void**) &data, 64, bytes);
                    memcpy(data, cookie, sizeof(cookie));
                    memcpy(data + bytes - sizeof(cookieEnd), cookieEnd, sizeof(cookieEnd));
                    tstart = dtime();
                    CHECK_COI_RESULT(COIBufferCreate(bytes, COI_BUFFER_PINNED, 0, data, 1, &proc, &buffer));
                    ttime = dtime() - tstart;
                    TSLOG_INFO("buffer create: " << bufferSize[i] << " " << iter->first << ", " << bufferSize[i]*toMB/(ttime) << " MB/s, " << ttime*1000 << " ms");

                    COIEVENT completion_event;
                    tstart = dtime();
                    CHECK_COI_RESULT(COIBufferWrite(buffer,0,data,bytes,COI_COPY_USE_DMA,0,NULL,&completion_event));
                    CHECK_COI_RESULT(COIEventWait(1,&completion_event,-1,true,NULL,NULL));
                    ttime = dtime() - tstart;
                    TSLOG_INFO("buffer write: " << bufferSize[i] << " " << iter->first << ", " << bufferSize[i]*toMB/(ttime) << " MB/s, " << ttime*1000 << " ms");

                    COI_ACCESS_FLAGS flag = COI_SINK_READ;
                    tstart = dtime();
                    CHECK_COI_RESULT(COIPipelineRunFunction(
                                     pipeline, func[0],
                                     1, &buffer, &flag,
                                     0, NULL,
                                     misc_data, strlength,
                                     return_value, strlength,
                                     &completion_event
                                     ));
                    CHECK_COI_RESULT(COIEventWait(
                                     1,
                                     &completion_event,
                                     -1,
                                     true,
                                     NULL, NULL
                                     ));
                    ttime = dtime() - tstart;
                    TSLOG_INFO("transferred: " << bufferSize[i] << " " << iter->first << ", " << bufferSize[i]*toMB/(ttime) << " MB/s, " << ttime*1000 << " ms");
                    COIBufferDestroy(buffer);
                    free(data);
                }
            }
        }
    }

    CHECK_COI_RESULT(COIPipelineDestroy(pipeline));
    printf("Destroyed pipeline\n");

    int8_t sink_return;
    uint32_t exit_reason;
    result = COIProcessDestroy(
                               proc,           // Process handle to be destroyed
                               -1,             // Wait indefinitely until main() (on sink side) returns
                               false,          // Don't force to exit. Let it finish executing
                                               // functions enqueued and exit gracefully
                               &sink_return,   // Don't care about the exit result.
                               &exit_reason
                               );

    if (result != COI_SUCCESS)
        {
        printf("COIProcessDestroy result %s\n", COIResultGetName(result));
        return -1;
    }

    printf("Sink process returned %d\n", sink_return);
    printf("Sink exit reason ");
    switch (exit_reason)
    {
    case 0:
        printf("SHUTDOWN OK\n");
        break;
    default:
        #ifndef _WIN32
        printf("Exit reason %d - %s\n",
               exit_reason,
               strsignal(exit_reason));
        #else
                printf("Exit reason %d\n",
                    exit_reason);

        #endif
    }

    return ((exit_reason == 0) ? 0 : -1);
}
