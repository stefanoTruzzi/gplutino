// #ifdef _USE_NVTX// SPDX-License-Identifier: BSD-3-Clause
#ifdef _USE_NVTX
#include <nvtx3/nvToolsExt.h>
#include <string.h>
#endif

void mynvtxstart_(const char *name, int color) {
#ifdef _USE_NVTX

   int hash = 0;

   nvtxEventAttributes_t eventAttrib = {0};
   eventAttrib.version = NVTX_VERSION;
   eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
   eventAttrib.colorType = NVTX_COLOR_ARGB;
   eventAttrib.color = color;
   eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
   eventAttrib.message.ascii = name;
   nvtxRangePushEx(&eventAttrib);
#endif
}

void mynvtxstop_() {
#ifdef _USE_NVTX
    nvtxRangePop();
#endif
}


/* 
STANDARD
CC      = nvc++
CFLAGS  = -c -O3 -cudalib=nvtx3 -D_USE_NVTX
LDFLAGS = -lm -cudalib=nvtx3

o
NVIDIA
CC      = nvc++
CFLAGS  = -c -acc -std=c++17 -gpu=rdc,managed,lineinfo -Minfo=accel -Minline=1000 -O3 -D_USE_NVTX
LDFLAGS = -lm -acc -gpu=rdc,managed -cudalib=nccl -lnccl -cudalib=nvtx3


nsys profile --trace=cuda,nvtx,openacc -o NOME_REPORT -f true ./pluto

nsys-ui

*/
