#ifndef _READ_FOCUS_H
#define _READ_FOCUS_H

#include "globals.h"
#include <netcdf.h>

extern void ReadFocusInts(char* output_file);
extern void ReadFocusArrays(char* output_file);
extern void WriteBoundaryNC(void);

#endif
