#ifndef _MULTI_FIL_H
#define _MULTI_FIL_H
#include "globals.h"
#include "bfield.h"

void CalculateBuildDirections(void);
void WriteBuildDirections(void);
void CalculateMultiFilaments(void);
void MultiFilField(void);
void MultiFilFieldSym(void);
void MultiFilFieldmpi(void);
void WriteMultiB(void);
void WriteMultiFilaments(void);
void CalculateFiniteBuild(void);
void WriteFiniteBuild(void);
void GatherFieldData(void);

#endif
 
