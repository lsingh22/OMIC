#ifndef _STARTUP_H
#define _STARTUP_H

#include "globals.h"

void OMICStartup(char* ext);
void ReadInputFile(void);
void Initialize(void);
void MPInit(void);
int CFileExists(char* file);

#endif
