#ifndef _BFIELD_H
#define _BFIELD_H
#include "globals.h"

void CalculateSingleField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz); 

void CalculateMultiFieldSym(double x, double y, double z, \
                          double* Bx, double* By, double* Bz); 

double cosnfp(int ip);  
double sinnfp(int ip); 

#endif
 
