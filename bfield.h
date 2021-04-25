#ifndef _BFIELD_H
#define _BFIELD_H
#include "globals.h"

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

//Globals 

double* xsurf; double* ysurf; double* zsurf; double* nsurfx; double* nsurfy; double* nsurfz;

int Ncoil; int iCoil; int Nzeta; int Nteta; int Nfp; double* currents; int Nseg;

double* sfilx; double* sfily; double* sfilz; int Nradfil; int Ntorfil; int Nfils;

double* mfilx; double* mfily; double* mfilz; int Ns;

//Functions

void CalculateSingleField(double x, double y, double z, \
                          double* Bx, double* By, double* Bz); 

void CalculateFieldAtPoint(double x, double y, double z, \
                    double* Bx, double* By, double* Bz); 

void CalculateFieldSerial(void);

//void CalculateFieldGPU(void);

double cosnfp(int ip);  
double sinnfp(int ip); 

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

#endif
 
