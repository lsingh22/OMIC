#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <netcdf.h>

//TODO: Organize by header file reference

// Input focus file

extern char* focus_output;
extern char* namelist;
// Store global variables

extern int DEBUG;
extern int Ncoils;
extern int Nfp;
extern int isSym;
extern int NFcoil;  // extern int* NFcoil after Focus update
extern int Ncoils;
extern int Ntor;
extern int Npol;
extern int Nteta;
extern int Nzeta;
extern double* coilamps;
extern double* cx;  //TODO: Once sfils are allowed to vary, make this COM instead according to Stuart
extern double* cy;
extern double* cz;
extern double* coilspace;
extern size_t size_coilspace;
extern size_t size_surf;
extern double* currents;
extern double* alphas;
extern int Nseg;
extern double Nturns;

extern double* xsurf;
extern double* ysurf;
extern double* zsurf;

extern double* nsurfx;
extern double* nsurfy;
extern double* nsurfz;

extern double* fbn;
extern double* fbx;
extern double* fby;
extern double* fbz;

extern double* sfilx;
extern double* sfily; 
extern double* sfilz;

extern double* Bsfilx;
extern double* Bsfily; 
extern double* Bsfilz;
extern double* Bsfiln;

extern double* Bmfilx;
extern double* Bmfily; 
extern double* Bmfilz;
extern double* Bmfiln;



extern double hwid; //normal half length 0.060
extern double hlen; //binormal half length 0.030

extern double* tx;
extern double* ty;
extern double* tz;

extern double* nx;
extern double* ny;
extern double* nz;

extern double* bx;
extern double* by;
extern double* bz;

extern double* finx;
extern double* finy;
extern double* finz;

extern double* mfilx;
extern double* mfily;
extern double* mfilz;

#endif
