#ifndef _GLOBALS_H
#define _GLOBALS_H


//TODO: Organize by header file reference

// Input focus file

extern char* focus_output;




// Store global variables

extern int Ncoils;
extern int Nfp;
extern int isSym;
extern int NFcoil;  // extern int* NFcoil after Focus update
extern int Ncoils;
extern int Ntor;
extern int Npol;
extern double* coilamps;
extern double* centroids;  //TODO: Once sfils are allowed to vary, make this COM instead according to Stuart
extern double* currents;
extern double* alphas;
extern int Nseg;
extern double Nturns;

extern double* sfilx;
extern double* sfily; 
extern double* sfilz;

extern double hwid; //normal half length 0.060
extern double hlen; //binormal half length 0.030
extern double norm;
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
