#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <netcdf.h>

///// INPUTS AND SETTINGS /////

extern char* focus_output;  // The input single filament focus file
extern char* alpha_input; // TODO: This will be in the namelist eventually
extern char* namelist; // TODO: The input namelist file

extern int Nthreads; // Threads to use in OpenMP field calculation
extern int DEBUG; // Option to suppress output
extern int case_alpha; // Determines where to get initial alpha amps ( 0 =all 0,  1 =all const, 2 =file)

///// PLASMA EQUILIBRIUM DATA  /////

extern int Nfp;
extern int isSym;
extern int Nteta;
extern int Nzeta;
extern size_t size_surf;

extern double* xsurf;
extern double* ysurf;
extern double* zsurf;

extern double* nsurfx;
extern double* nsurfy;
extern double* nsurfz;
extern double* nsurfn; //Jacocbian

///// SINGLE FILAMENT COIL DATA /////

extern int Ncoils;
extern int NFcoil;
extern double* coilspace;
extern size_t size_coilspace;
extern double* coilamps;
extern double* currents;
extern int Nseg;
int* ind_arr;

// Single filament XYZ points

extern double* sfilx;
extern double* sfily; 
extern double* sfilz;

// MULTI single filament magnetics

extern double* Bsfilx;
extern double* Bsfily; 
extern double* Bsfilz;
extern double* Bsfiln;
extern double* Bsfil;

// Single filament centroids

extern double* cx;  
extern double* cy;
extern double* cz;

///// MULTI-FILAMENT COIL DATA /////

// Winding pack dimensions and packing data

extern int Nradfil;
extern int Ntorfil;
extern double wid; 
extern double len; 

// Local basis data found from Gram-Schmidt and Frenet-Serret formulae

// Frenet unit tangent vectors

extern double* tx;
extern double* ty;
extern double* tz;

// Gram-Schmidt normal vectors

extern double* nx;
extern double* ny;
extern double* nz;

// Binormal vectors from cross product of unit tangent and normal vectors

extern double* bx;
extern double* by;
extern double* bz;

// Alpha rotation parameters

extern int NFalpha;
extern double* alpampsinit; // Stores the namelist alpha amplitudes
extern double  alp_const;
extern double* alpamps; // Stores the objective function alpha amplitudes
extern double* alp; // Stores the real space alpha values for calculating multifilaments

// Multi-Filament XYZ points

extern double* mfilx;
extern double* mfily;
extern double* mfilz;

// Multi-Filament magnetics

extern double* Bmfilx;
extern double* Bmfily; 
extern double* Bmfilz;
extern double* Bmfiln;
extern double* Bmfil;

// Finite build coil XYZ points

extern double* finx;
extern double* finy;
extern double* finz;

// Derivatives (only field rn)

extern double* derivs;
extern double* descent_dir;


///// BENCHMARKING VARIABLES /////

// FOCUS single filament magnetics

extern double* fbn;
extern double* fbx;
extern double* fby;
extern double* fbz;

#endif
