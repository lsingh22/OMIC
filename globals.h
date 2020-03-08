#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <netcdf.h>

///// INPUTS AND SETTINGS /////

extern char* focus_output;  // The input single filament focus file
extern char* multi_output; // TODO: This will be in the namelist eventually
extern char* namelist; // TODO: The input namelist file
extern char* mfil_output; // The multfilament coils file
extern char* sfil_output; // The single filament coils file


extern int Nthreads; // Threads to use in OpenMP field calculation
extern int DEBUG; // Option to suppress output
extern int case_alpha; // Determines where to get initial alpha amps ( 0 =all 0,  1 =all const, 2 =file)
extern int case_opt;
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
extern double surface_area;

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
extern double len_rad; 
extern double len_tor; 

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
extern double* malp; // Stores the namelist alpha amplitudes
extern double  alp_const;
extern double* alpamps; // Stores the objective function alpha amplitudes
extern double* alp; // Stores the real space alpha values for calculating multifilaments

// Multi-Filament XYZ points

extern double* mfilx;
extern double* mfily;
extern double* mfilz;

// Finite build surface xyz points

extern double* ffilx;
extern double* ffily;
extern double* ffilz;

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

// Complexity optimization
extern double weight_comp;
extern double nvals_scaling;

// Optimization settings
extern int niter;
extern double deriv;
extern double multi_error_init;
extern double comp_penalty_init; 

///// BENCHMARKING VARIABLES /////

// FOCUS single filament magnetics

extern double* fbn;
extern double* fbx;
extern double* fby;
extern double* fbz;

#endif
