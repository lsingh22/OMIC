#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <netcdf.h>

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// MPI PARAMETERS
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

extern int pn;
extern int nproc;
extern int* startind;
extern int* endind;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// INPUTS, FILES, SETTINGS 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

extern char* focus_output;  // The input single filament focus file
extern char* multi_output; 
extern char* omic_input; 
extern char* mfil_output; // The multfilament coils file
extern char* sfil_output; // The single filament coils file

extern int case_alpha; // Determines where to get initial alpha amps ( 0 =all 0,  1 =all const, 2 =file)
extern int case_objfun;
extern int case_optimize; 

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// PLASMA EQUILIBRIUM DATA 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

extern int Nfp;
extern int isSym;
extern int isStellSym;
extern int Nteta;
extern int Nzeta;
extern size_t size_surf;
extern int size_fp;
extern int Ns;

extern double* xsurf;
extern double* ysurf;
extern double* zsurf;

extern double* nsurfx;
extern double* nsurfy;
extern double* nsurfz;
extern double* nsurfn; 
extern double surface_area;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// SINGLE FILAMENT COIL DATA 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

extern int Ncoil;
extern int iCoil;
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

// Single filament centroids

extern double* cx;  
extern double* cy;
extern double* cz;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// MULTI-FILAMENT COIL DATA 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

// Winding pack dimensions and packing data

extern int Nfils;
extern int Nradfil;
extern int Ntorfil;
extern double len_rad; 
extern double len_tor; 

// Local basis data found from Gram-Schmidt and Frenet-Serret formulae

extern double* tx;
extern double* ty;
extern double* tz;

extern double* nx;
extern double* ny;
extern double* nz;

extern double* bx;
extern double* by;
extern double* bz;

// Rotation functions and degrees of freedom

extern int NFalpha;
extern double* malp; 
extern double  alp_const;
extern double* alpamps; // Stores the degrees of freedom for rotation function
extern double* alp; // Stores the rotation function values for all coils for each segment
extern int size_alpamp;

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

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// OPTIMIZATION PARAMETERS
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

// Gradient information

extern double* derivs;
extern double* descent_dir;

// Complexity optimization options

extern double weight_comp;
extern double nvals_scaling;

// Optimization settings

extern int niter;
extern double deriv;

// GPU related
extern int isGPU;

// OpenMP related
extern int Nthreads;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
// BENCHMARKING VARIABLES 
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

// OMIC single filament magnetics

extern double* Bsfilx;
extern double* Bsfily; 
extern double* Bsfilz;
extern double* Bsfiln;
extern double* Bsfil;


// FOCUS single filament magnetics

extern double* fbn;
extern double* fbx;
extern double* fby;
extern double* fbz;

#endif
