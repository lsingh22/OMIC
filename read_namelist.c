//This is for initializing user inputs

#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
// GLOBALS SCOPED IN SOURCE FILE
double wid;
double len;
double Nturns;
int DEBUG;
char* focus_output;
double* alp;
int Nradfil;
int Ntorfil;
int Nseg;
int isVaryRotation;
int Nthreads;
int case_alpha;
int case_opt;
int NFalpha;
double alp_const;
int niter;

// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){
  
   //focus_output = "./inputfiles/focus_circular_01.nc";
   focus_output = "./inputfiles/focus_hsx.m12_07.nc";
   multi_output = "./outputfiles/output.nc";

//   wid = 0.120;
//   len = 0.060;
   case_alpha = 0;
   case_opt = 0;
   NFalpha = 1;
   alp_const = 0.5;
   Nseg = 128;
   wid = 0.055;
   len = 0.105;
   //wid = 0.100;
   //len = 0.200;
   Nturns = 0;
   Nthreads = 16;
   DEBUG = 0; 
  
   Nradfil = 3;
   Ntorfil = 2;
   isVaryRotation = 1;
   niter = 1;
   //unpack_alpha(case_alpha);
 
}


//TODO: Implement F90-like namelist 

/*
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
