//This is for initializing user inputs

#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
// GLOBALS SCOPED IN SOURCE FILE
double len_rad;
double len_tor;
double Nturns;
int DEBUG;
char* focus_output;
char* multi_output;
char* sfil_output;
char* mfil_output;
double* alp;
int Nradfil;
int Ntorfil;
int Nseg;
int Nthreads;
int case_alpha;
int NFalpha;
double alp_const;
int niter;

// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){

   //focus_output = "./inputfiles/focus_hsx.m12_07_12864.nc";  
   //focus_output = "./inputfiles/focus_circular_01.nc";
   focus_output = "./inputfiles/focus_hsx.m12_07.nc";
   multi_output = "./outputfiles/output_hsx.m12_07_test.nc";
   //multi_output = "./outputfiles/output_hsx.m12_07_53.nc";
  
   mfil_output = "./outputfiles/output_hsx.m12_07_test";
   //mfil_output = "./outputfiles/coils.hsx.m12_07_53_mfil";

//   wid = 0.120;
//   len = 0.060;
   case_alpha = 0;
   NFalpha = 2;
   alp_const = 0.000;
   Nseg = 128;
   len_rad = 0.060; //radial build
   len_tor = 0.030; //toroidal build
   niter = 0;   

   Nthreads = 16;
   DEBUG = 0; 
  
   Nradfil = 4;
   Ntorfil = 3;
 
}


//TODO: Implement F90-like namelist 

/*
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
