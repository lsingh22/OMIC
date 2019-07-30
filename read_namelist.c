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
int Nthreads;
int case_alpha;
int NFalpha;
double alp_const;

// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){
  
   //focus_output = "./inputfiles/focus_circular_01.nc";
   focus_output = "./inputfiles/focus_hsx.m12_07.nc";
 
//   wid = 0.120;
//   len = 0.060;
   case_alpha = 1;
   NFalpha = 1;
   alp_const = 0.001;
   Nseg = 128;
   wid = 0.120;
   len = 0.060;
   Nturns = 0;
   Nthreads = 128;
   DEBUG = 1; 
  
   Nradfil = 2;
   Ntorfil = 2;

   //unpack_alpha(case_alpha);
 
}


//TODO: Implement F90-like namelist 

/*
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
