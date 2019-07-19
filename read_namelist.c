//This is for initializing user inputs

#include "read_namelist.h"
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
// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){
  
   //focus_output = "./inputfiles/focus_circular_01.nc";
   focus_output = "./inputfiles/focus_hsx.m12_07_128.nc";
 
//   wid = 0.120;
//   len = 0.060;
   Nseg = 128;
   wid = 0.120;
   len = 0.060;
   Nturns = 0;
   Nthreads = 4;
   DEBUG = 1;
 //  alp = (double *) malloc(Ncoils*Nseg*sizeof(double));  //set to 0 for now
   

   Nradfil = 3;
   Ntorfil = 2; 
}


//TODO: Implement F90-like namelist 

/*
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
