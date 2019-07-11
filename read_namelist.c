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
// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){
   wid = 0.120;
   len = 0.060;
   Nturns = 0;
   focus_output = "./example/focus_hsx.m12_07.nc";
   DEBUG = 0;
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
