//This is for initializing user inputs

#include "read_namelist.h"

// GLOBALS SCOPED IN SOURCE FILE
double hwid;
double hlen;
double Nturns;
char* focus_output;

// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){
   hwid = 0.060;
   hlen = 0.030;
   Nturns = 0;
   focus_output = "./example/focus_hsx.m12_07.nc";
}


//TODO: Implement F90-like namelist 

/*
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
