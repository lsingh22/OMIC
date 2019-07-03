//This is for initializing user inputs

#include "read_namelist.h"

double hwid;
double hlen;
double Nturns;
char* focus_output;

//This function updates the globals in globals.h
void SetInputs(void){
   hwid = 0.060;
   hlen = 0.030;
   Nturns = 0;
   focus_output = "/home/luquants/multi/example/focus_hsx.m12_07.nc";
}


//TODO: Implement F90-like namelist 

/*
//This function reads the input file and updates the globals in globals.h
void SetDefault(char *namelist){
  namelist = argv;
 }
*/
