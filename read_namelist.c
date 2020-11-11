#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
// GLOBALS SCOPED IN SOURCE FILE
double len_rad; double len_tor; double Nturns;
char* focus_output; char* multi_output; char* sfil_output; char* mfil_output;

double* alp; int Nradfil; int Ntorfil; int Nseg; int isVaryRotation; int case_alpha; int case_objfun;
int NFalpha; double alp_const; int niter; double surface_area; 

int nproc; int* startind; int* endind;
int isStellSym; int iCoil; int size_alpamp; int size_fp; int Nfils; int Ns;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void SetInputs(void){
//----------------------------------------------------------------------------------------------------
// Write inputs to globals
//----------------------------------------------------------------------------------------------------

   niter = 0;
 
   case_alpha = 0;
   NFalpha = 1;
   alp_const = 0.000;
 
   isStellSym = 0;
   Nseg = 128;

   weight_comp = 0.01; //complexity weighting
   case_objfun = 1; //0 for fbn , 1 for both fbn and fc
   nvals_scaling = 2; // the beta in the complexity formulation   

   Nradfil = 7;
   Ntorfil = 2;

   len_rad = 0.120; // 1/2 hsx 
   len_tor = 0.060;
   //len_rad = 0.1500; // wista
   //len_tor = 0.0750;
   //len_rad = 0.3500; // ellipse
   //len_tor = 0.7000;
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void Initialize(void){

   if(isStellSym == 1)
   {
      iCoil = Ncoil / ( 2 * Nfp );
      size_fp = Nteta * Nzeta / ( 2 * Nfp ); 
      Ns = 1; 
   }
   else if(isStellSym == 0)
   {
      iCoil = Ncoil / Nfp;
      size_fp = Nteta * Nzeta / Nfp; 
      Ns = 0;
   }

   size_alpamp = iCoil * ( 2 * NFalpha + 1 ); 
   Nfils = Nradfil * Ntorfil; 
}

#define MAXCHAR 32

void ReadInputs(void){

   FILE *fp;
   char str[MAXCHAR];
   char* filename = "sample.input";

   const char s[2] = "=";
   const char t[2] = ";";
   char *token;
   
   fp = fopen(filename, "r");
   
   if (fp == NULL)
   {
      printf("Could not open file %s",filename);
      return;
   }
   
   char *cut = "=";
   int index;
   char *strcopy = str;

   while (fgets(str, MAXCHAR, fp) != NULL)
   {
      token = strtok(str, " ");
      token = strtok(NULL, " "); 
      //printf( "%s\n", token);
   }

   fclose(fp);
}


void MPISetup(void){
//----------------------------------------------------------------------------------------------------
// Allocate and write some useful arrays for mpi implementation
//----------------------------------------------------------------------------------------------------
 
   startind = (int*) malloc(nproc*sizeof(int));
     endind = (int*) malloc(nproc*sizeof(int)); 

   int size_fp = Nteta*Nzeta / Nfp;
   int floor_points_per_proc = (size_fp - (size_fp % 2)) / nproc; //returns even average of points per processor
   int endpn = nproc - 1; //the index value of last processor

   for(int i=0;i<nproc;i++)
   {
      *(startind+i) = 0 + i*floor_points_per_proc;
      *(endind+i)   = (i+1)*floor_points_per_proc - 1; 
   }
   
   //increment the last index if there are an odd number of points
   if(size_fp % 2 == 1)
   {
      *(endind+endpn) = *(endind+endpn) + 1;
   }
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

