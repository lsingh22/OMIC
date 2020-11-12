#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
// GLOBALS SCOPED IN SOURCE FILE
double len_rad; double len_tor; double Nturns;
char* focus_output; char* multi_output; char* sfil_output; char* mfil_output;

double* alp; int Nradfil; int Ntorfil; int Nseg; int isVaryRotation; int case_alpha; int case_objfun; double nvals_scaling;
int NFalpha; double alp_const; int niter; double surface_area; 

int nproc; int* startind; int* endind;
int isStellSym; int iCoil; int size_alpamp; int size_fp; int Nfils; int Ns;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
#define MAXCHAR 32

void ReadInputs(void){
//----------------------------------------------------------------------------------------------------
// Read in user-specified inputs
//----------------------------------------------------------------------------------------------------
 
   FILE *fp;
   char str[MAXCHAR];
   char* filename = omic_input;

   const char s[2] = "=";
   const char t[2] = ";";
   char *token;
   
   fp = fopen(filename, "r");
   
   if (fp == NULL)
   {
      printf("Could not open file %s",filename);
      return;
   }

   register int i; i=0;
   
   while (fgets(str, MAXCHAR, fp) != NULL)
   {
      i++;
      token = strtok(str, " ");
      token = strtok(NULL, " "); 
   
      if(i==1)
      {
         sscanf(token, "%d", &niter);         
      }
      if(i==2)
      {
         sscanf(token, "%d", &NFalpha);         
      } 
      if(i==3)
      {
         sscanf(token, "%d", &isStellSym);         
      }
      if(i==4)
      {
         sscanf(token, "%d", &case_alpha);         
      }
      if(i==5)
      {
         sscanf(token, "%lf", &alp_const);         
      }
      if(i==6)
      {
         sscanf(token, "%lf", &weight_comp);         
      }
      if(i==7)
      {
         sscanf(token, "%d", &case_objfun);         
      }
      if(i==8)
      {
         sscanf(token, "%lf", &nvals_scaling);         
      }
      if(i==9)
      {
         sscanf(token, "%d", &Nseg);         
      }
      if(i==10)
      {
         sscanf(token, "%d", &Nradfil);         
      }
      if(i==11)
      {
         sscanf(token, "%d", &Ntorfil);         
      }
      if(i==12)
      {
         sscanf(token, "%lf", &len_rad);         
      }
      if(i==13)
      {
         sscanf(token, "%lf", &len_tor);         
      }
   }
   fclose(fp);
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

