#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>
#include <stdio.h>
//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
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
int isVaryRotation;
int Nthreads;
int case_alpha;
int case_opt;
int NFalpha;
double alp_const;
int niter;
double surface_area;
int nproc;
int* startind;
int* endind;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void SetInputs(void){
//----------------------------------------------------------------------------------------------------
// Write inputs to globals
//----------------------------------------------------------------------------------------------------

   focus_output = "./inputfiles/focus_hsx.m12_07_51264.nc";
   //focus_output = "./inputfiles/focus_wista_53.nc";   
   //focus_output = "./inputfiles/focus_ellipse.nc";  
   //focus_output = "./inputfiles/focus_hsx.m12_07_debug.nc";
   
   multi_output = "./runs/hsx/multi_hsx.m12_07_mpi3.nc";
   //multi_output = "./runs/wista_ss/multi_wista_ss_00.nc";
   //multi_output = "./runs/ell/multi_ell_01.nc";
   //multi_output = "./runs/multi_debug2.nc";
   
   mfil_output = "./runs/hsx/coils.hsx.m12_07_mpi3";
   //mfil_output = "./runs/wista_ss/coils.wista_ss_00";
   //mfil_output = "./runs/ell/coils.ell_01";
   //mfil_output = "./runs/coils.test";

   niter = 0;
   case_alpha = 0;
   NFalpha = 5;
   alp_const = 0.000;
   Nseg = 128;

   weight_comp = 0.01; //complexity weighting
   case_opt = 1; //0 for fbn , 1 for both fbn and fc
   nvals_scaling = 2; // the beta in the complexity formulation   

   Nradfil = 2;
   Ntorfil = 2;

   len_rad = 0.120; // 1/2 hsx 
   len_tor = 0.060;
   //len_rad = 0.1500; // wista
   //len_tor = 0.0750;
   //len_rad = 0.3500; // ellipse
   //len_tor = 0.7000;

   Nthreads = 32;
//   DEBUG = 0; 
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

void MPInit(void){
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

