#include "read_namelist.h"
#include "alpha.h"
#include <stdlib.h>

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

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----
 
void SetInputs(void){
//----------------------------------------------------------------------------------------------------
// Write inputs to globals
//----------------------------------------------------------------------------------------------------

   //focus_output = "./inputfiles/focus_hsx.m12_07_51264.nc";
   focus_output = "./inputfiles/focus_wista_53.nc";   
   //focus_output = "./inputfiles/focus_ellipse.nc";  

   //multi_output = "./runs/hsx/multi_hsx.m12_07_test.nc";
   multi_output = "./runs/wista_ss/multi_wista_ss_00.nc";
   //multi_output = "./runs/ell/multi_ell_01.nc";
   
   //mfil_output = "./runs/hsx/coils.hsx.m12_07_test";
   mfil_output = "./runs/wista_ss/coils.wista_ss_00";
   //mfil_output = "./runs/ell/coils.ell_01";

   niter = 0;
   case_alpha = 0;
   NFalpha = 5;
   alp_const = 0.000;
   Nseg = 128;

   weight_comp = 0.000; //complexity weighting
   case_opt = 1; //0 for fbn , 1 for both fbn and fc
   nvals_scaling = 2; // the beta in the complexity formulation   

   Nradfil = 7;
   Ntorfil = 2;

   //len_rad = 0.120; // 1/2 hsx 
   //len_tor = 0.060;
   len_rad = 0.1500; // wista
   len_tor = 0.0750;
   //len_rad = 0.3500; // ellipse
   //len_tor = 0.7000;

   Nthreads = 32;
//   Nturns = 0;  not currently supported
//   DEBUG = 0; 
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//---- 
