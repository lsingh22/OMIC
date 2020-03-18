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
int isVaryRotation;
int Nthreads;
int case_alpha;
int case_opt;
int NFalpha;
double alp_const;
int niter;
double surface_area;
// UPDATES GLOBALS TO USER INPUT
void SetInputs(void){

   focus_output = "./inputfiles/focus_hsx.m12_07_debug.nc";
   //focus_output = "./inputfiles/focus_wista_08_51264.nc";   
   //focus_output = "./inputfiles/focus_ellipse.nc";  

   multi_output = "./runs/hsx/multi_hsx.m12_07_debug.nc";
   //multi_output = "./runs/wista/multi_wista_08_25.nc";
   //multi_output = "./runs/ell/multi_ell_00.nc";
   
   mfil_output = "./runs/hsx/coils.hsx.m12_07_debug";
   //mfil_output = "./runs/wista/coils.wista_08_25";
   //mfil_output = "./runs/ell/coils.ell_00";

   niter = 5;
   case_alpha = 0;
   NFalpha = 1;
   alp_const = 0.000;
   Nseg = 16;

   weight_comp = 0.015; //complexity weighting
   case_opt = 0; //0 for fbn , 1 for both fbn and fc
   nvals_scaling = 2; // the beta in the complexity formulation   

   Nradfil = 7;
   Ntorfil = 2;

   len_rad = 0.120; // 1/2 hsx 
   len_tor = 0.060;
   //len_rad = 0.120; // wista
   //len_tor = 0.060;
   //len_rad = 0.240; // ellipse
   //len_tor = 0.120;

   Nthreads = 32;
//   Nturns = 0;  not currently supported
//   DEBUG = 0; 
}
