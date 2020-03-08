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

   focus_output = "./inputfiles/focus_hsx.m12_07_51264.nc";
   //focus_output = "./inputfiles/focus_wista_08_51264.nc";   
   //focus_output = "./inputfiles/focus_ellipse.nc";  

   multi_output = "./runs/hsx/multi_hsx.m12_07_1_1_nocomp.nc";
   //multi_output = "./runs/wista/multi_wista_08_01_04.nc";
   //multi_output = "./outputfiles/ellipse/output_ellipse_00.nc";
   
   mfil_output = "./runs/hsx/coils.hsx.m12_07_1_10_nocomp";
   //mfil_output = "./runs/wista/wista_08_01_04/coils.wista_08_01_01_04";
   //mfil_output = "./outputfiles/ellipse/output_ellipse_00";

   niter = 1;
   case_alpha = 0;
   NFalpha = 10;
   alp_const = 0.000;
   Nseg = 128;

   weight_comp = 0.000010; //complexity weighting
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

