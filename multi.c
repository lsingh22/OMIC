#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_namelist.h"  
#include "read_focus.h"
#include "single_fil.h"
#include "multi_fil.h"
#include "bfield.h"

//THIS IS THE MAIN FOR THE MULTIFILAMENT OPTIMIZATION CODE
int main(int argc, char **argv){
   
   SetInputs();
   ReadFocusInts(focus_output);
   ReadFocusArrays(focus_output);
   UnpackSingleFilaments();
   CalculateBuildDirections();
   CalculateMultiFilaments();
   WriteMultiFilaments(); 

   //CalcSingleFilsB();
   //printf("%f\n", mfilx[93]);
  // printf("%.15f\n", Bsfilx[1]*nsurfx[1]+Bsfily[1]*nsurfy[1]+Bsfilz[1]*nsurfz[1]);
   //printf("%f\n", fbx[128*128]); 

   // printf("%f\n", coilspace[1]);
   if ( DEBUG == 1)
   {
   
      printf("The first entry of coilspace is:   %f\n", coilspace[0]);
      printf("The last entry of xsurf is:   %f\n", xsurf[16383]);
      printf("The first entry of ysurf is:   %f\n", ysurf[16383]);
      printf("The first entry of zsurf is:   %f\n", zsurf[0]);
      printf("The first entry of nx is:   %f\n", nsurfx[0]);
      printf("The last entry of ny is:   %f\n", nsurfy[16383]);
      printf("The first entry of Bn is:   %f\n", fbn[0]);
      printf("The second entry of Bx is:   %f\n", fbx[1]);
      printf("The first x coordinate of the first coil is:   %f\n", sfilx[0]);
      printf("The current of the third coil is:   %f\n", currents[2]);
      printf("The centroid of the first coil is:   %f   %f   %f\n", cx[0],cy[0],cz[0]);
   
      //WriteSingleFilaments();


   }
/* 
   ReadInputs();
   ReadFocus();
   LoadSingleFilaments();
   LoadBoundary();
   UnpackSingleFilaments();
   CalculateBuildDirections();
   WriteMultifilaments();
   CalculateField();
   

   Debug:
   WriteCoilFiles();
   DiagnoseField();
  
   //Allocate pointers and call function to load wout file, load_wout.c
   //
   //Allocate pointers and call function to load .input file, load_input.c
   //
   //Allocate pointers and call function to load .focus file, load_focus.c
   //
   //Call function to produce XYZ points for single filament, single_fil.c
   //
   //Write python code to plot filaments, single or multi, plotFil.py
   //
   //Hardcode alphas to 0 and call function to produce XYZ points for multifilament, multi_fil.c
   //
   //write python code to plot finite build coils, plotBuild.py 
*/  
}
