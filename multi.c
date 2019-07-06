#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "read_namelist.h"  
#include "read_focus.h"

//THIS IS THE MAIN FOR THE MULTIFILAMENT OPTIMIZATION CODE
int main(int argc, char **argv) {
   
   SetInputs();
   ReadFocusInts(focus_output);
   ReadFocusArrays(focus_output);
 
//DEBUG 
/*
   // printf("%f\n", coilspace[1]);
   printf("%f\n", coilspace[1]);
   printf("%f\n", xsurf[16383]);
   printf("%f\n", ysurf[16383]);
   printf("%f\n", zsurf[16383]);
   printf("%f\n", nsurfx[16383]);
   printf("%f\n", nsurfy[16383]);
   printf("%f\n", nsurfz[16383]);
   printf("%f\n", fbn[16383]);
   printf("%f\n", fbx[16383]);
   printf("%f\n", fby[16383]);
   printf("%f\n", fbz[16383]);
    
*/
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
