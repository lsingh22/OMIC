#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "read_namelist.h"  
#include "read_focus.h"
#include "single_fil.h"
#include "multi_fil.h"
#include "bfield.h"
#include "alpha.h"
#include "output.h"
#include "solvers.h"

//THIS IS THE MAIN FOR THE MULTIFILAMENT OPTIMIZATION CODE
int main(int argc, char **argv){
   
   double tot_time;
   clock_t start, end;
   // At some point use higher resolution timer 
   
   start = clock();

  
 
   //printf("This is 1 \n");
   SetInputs();
   //printf("This is 2 \n");
   ReadFocusInts(focus_output);
   //printf("This is 3 \n");
   ReadFocusArrays(focus_output);
   //printf("This is 4 \n");
   UnpackSingleFilaments();
   Init_alpha(case_alpha);
   //printf("This is 6 \n");
   CalculateMultiFilaments();

   printf("This is 7 \n");
   int iter = 5;
   for(int i=0;i<iter;i++){
      Central_diff();
      printf("This is 8 \n");
      Steepest_descent();
      printf("This is 9 \n");
      Forward_track();
      printf("Done with iteration: %d\n",i);
   }
   printf("This is 10 \n");
 
   WriteMultiFilaments(); 

   //SingleFilField();
   //WriteBoundaryNC();
   //MultiFilFieldSym();
   
   //Check to see if this takes majority of time 
   //WriteSingleB();
   //WriteMultiB();
   

 
   
   
   
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
      for(int i=0;i<Ncoils;i++){
      //   printf("The current of the %d th coil is:   %f\n",i, currents[i]);
      }
      printf("The centroid of the first coil is:   %f   %f   %f\n", cx[0],cy[0],cz[0]);
   
      WriteSingleFilaments();
      WriteBoundary();
    //  WriteBoundaryNC();
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

end = clock();
tot_time = ((double) (end-start)) / CLOCKS_PER_SEC;
printf("\nTotal time taken is: %f\n", tot_time); 
}
