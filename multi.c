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
   int i;
   double tot_time;
   clock_t start, end;
   // At some point use higher resolution timer 
   
   start = clock(); 
   //printf("This is 1 \n");
   SetInputs();
   //printf("This is 2 \n");
   ReadFocusInts(focus_output);

   printf("\nThis is OMIC...\n\n");
   printf("The number of iterations is: %d\n",niter);
   printf("The number of alpha harmonics is: %d\n",NFalpha);
   printf("The number of segments is: %d\n", Nseg);
   printf("The dimensions (radxtor) are : %.4f x %.4f \n",len_rad,len_tor);
   printf("The filaments (radxtor) are : %d x %d\n",Nradfil, Ntorfil);
   printf("The spectral weighting factor is: %.6f \n", nvals_scaling);
   printf("The complexity weighting is: %.6f \n", weight_comp);  


//   printf("This is 3 \n");
   ReadFocusArrays(focus_output);
//   printf("This is 4 \n");
   UnpackSingleFilaments();
//   printf("This is 5 \n");
   Init_alpha(case_alpha);
//   printf("This is 6 \n");
   CalculateMultiFilaments();
   
   MultiFilFieldSym();
   //printf("This is 7 \n");

   double sfil_error = SingleFieldError();
   printf("The single-filament error is %.15f\n", sfil_error);
    
   //printf("This is 8 \n");  
   double multi_error_init = MultiFieldError();
   printf("The initial multi-filament error is %.15f\n", multi_error_init);

   
   double comp_penalty_init = ComplexityPenalty();
   printf("The initial complexity error is %.15f\n", comp_penalty_init);

   for(int i=0;i<niter;i++){
      Central_diff(alpamps,multi_error_init);
      //printf("This is 9 \n");
      Steepest_descent();
      //printf("This is 10 \n");
      Forward_track(multi_error_init);
      printf("Done with iteration: %d\n",i+1);
      WriteOutputNC();
   }

   surface_area = SurfaceArea();

   printf("The surface area is %.15f\n", surface_area);

   WriteMultiFilaments(); 
  
   WriteOutputNC(); 
     
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

   end = clock();
   tot_time = ((double) (end-start)) / CLOCKS_PER_SEC;
   printf("\nTotal time taken is: %f\n", tot_time); 

}

