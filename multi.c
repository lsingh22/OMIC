#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "read_namelist.h"  
#include "read_focus.h"
#include "single_fil.h"
#include "multi_fil.h"
#include "bfield.h"
#include "alpha.h"
#include "output.h"
#include "solvers.h"

int pn, nproc;

int main(int argc, char **argv){
   int i, ierr;
   double tot_time;
   clock_t start, end;
  // int pn, nproc;

   //Initialize and prepare MPI stuff
   MPI_Init(&argc, &argv);   
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pn);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);   

   start = clock();

   if(pn==0){printf("This is 0 \n");}
   SetInputs();

   if(pn==0){printf("This is 1 \n");}
   ReadFocusInts(focus_output);

   if(pn==0){printf("This is 2 \n");}      
   MPInit();

   if(pn==0)
   {
      printf("\nThis is the stellarator coil tool OMIC...\n\n");
      printf("The number of iterations is: %d\n",niter);
      printf("The number of alpha harmonics is: %d\n",NFalpha);
      printf("The number of segments is: %d\n", Nseg);
      printf("The dimensions (radxtor) are : %.4f x %.4f \n",len_rad,len_tor);
      printf("The filaments (radxtor) are : %d x %d\n",Nradfil, Ntorfil);
      printf("The spectral weighting factor is: %.6f \n", nvals_scaling);
      printf("The complexity weighting is: %.6f \n", weight_comp);  
   }

   if(pn==0){printf("This is 3 \n");}
   ReadFocusArrays(focus_output);
   
   if(pn==0){printf("This is 4 \n");}
   UnpackSingleFilaments();
   
   if(pn==0){printf("This is 5 \n");}
   Init_alpha(case_alpha);
   
   if(pn==0){printf("This is 6 \n");}
   CalculateMultiFilaments();
   
   if(pn==0){printf("This is 7 \n");}
   MultiFilFieldSym();

   if(pn==0){printf("This is 8 \n");}
   if(nproc > 1){GatherFieldData();}

   if(pn==0){printf("This is 9 \n");}
   
   double sfil_error = SingleFieldError();
   if(pn==0){printf("The single-filament error is %.15f\n", sfil_error);} 
   MPI_Bcast(&sfil_error, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   double multi_error_init = MultiFieldError();
   if(pn==0){printf("The initial multi-filament error is %.15f\n", multi_error_init);}
   MPI_Bcast(&multi_error_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
   double comp_penalty_init = ComplexityPenalty();
   if(pn==0){printf("The initial complexity error is %.15f\n", comp_penalty_init);}
   MPI_Bcast(&comp_penalty_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   //TODO: This should really just be calling something like OptimizeCoils or something to match the rest of the file...

   if(pn==0){printf("This is 10 \n");}
   for(int i=0;i<niter;i++)
   {
      if(pn==0){printf("This is 11 \n");}
      Central_diff(alpamps,multi_error_init);
      
      if(pn==0){printf("This is 12 \n");}
      Steepest_descent();
    
      if(pn==0){printf("This is 13 \n");}
      Forward_track(multi_error_init);
         
      printf("Done with iteration: %d\n",i+1);
      //WriteOutputNC();
   }
 

   surface_area = SurfaceArea();
   if(pn==0){printf("The surface area is %.15f\n", surface_area);}

   if(pn==0)
   {
      WriteMultiFilaments(); 
      WriteOutputNC(); 
   }  

   end = clock();
   tot_time = ((double) (end-start)) / CLOCKS_PER_SEC;
   if(pn==0){printf("\nTotal time taken is: %f\n", tot_time);} 

   ierr = MPI_Finalize();
}

