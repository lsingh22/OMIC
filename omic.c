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
#include "startup.h"

int pn, nproc;

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

int main(int argc, char **argv){
   
   printf("\n*********************************************************"); 
   printf("\n         This is the stellarator coil tool OMIC.\n          Authors: Luquant Singh, Thomas Kruger");
   printf("\n*********************************************************\n"); 
 
   if(argc==1){printf("\nERROR: File extension for 'prefix.input' not specified.\n\n"); exit(0);}
   char *ext = *(argv+1);

   OMICStartup(ext); 
  
   double tot_time;
   double multi_error_init, comp_penalty_init;
   double start, end;
      
   //Initialize MPI, store node configuration and current node
   MPI_Init(&argc, &argv);   
   MPI_Comm_rank(MPI_COMM_WORLD, &pn);
   MPI_Comm_size(MPI_COMM_WORLD, &nproc);   
   start = MPI_Wtime();
   
//   SetInputs();

   ReadInputs();

   ReadFocusInts(focus_output);
   
   MPISetup();

   if(pn==0)
   {
      printf("\n                      USER INPUTS               ");
      printf("\n*********************************************************\n\n"); 
      printf("The number of iterations is:                           %d\n",niter);
      printf("The number of alpha harmonics is:                      %d\n",NFalpha);
      printf("The number of segments is:                            %d\n", Nseg);
      printf("The dimensions are:                        %.1fcm x %.1fcm \n",len_rad*100.0,len_tor*100.0);
      printf("The filaments are:                                  %d x %d\n",Nradfil, Ntorfil);
      printf("The spectral weighting factor is:                %.6f \n", nvals_scaling);
      printf("The complexity weighting is:                     %.6f \n\n", weight_comp);  
      printf("*********************************************************\n"); 
   }

   ReadFocusArrays(focus_output);
   
   Initialize();
   
   UnpackSingleFilaments();
   
//   SingleFilField();
  
   Init_alpha(case_alpha);
   
   CalculateMultiFilaments();
   
   MultifilamentField();

   if(nproc > 1){GatherFieldData();}
      
   multi_error_init = MultiFieldError();

   if(pn==0)
   {
      printf("\n                      INITIAL VALUES             ");
      printf("\n*********************************************************\n\n"); 
      multi_error_init = MultiFieldError();
      printf("Multi-filament error is:                %.15f\n", multi_error_init);
      comp_penalty_init = ComplexityPenalty();
      printf("Complexity error is:                    %.15f\n", comp_penalty_init);
      printf("\n*********************************************************\n\n"); 
      surface_area = SurfaceArea();
      printf("The surface area is %.15f\n", surface_area);
   }
   
//TODO: The following should be calling something like OptimizeCoils on each iteration.
   
   for(int i=0;i<niter;i++)
   {
      if(pn==0){printf("\n********************** Iteration %d **********************\n",i);} 
      if(pn==0){printf("\nCalculating derivatives...");}
   
      Central_diff(alpamps,multi_error_init);
      if(pn==0){printf("Done.\n\n");}
  
      Steepest_descent();
    
      Forward_track(multi_error_init);      
      if(pn==0){printf("Done with iteration: %d\n",i+1);}
   }
 
   if(pn==0)
   {
      WriteMultiFilaments(); 
      WriteOutputNC(); 
   }  

   end = MPI_Wtime();
   tot_time = end - start;
   if(pn==0){printf("\nTotal time taken is: %f\n\n", tot_time);} 

   MPI_Finalize();
}

//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----//----

