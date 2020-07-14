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
   double sfil_error, multi_error_init, comp_penalty_init;
   double start, end;
   // int pn, nproc;

   //Testbed for speed improvements
   double sum=0;
   double ok = 5;
   double t1,t2;
   int y;
/*
   t1 = MPI_Wtime();
   
   t2 = MPI_Wtime();

   if(pn==0){printf("\nTotal time for pow is: %f\n\n", t2-t1);} 

   t1 = MPI_Wtime();
   sum = pow(ok,ok);
   t2 = MPI_Wtime();

   if(pn==0){printf("\nTotal time for explicit is: %f\n\n", 1000*(t2-t1));} 
*/

   //Initialize and prepare MPI stuff
   MPI_Init(&argc, &argv);   
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &pn);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &nproc);   

   start = MPI_Wtime();

//   if(pn==0){printf("This is 0 \n");}
   SetInputs();

//   if(pn==0){printf("This is 1 \n");}
   ReadFocusInts(focus_output);

//   if(pn==0){printf("This is 2 \n");}      
   MPInit();

   if(pn==0)
   {
      printf("\n*********************************************************"); 
      printf("\n         This is the stellarator coil tool OMIC.\n          Authors: Luquant Singh, Thomas Kruger");
      printf("\n*********************************************************\n"); 
      printf("\n*********************************************************"); 
      printf("\n                      USER INPUTS               ");
      printf("\n*********************************************************\n\n"); 
      printf("The number of iterations is:                           %d\n",niter);
      printf("The number of alpha harmonics is:                      %d\n",NFalpha);
      printf("The number of segments is:                            %d\n", Nseg);
      printf("The dimensions are:                        %.1fcm x %.1fcm \n",len_rad*100,len_tor*100);
      printf("The filaments are:                                  %d x %d\n",Nradfil, Ntorfil);
      printf("The spectral weighting factor is:                %.6f \n", nvals_scaling);
      printf("The complexity weighting is:                     %.6f \n\n", weight_comp);  
      printf("*********************************************************\n"); 
   }

//   if(pn==0){printf("This is 3 \n");}
   ReadFocusArrays(focus_output);
   
//   if(pn==0){printf("This is 4 \n");}
   UnpackSingleFilaments();
   
//   SingleFilField();
  
//   if(pn==0){printf("This is 5 \n");}
   Init_alpha(case_alpha);
   
//   if(pn==0){printf("This is 6 \n");}
   CalculateMultiFilaments();
   
//   if(pn==0){printf("This is 7 \n");}
   MultiFilFieldSym();

//   if(pn==0){printf("This is 8 \n");}
   if(nproc > 1){GatherFieldData();}

//   if(pn==0){printf("This is 9 \n");}

   if(pn==0)
   {
      printf("\n*********************************************************"); 
      printf("\n                      INITIAL VALUES             ");
      printf("\n*********************************************************\n\n"); 
//      sfil_error = SingleFieldError();
//      printf("Single-filament error is:               %.15f\n", sfil_error);
      multi_error_init = MultiFieldError();
      printf("Multi-filament error is:                %.15f\n", multi_error_init);
      comp_penalty_init = ComplexityPenalty();
      printf("Complexity error is:                    %.15f\n", comp_penalty_init);
      printf("\n*********************************************************\n\n"); 
      surface_area = SurfaceArea();
      printf("\nThe surface area is %.15f\n", surface_area);
   }
   
//   MPI_Bcast(&multi_error_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   

//TODO: This should really just be calling something like OptimizeCoils or something to match the rest of the file...

//   if(pn==0){printf("This is 10 \n");} 
   
   for(int i=0;i<niter;i++)
   {
      if(pn==0){printf("\n********************** Iteration %d **********************\n",i);} 
//      if(pn==0){printf("This is 11 \n");}
      if(pn==0){printf("\nCalculating derivatives...");}
      Central_diff(alpamps,multi_error_init);
      if(pn==0){printf("Done.\n\n");}
//      if(pn==0){printf("This is 12 \n");}
      Steepest_descent();
    
//      if(pn==0){printf("This is 13 \n");}
      Forward_track(multi_error_init);
         
      if(pn==0){printf("Done with iteration: %d\n",i+1);}
      //WriteOutputNC();
   }
 
   if(pn==0)
   {
      WriteMultiFilaments(); 
      WriteOutputNC(); 
   }  

//   if(pn==0){printf("This is end \n");} 

   end = MPI_Wtime();
   tot_time = ((double) (end-start)) / CLOCKS_PER_SEC;
   if(pn==0){printf("\nTotal time taken is: %f\n\n", tot_time);} 

   ierr = MPI_Finalize();
}

